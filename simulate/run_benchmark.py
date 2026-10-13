#!/usr/bin/env python3
"""
Run a trapo-seq checkout on a simulated dataset and score it with evaluate.py.

    python simulate/run_benchmark.py --checkout PATH_TO_TRAPO_SEQ \\
        --sim_dir SIM_DIR --work_dir WORK_DIR --out_dir EVAL_DIR [--profile v0.1.0|dev]

What it does
------------
1. Builds WORK_DIR/data/<barcode>/<barcode>_raw.fastq from SIM_DIR (read IDs
   optionally rewritten '_' -> '-', see --keep_read_ids) and copies the plasmid,
   genome and IS FASTA into WORK_DIR/ref/ (fresh BLAST DBs every run, AUDIT B17).
2. Runs the pipeline steps through the checkout's trapo-seq.py, with the
   checkout as working directory (v0.1.0 resolves helper scripts relative to
   the cwd, AUDIT B50):
       filter (per barcode) -> map -> insert_finder_batch -> blast_annot_batch -> dr_finder
   Each step's stdout/stderr goes to WORK_DIR/logs/<step>.log; a failing step is
   recorded in run_log.json and the remaining steps still run, so that partial
   output is scored (missing outputs become N/A in the evaluation).
3. Runs simulate/evaluate.py (from this repository) on WORK_DIR/data.

Each step record in run_log.json has "seconds" (wall time) and "cpu_seconds"
(user + system CPU time of the step's child processes, from
resource.getrusage(RUSAGE_CHILDREN) deltas; it includes every descendant the
shell waited for, e.g. minimap2 | samtools, and counts all their threads).

Reusing upstream work (--reuse_work, --from_step, --snapshot_after)
--------------------------------------------------------------------
Every run writes WORK_DIR/.run_benchmark_stamp.json: the identity of the run
(sim dir and a fingerprint of its truth.tsv/clones.tsv/simulation.json,
checkout path, commit and a fingerprint of uncommitted code changes, profile,
step templates, reference FASTA hashes, read-ID/B01/B52 options) and, after
each step, the step's records and a manifest of the work dir (path, size,
mtime). With --reuse_work the work dir is kept, and with --from_step STEP only
STEP and the steps after it run (--from_step implies --reuse_work; --reuse_work
alone resumes after the last completed step). Before anything runs the runner
checks that
  * the stamp exists and its identity matches the current invocation (the
    templates are compared only for the steps that are reused),
  * every reused step completed, and the work dir still matches the manifest
    written after the last reused step (no file missing or changed);
files that later steps created are deleted, so the next step sees exactly the
state it would see in a run from scratch. On a mismatch the runner exits with
status 3 (REUSE_MISMATCH_RC) and changes nothing. The reused steps' records
(wall and CPU time) are copied into run_log.json with "reused": true, so
runtime totals still describe a full pipeline run.
--snapshot_after STEP DIR copies the work dir (stamp included) to DIR right
after STEP; run_grid.py uses it to run map + insert_finder once and every
annotation variant from blast_annot_batch on a copy.

Step commands are plain templates in PROFILES below; adjust them there (or pass
--steps_json FILE with a {"step": "template"} mapping) when the CLI changes.
Template fields: {py} {cli} {data} {plasmid} {genome} {is_fasta} {temp}
{threads} {min_len} {fastq_in} {fastq_out} {barcode}.

Notes for v0.1.0
----------------
* The map step writes <barcode>.fasta with GNU-only `sed -n '1~4s/^@/>/p;2~4p'`
  (AUDIT B01). With --fix_b01 (default for profile v0.1.0) the runner
  regenerates that file after mapping with a byte-for-byte Python equivalent
  of the GNU sed command whenever it is missing or empty. It is skipped when
  GNU sed is the `sed` on PATH and produced a non-empty file.
* v0.1.0 imports Bio.Blast.Applications (AUDIT B52), which Biopython >= 1.86
  removed. The import is unused (BLAST runs through os.system), so with
  --shim_b52 (default for profile v0.1.0) the runner puts a sitecustomize.py
  on PYTHONPATH that registers an empty stand-in module. No pipeline code
  is changed.
* v0.1.0 builds commands with os.system and unquoted paths (AUDIT B18), so the
  work directory must not contain spaces; the runner refuses such paths for
  that profile.
"""

import argparse
import hashlib
import json
import os
import resource
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

PROFILES = {
    # The original release (main @ 980d5e6).
    "v0.1.0": {
        "filter": "{py} {cli} filter -f {fastq_in} -o {fastq_out} -l {min_len}",
        "map": "{py} {cli} map -d {data} -p {plasmid} -g {genome} --force",
        "insert_finder_batch": "{py} {cli} insert_finder_batch -d {data}",
        "blast_annot_batch": "{py} {cli} blast_annot_batch -d {data} -g {genome} --is_fasta {is_fasta} --temp {temp} --threads {threads}",
        "dr_finder": "{py} {cli} dr_finder -d {data} -p {plasmid}",
    },
}
# The development branch starts from the same CLI; edit here as it diverges.
PROFILES["dev"] = dict(PROFILES["v0.1.0"])
PROFILES["dev"]["map"] = "{py} {cli} map -d {data} -p {plasmid} -g {genome} --force --threads {threads}"
PROFILES["dev"]["dr_finder"] = "{py} {cli} dr_finder -d {data} -p {plasmid} --is_fasta {is_fasta} --threads {threads}"
SHIM_B52 = """# Written by simulate/run_benchmark.py (AUDIT B52): v0.1.0 imports the unused
# Bio.Blast.Applications.NcbiblastnCommandline, removed in Biopython >= 1.86.
import sys, types
try:
    import Bio.Blast.Applications  # noqa: F401
except Exception:
    try:
        import Bio.Blast
        _m = types.ModuleType("Bio.Blast.Applications")
        class NcbiblastnCommandline:  # never called by v0.1.0
            def __init__(self, *a, **k):
                raise RuntimeError("stand-in only")
        _m.NcbiblastnCommandline = NcbiblastnCommandline
        sys.modules["Bio.Blast.Applications"] = _m
        Bio.Blast.Applications = _m
    except Exception:
        pass
"""
STEP_ORDER = ["filter", "map", "insert_finder_batch", "blast_annot_batch", "dr_finder"]
STAMP = ".run_benchmark_stamp.json"
STAMP_VERSION = 1
REUSE_MISMATCH_RC = 3
CODE_EXT = (".py", ".sh", ".pl", ".R", ".r")


def fastq_to_fasta_gnu_sed(fq, fa):
    """Exact equivalent of GNU `sed -n '1~4s/^@/>/p;2~4p' fq > fa`."""
    with open(fq) as fin, open(fa, "w") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                if line.startswith("@"):
                    fout.write(">" + line[1:])
            elif i % 4 == 1:
                fout.write(line)


def sha1_file(path, h=None):
    h = h or hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h


def sim_fingerprint(sim_dir):
    """Hash of the simulation's truth files (cheap; the FASTQ follows from them,
    and run_grid may already have removed it)."""
    h = hashlib.sha1()
    for name in ("simulation.json", "truth.tsv", "clones.tsv"):
        p = os.path.join(sim_dir, name)
        h.update(name.encode())
        if os.path.exists(p):
            sha1_file(p, h)
    return h.hexdigest()[:16]


def code_fingerprint(checkout):
    """(commit, hash of uncommitted code changes or ''). Uncommitted changes are
    the tracked diff against HEAD plus untracked code files; compiled files
    are ignored (tracked in v0.1.0, they change whenever it runs)."""
    def g(*a):
        p = subprocess.run(["git", "-C", checkout, *a], capture_output=True)
        return p.stdout if p.returncode == 0 else b""
    commit = g("rev-parse", "HEAD").decode().strip() or None
    excl = ["--", ".", ":(exclude)*.pyc", ":(exclude)*/__pycache__/*"]
    h = hashlib.sha1(g("diff", "HEAD", *excl))
    dirty = bool(g("status", "--porcelain", "--untracked-files=no", *excl).strip())
    for rel in sorted(g("ls-files", "--others", "--exclude-standard").decode().splitlines()):
        if rel.endswith(CODE_EXT) and os.path.isfile(os.path.join(checkout, rel)):
            dirty = True
            h.update(rel.encode())
            sha1_file(os.path.join(checkout, rel), h)
    return commit, h.hexdigest()[:16] if dirty else ""


def manifest(work):
    """{relative path: [size, mtime_ns]} of every file in the work dir except
    the stamp, and {relative dir + '/': ['dir']} for every sub-directory."""
    out = {}
    for root, dirs, files in os.walk(work):
        dirs.sort()
        for d in dirs:
            out[os.path.relpath(os.path.join(root, d), work) + "/"] = ["dir"]
        for f in sorted(files):
            p = os.path.join(root, f)
            rel = os.path.relpath(p, work)
            if rel in (STAMP, STAMP + ".tmp"):
                continue
            st = os.lstat(p)
            out[rel] = [st.st_size, st.st_mtime_ns]
    return out


def write_stamp(work, stamp):
    tmp = os.path.join(work, STAMP + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(stamp, fh, indent=1)
    os.replace(tmp, os.path.join(work, STAMP))


def copy_work(src, dst):
    """Copy a work dir preserving mtimes (manifest check); APFS clone when possible."""
    if os.path.exists(dst):
        shutil.rmtree(dst)
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    if sys.platform == "darwin":
        if subprocess.run(["cp", "-cRp", src, dst], capture_output=True).returncode == 0:
            return dst
        shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst, copy_function=shutil.copy2, symlinks=True)
    return dst


class ReuseError(Exception):
    pass


def check_reuse(work, identity, templates, from_step, wanted=STEP_ORDER):
    """Validate WORK_DIR for a run that starts at `from_step`.

    Returns (stamp, reused step names). Raises ReuseError on any mismatch;
    on success deletes the files created after the last reused step."""
    path = os.path.join(work, STAMP)
    if not os.path.exists(path):
        raise ReuseError(f"no {STAMP} in {work}: nothing to reuse")
    with open(path) as fh:
        stamp = json.load(fh)
    if stamp.get("version") != STAMP_VERSION:
        raise ReuseError(f"stamp version {stamp.get('version')} != {STAMP_VERSION}")
    old = stamp.get("identity", {})
    diff = sorted(k for k in set(old) | set(identity) if old.get(k) != identity.get(k))
    if diff:
        raise ReuseError("work dir was produced by a different run: " +
                         "; ".join(f"{k}: {old.get(k)!r} != {identity.get(k)!r}" for k in diff))
    done = stamp.get("done", [])
    if from_step is None:  # resume after the last completed step
        later = [s for s in STEP_ORDER if s not in done]
        from_step = later[0] if later else None
    cut = STEP_ORDER.index(from_step) if from_step else len(STEP_ORDER)
    reused = [s for s in done if STEP_ORDER.index(s) < cut]
    missing = [s for s in STEP_ORDER[:cut] if s in wanted and s not in reused]
    if missing:
        raise ReuseError(f"upstream step(s) {missing} never completed in this work dir")
    for s in reused:
        if stamp.get("templates", {}).get(s) != templates.get(s):
            raise ReuseError(f"step '{s}' ran with a different command template")
    state = reused[-1] if reused else "prepare"
    want = stamp.get("manifests", {}).get(state)
    if want is None:
        raise ReuseError(f"no manifest for state '{state}' in the stamp")
    have = manifest(work)
    bad = sorted(rel for rel, v in want.items() if have.get(rel) != v)
    if bad:
        raise ReuseError(f"{len(bad)} file(s) from step '{state}' missing or changed since, e.g. {bad[:3]}")
    extra = set(have) - set(want)
    for rel in sorted(r for r in extra if not r.endswith("/")):
        os.remove(os.path.join(work, rel))
    for rel in sorted((r for r in extra if r.endswith("/")), reverse=True):
        shutil.rmtree(os.path.join(work, rel), ignore_errors=True)
    stamp["done"] = reused
    stamp["manifests"] = {k: v for k, v in stamp.get("manifests", {}).items() if k in reused + ["prepare"]}
    stamp["records"] = {k: v for k, v in stamp.get("records", {}).items() if k in reused}
    return stamp, reused, from_step


def prepare(sim_dir, work, refs, rename):
    data = os.path.join(work, "data")
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(data)
    ref_dir = os.path.join(work, "ref")
    os.makedirs(ref_dir)
    local = {}
    for key, src in refs.items():
        dst = os.path.join(ref_dir, os.path.basename(src))
        shutil.copy(src, dst)
        local[key] = dst
    barcodes = sorted(d for d in os.listdir(sim_dir)
                      if os.path.isdir(os.path.join(sim_dir, d)) and os.path.exists(os.path.join(sim_dir, d, f"{d}_raw.fastq")))
    for bc in barcodes:
        os.makedirs(os.path.join(data, bc))
        src = os.path.join(sim_dir, bc, f"{bc}_raw.fastq")
        dst = os.path.join(data, bc, f"{bc}_raw.fastq")
        with open(src) as fin, open(dst, "w") as fout:
            for i, line in enumerate(fin):
                if rename and i % 4 == 0:
                    rid, sep, rest = line[1:].rstrip("\n").partition(" ")
                    line = "@" + rid.replace("_", "-") + (sep + rest if sep else "") + "\n"
                fout.write(line)
    return data, local, barcodes


def children_cpu():
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_utime + ru.ru_stime


def run(cmd, cwd, log, env):
    t0 = time.time()
    c0 = children_cpu()
    with open(log, "w") as fh:
        fh.write(f"$ {cmd}\n# cwd: {cwd}\n\n")
        fh.flush()
        p = subprocess.run(cmd, shell=True, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    return {"cmd": cmd, "returncode": p.returncode, "seconds": round(time.time() - t0, 1),
            "cpu_seconds": round(children_cpu() - c0, 2), "log": log}


def log_has_traceback(path):
    try:
        with open(path, errors="replace") as fh:
            return "Traceback (most recent call last)" in fh.read()
    except OSError:
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkout", required=True, help="trapo-seq checkout containing trapo-seq.py")
    ap.add_argument("--sim_dir", required=True, help="Output of simulate_traposome.py")
    ap.add_argument("--work_dir", required=True, help="Scratch directory for the run (wiped first unless --reuse_work/--from_step)")
    ap.add_argument("--out_dir", required=True, help="Where evaluation results and run_log.json go")
    ap.add_argument("--profile", default="dev", choices=sorted(PROFILES))
    ap.add_argument("--steps_json", help="JSON file overriding step command templates")
    ap.add_argument("--steps", nargs="*", default=STEP_ORDER, help="Subset of steps to run")
    ap.add_argument("-p", "--plasmid", default=None, help="Plasmid FASTA (default: from simulation.json)")
    ap.add_argument("-g", "--genome", default=None)
    ap.add_argument("-i", "--is_fasta", default=None)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--min_len", type=int, default=2000, help="filter -l (default 2000, as in test_script.sh)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tol", type=int, default=20)
    ap.add_argument("--keep_read_ids", action="store_true",
                    help="Do not rewrite '_' in read IDs (v0.1.0 dr_finder then fails, AUDIT B34)")
    ap.add_argument("--fix_b01", choices=["auto", "yes", "no"], default="auto",
                    help="Regenerate <barcode>.fasta after map if empty (auto: yes for v0.1.0)")
    ap.add_argument("--shim_b52", choices=["auto", "yes", "no"], default="auto",
                    help="Provide a stand-in Bio.Blast.Applications module (auto: yes for v0.1.0)")
    ap.add_argument("--reuse_work", action="store_true",
                    help="Keep the existing work dir (checked against its stamp) and resume after the "
                         "last completed step, or at --from_step")
    ap.add_argument("--from_step", choices=STEP_ORDER, default=None,
                    help="Run only this step and the later ones on the reused work dir (implies --reuse_work)")
    ap.add_argument("--snapshot_after", nargs=2, metavar=("STEP", "DIR"), default=None,
                    help="Copy the work dir (with its stamp) to DIR right after STEP finished")
    args = ap.parse_args(argv)
    if args.snapshot_after and args.snapshot_after[0] not in STEP_ORDER:
        ap.error(f"--snapshot_after STEP must be one of {STEP_ORDER}")
    reuse = args.reuse_work or args.from_step is not None

    checkout = os.path.abspath(args.checkout)
    sim_dir = os.path.abspath(args.sim_dir)
    work = os.path.abspath(args.work_dir)
    out_dir = os.path.abspath(args.out_dir)
    if args.profile == "v0.1.0" and " " in work:
        sys.exit("v0.1.0 cannot handle spaces in paths (AUDIT B18); choose a --work_dir without spaces")
    sim_meta = {}
    if os.path.exists(os.path.join(sim_dir, "simulation.json")):
        with open(os.path.join(sim_dir, "simulation.json")) as fh:
            sim_meta = json.load(fh)

    def resolve_ref(given, key):
        path = given or sim_meta.get(key)
        if path and not os.path.isabs(path) and not os.path.exists(path):
            path = os.path.join(os.path.dirname(HERE), path)  # simulation.json paths are relative to the repo root
        if not path or not os.path.exists(path):
            sys.exit(f"Cannot find {key} FASTA ({path}); pass it explicitly")
        return os.path.abspath(path)

    refs = {"plasmid": resolve_ref(args.plasmid, "plasmid"), "genome": resolve_ref(args.genome, "genome"),
            "is_fasta": resolve_ref(args.is_fasta, "is_fasta")}
    steps = dict(PROFILES[args.profile])
    if args.steps_json:
        with open(args.steps_json) as fh:
            steps.update(json.load(fh))
    fix_b01 = args.fix_b01 == "yes" or (args.fix_b01 == "auto" and args.profile == "v0.1.0")
    shim_b52 = args.shim_b52 == "yes" or (args.shim_b52 == "auto" and args.profile == "v0.1.0")
    commit, code_fp = code_fingerprint(checkout)
    identity = {"sim_dir": sim_dir, "sim_fingerprint": sim_fingerprint(sim_dir), "checkout": checkout,
                "commit": commit, "uncommitted_code": code_fp, "profile": args.profile,
                "python": args.python, "min_len": args.min_len, "read_ids_rewritten": not args.keep_read_ids,
                "fix_b01": fix_b01, "shim_b52": shim_b52,
                "refs": {k: [os.path.basename(v), sha1_file(v).hexdigest()[:16]] for k, v in sorted(refs.items())}}
    templates = {s: steps[s] for s in STEP_ORDER}

    reused, from_step = [], args.from_step
    if reuse:
        try:
            stamp, reused, from_step = check_reuse(work, identity, templates, args.from_step, args.steps)
        except ReuseError as exc:
            print(f"--reuse_work refused: {exc}", file=sys.stderr, flush=True)
            sys.exit(REUSE_MISMATCH_RC)
        data = os.path.join(work, "data")
        local = {k: os.path.join(work, "ref", os.path.basename(v)) for k, v in refs.items()}
        barcodes = stamp["barcodes"]
    else:
        data, local, barcodes = prepare(sim_dir, work, refs, rename=not args.keep_read_ids)
        stamp = {"version": STAMP_VERSION, "identity": identity, "templates": templates, "barcodes": barcodes,
                 "done": [], "records": {}, "manifests": {}}
    stamp["templates"] = dict(stamp.get("templates", {}), **{s: templates[s] for s in STEP_ORDER
                                                             if s not in reused})
    logs = os.path.join(work, "logs")
    os.makedirs(logs, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1")
    if shim_b52:
        shim = os.path.join(work, "shim")
        if not reuse:
            os.makedirs(shim)
            with open(os.path.join(shim, "sitecustomize.py"), "w") as fh:
                fh.write(SHIM_B52)
        env["PYTHONPATH"] = shim + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    if not reuse:
        stamp["manifests"]["prepare"] = manifest(work)
    write_stamp(work, stamp)
    q = shlex.quote
    fields = {"py": q(args.python), "cli": q(os.path.join(checkout, "trapo-seq.py")), "data": q(data),
              "plasmid": q(local["plasmid"]), "genome": q(local["genome"]), "is_fasta": q(local["is_fasta"]),
              "temp": q(os.path.join(work, "blast_temp")), "threads": args.threads, "min_len": args.min_len}

    tools = {}
    for t in ("minimap2", "samtools", "blastn", "makeblastdb", "mafft", "sed"):
        tools[t] = shutil.which(t, path=env.get("PATH"))
    git = subprocess.run(["git", "-C", checkout, "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    record = {"checkout": checkout, "commit": git.stdout.strip() or None, "uncommitted_code": code_fp,
              "profile": args.profile, "sim_dir": sim_dir, "work_dir": work, "refs": refs, "barcodes": barcodes,
              "read_ids_rewritten": not args.keep_read_ids, "fix_b01": fix_b01, "shim_b52": shim_b52, "tools": tools,
              "python": args.python, "steps": [], "workarounds": [],
              "reused_steps": reused, "from_step": from_step if reuse else None}
    for s in reused:
        for r in stamp["records"].get(s, []):
            record["steps"].append(dict(r, reused=True))
        record["workarounds"] += stamp.get("workarounds", {}).get(s, [])

    cut = STEP_ORDER.index(from_step) if (reuse and from_step) else (len(STEP_ORDER) if reuse else 0)
    to_run = [s for s in STEP_ORDER if s in args.steps and STEP_ORDER.index(s) >= cut]
    for step in to_run:
        tmpl = steps[step]
        recs, notes = [], []
        if step == "filter":
            for bc in barcodes:
                f = dict(fields, barcode=bc, fastq_in=q(os.path.join(data, bc, f"{bc}_raw.fastq")),
                         fastq_out=q(os.path.join(data, bc, f"{bc}.fastq")))
                r = run(tmpl.format(**f), checkout, os.path.join(logs, f"filter_{bc}.log"), env)
                r["step"] = f"filter:{bc}"
                recs.append(r)
                # v0.1.0 map looks for exactly <bc>/<bc>.fastq; the raw file must not shadow it
        else:
            r = run(tmpl.format(**fields), checkout, os.path.join(logs, f"{step}.log"), env)
            r["step"] = step
            r["traceback_in_log"] = log_has_traceback(r["log"])
            recs.append(r)
            print(f"[{step}] rc={r['returncode']} traceback={r['traceback_in_log']} {r['seconds']}s "
                  f"cpu={r['cpu_seconds']}s", flush=True)
        if step == "map" and fix_b01:
            for bc in barcodes:
                fa = os.path.join(data, bc, f"{bc}.fasta")
                fq = os.path.join(data, bc, "04_mapped_to_plasmid.fastq")
                if os.path.exists(fq) and (not os.path.exists(fa) or os.path.getsize(fa) == 0):
                    fastq_to_fasta_gnu_sed(fq, fa)
                    notes.append(
                        f"B01: regenerated {bc}/{bc}.fasta from 04_mapped_to_plasmid.fastq "
                        "(Python equivalent of GNU sed -n '1~4s/^@/>/p;2~4p')")
            # B07: 08_mapping_stats.txt is written after the sed line; it is informational only.
        record["steps"] += recs
        record["workarounds"] += notes
        stamp["done"] = [s for s in stamp["done"] if s != step] + [step]
        stamp["records"][step] = recs
        stamp.setdefault("workarounds", {})[step] = notes
        stamp["manifests"][step] = manifest(work)
        write_stamp(work, stamp)
        if args.snapshot_after and args.snapshot_after[0] == step:
            copy_work(work, os.path.abspath(args.snapshot_after[1]))

    with open(os.path.join(out_dir, "run_log.json"), "w") as fh:
        json.dump(record, fh, indent=2)

    ev = [args.python, os.path.join(HERE, "evaluate.py"), "-s", sim_dir, "-r", data, "-o", out_dir,
          "--tol", str(args.tol), "--min_read_len", str(args.min_len)]
    subprocess.run(ev, check=True)
    # Keep the raw pipeline tables next to the scores for later inspection.
    raw = os.path.join(out_dir, "pipeline_outputs")
    os.makedirs(raw, exist_ok=True)
    for bc in barcodes:
        for name in ("insertion_from_bam.tab", "best_alignment.tab", "08_mapping_stats.txt"):
            src = os.path.join(data, bc, name)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(raw, f"{bc}_{name}"))
    ins = os.path.join(data, "insertions")
    if os.path.isdir(ins):
        for f in os.listdir(ins):
            shutil.copy(os.path.join(ins, f), os.path.join(raw, f))
    shutil.copytree(logs, os.path.join(out_dir, "logs"), dirs_exist_ok=True)


if __name__ == "__main__":
    main()
