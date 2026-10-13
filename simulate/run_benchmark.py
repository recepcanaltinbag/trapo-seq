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
import json
import os
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


def fastq_to_fasta_gnu_sed(fq, fa):
    """Exact equivalent of GNU `sed -n '1~4s/^@/>/p;2~4p' fq > fa`."""
    with open(fq) as fin, open(fa, "w") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                if line.startswith("@"):
                    fout.write(">" + line[1:])
            elif i % 4 == 1:
                fout.write(line)


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


def run(cmd, cwd, log, env):
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(f"$ {cmd}\n# cwd: {cwd}\n\n")
        fh.flush()
        p = subprocess.run(cmd, shell=True, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    return {"cmd": cmd, "returncode": p.returncode, "seconds": round(time.time() - t0, 1), "log": log}


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
    ap.add_argument("--work_dir", required=True, help="Scratch directory for the run (wiped first)")
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
    args = ap.parse_args(argv)

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

    data, local, barcodes = prepare(sim_dir, work, refs, rename=not args.keep_read_ids)
    logs = os.path.join(work, "logs")
    os.makedirs(logs)
    os.makedirs(out_dir, exist_ok=True)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1")
    shim_b52 = args.shim_b52 == "yes" or (args.shim_b52 == "auto" and args.profile == "v0.1.0")
    if shim_b52:
        shim = os.path.join(work, "shim")
        os.makedirs(shim)
        with open(os.path.join(shim, "sitecustomize.py"), "w") as fh:
            fh.write(SHIM_B52)
        env["PYTHONPATH"] = shim + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    q = shlex.quote
    fields = {"py": q(args.python), "cli": q(os.path.join(checkout, "trapo-seq.py")), "data": q(data),
              "plasmid": q(local["plasmid"]), "genome": q(local["genome"]), "is_fasta": q(local["is_fasta"]),
              "temp": q(os.path.join(work, "blast_temp")), "threads": args.threads, "min_len": args.min_len}

    tools = {}
    for t in ("minimap2", "samtools", "blastn", "makeblastdb", "mafft", "sed"):
        tools[t] = shutil.which(t, path=env.get("PATH"))
    git = subprocess.run(["git", "-C", checkout, "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    record = {"checkout": checkout, "commit": git.stdout.strip() or None, "profile": args.profile,
              "sim_dir": sim_dir, "work_dir": work, "refs": refs, "barcodes": barcodes,
              "read_ids_rewritten": not args.keep_read_ids, "fix_b01": fix_b01, "shim_b52": shim_b52, "tools": tools,
              "python": args.python, "steps": [], "workarounds": []}

    for step in [s for s in STEP_ORDER if s in args.steps]:
        tmpl = steps[step]
        if step == "filter":
            for bc in barcodes:
                f = dict(fields, barcode=bc, fastq_in=q(os.path.join(data, bc, f"{bc}_raw.fastq")),
                         fastq_out=q(os.path.join(data, bc, f"{bc}.fastq")))
                r = run(tmpl.format(**f), checkout, os.path.join(logs, f"filter_{bc}.log"), env)
                r["step"] = f"filter:{bc}"
                record["steps"].append(r)
                # v0.1.0 map looks for exactly <bc>/<bc>.fastq; the raw file must not shadow it
            continue
        r = run(tmpl.format(**fields), checkout, os.path.join(logs, f"{step}.log"), env)
        r["step"] = step
        r["traceback_in_log"] = log_has_traceback(r["log"])
        record["steps"].append(r)
        print(f"[{step}] rc={r['returncode']} traceback={r['traceback_in_log']} {r['seconds']}s", flush=True)
        if step == "map" and fix_b01:
            for bc in barcodes:
                fa = os.path.join(data, bc, f"{bc}.fasta")
                fq = os.path.join(data, bc, "04_mapped_to_plasmid.fastq")
                if os.path.exists(fq) and (not os.path.exists(fa) or os.path.getsize(fa) == 0):
                    fastq_to_fasta_gnu_sed(fq, fa)
                    record["workarounds"].append(
                        f"B01: regenerated {bc}/{bc}.fasta from 04_mapped_to_plasmid.fastq "
                        "(Python equivalent of GNU sed -n '1~4s/^@/>/p;2~4p')")
            # B07: 08_mapping_stats.txt is written after the sed line; it is informational only.

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
