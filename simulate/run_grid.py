#!/usr/bin/env python3
"""
Benchmark grid: simulations x pipeline versions x pipeline settings.

    python simulate/run_grid.py --preset full  [--out DIR] [--jobs 2] [--threads 4]
    python simulate/run_grid.py --spec grid.json --dry_run
    python simulate/run_grid.py --out DIR --aggregate_only

For every *cell* (one simulated dataset + one pipeline variant) the driver
  1. simulates the dataset with simulate_traposome.py (once per distinct set of
     simulation parameters; shared by all pipeline variants of that dataset),
  2. runs simulate/run_benchmark.py on it (pipeline + evaluate.py),
  3. deletes the large intermediates (run_benchmark work dir with BAM/FASTQ;
     the simulated FASTQ once every cell of that dataset is done) unless --keep,
and finally aggregates all finished cells into tidy tables and figures.

Grid specification
------------------
A Python dict (PRESETS below) or a JSON file with the same keys:

  name          free text
  refs          {"plasmid", "genome", "is_fasta"}: FASTA paths (relative to the
                repository root or absolute)
  sim_defaults  simulate_traposome.py options shared by every dataset (option
                name without dashes; lists for nargs options; null = not passed).
                read_accuracy_sd may be "scaled:F", meaning F * (1 - mean accuracy),
                so the per-read spread shrinks as reads get better (the r9.4.1
                preset is 0.04 at 0.92 = 0.5 * (1 - 0.92)).
  seeds         replicate seeds (each dataset is simulated once per seed)
  pipelines     {variant_name: {"checkout": "dev" | "baseline" | PATH,
                 "profile": run_benchmark profile, "settings": {option: [values]},
                 "default_settings": {option: value}, "tag": optional label}}
                settings are blast_annot_batch options (evalue, blast_task) and
                are appended to that step's command; a variant without settings
                (v0.1.0) runs with the checkout's defaults.
  sweeps        list of {"name", "set": {fixed sim overrides}, "sim": {option:
                [values]} (factorial within the sweep), "pipelines": "all" |
                "default"}. "default" runs every variant at default_settings only;
                "all" runs the full settings product. Identical cells from
                different sweeps are run once and list every sweep they belong to.
  run           {"min_len": filter length, "tol": evaluation tolerance}

Cell identity and resume
------------------------
sim_id  = sha1 of the canonical JSON of the resolved simulation parameters
          (refs, defaults, overrides, seed), first 12 hex digits.
cell_id = sha1 of {sim parameters, pipeline variant (name, profile, checkout key,
          settings, tag), run options}, first 12 hex digits. Thread counts and
          output paths are not part of the identity. The commit is NOT part of
          it either: it is recorded per cell (pipe.commit, pipe.dirty) and
          aggregation warns when one variant mixes commits. Put a "tag" on the
          variant (or use a fresh --out) to force a re-run after code changes.

Layout under --out (default publication/benchmark/grid/):
  manifest.json             spec, tool versions, checkouts, every cell + status
  sims/<sim_id>/            simulation (truth.tsv, clones.tsv, simulation.json,
                            sim_params.json, sim_stats.json; FASTQ removed later)
  cells/<cell_id>/          run_benchmark output (metrics.json, summary.tsv,
                            per_read.tsv, sites.tsv, run_log.json, logs/,
                            pipeline_outputs/) + cell.json (parameters, commit,
                            tool versions, timings). A cell runs in
                            cells/<cell_id>.partial/ and is renamed when done,
                            so cells/<cell_id>/metrics.json marks a finished
                            cell; such cells are skipped. Failed cells keep
                            cells/<cell_id>.failed/ (and their work dir) and are
                            retried on the next invocation.
  grid_results.tsv          tidy: one row per cell, parameters + every numeric
                            metric of metrics.json["overall"] (dotted names as
                            in summary.tsv) + derived.* + runtime.*
  summary_long.tsv          mean, sd, n over seeds for every metric
  summary_key.tsv           "mean ± sd" for the key paper metrics
  figures/                  plot_grid.py output (PDF + PNG)

Work directories go to --work_root (default: a folder in the system temp dir).
It must not contain spaces when v0.1.0 is in the grid (AUDIT B18).

derived.* metrics (computed here from per_read.tsv, so they do not depend on
evaluate.py's summary): among true insertion reads that the insert finder
called, the % whose primary annotation is IS_DB / Genome / Contamination /
Unannotated (no annotation row, or the note "Unannotated"), and the IS-name
accuracy with all called IS reads as denominator.
"""

import argparse
import concurrent.futures as cf
import copy
import csv
import datetime
import hashlib
import importlib.util
import itertools
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_OUT = os.path.join(REPO, "publication", "benchmark", "grid")
BASELINE_COMMIT = "980d5e6"
INSERTION_CATS = ("is_insertion", "novel_insertion")
STEP_NAMES = ("filter", "map", "insert_finder_batch", "blast_annot_batch", "dr_finder")
NA = "NA"

ACCURACIES = [0.99, 0.97, 0.95, 0.92, 0.90, 0.87, 0.85]
EVALUES = [1e-5, 1e-10, 1e-20, 1e-50]
BLAST_TASKS = ["megablast", "dc-megablast"]

_COMMON = {
    "refs": {"plasmid": "test-data/plasmid.fasta", "genome": "test-data/genome.fasta",
             "is_fasta": "test-data/IS_curated.fasta"},
    "pipelines": {
        "v0.1.0": {"checkout": "baseline", "profile": "v0.1.0"},
        "dev": {"checkout": "dev", "profile": "dev",
                "settings": {"evalue": EVALUES, "blast_task": BLAST_TASKS},
                "default_settings": {"evalue": 1e-10, "blast_task": "megablast"}},
    },
    "run": {"min_len": 2000, "tol": 20},
}

PRESETS = {
    # One-factor-at-a-time sweeps around a centre point (ONT r9.4.1-like, 0.92),
    # plus the accuracy x annotation-settings factorial for the B10 figure.
    "full": dict(copy.deepcopy(_COMMON), name="full", seeds=[1, 2, 3],
                 sim_defaults={"barcodes": 1, "reads": 2000, "insertion_rate": 0.3, "clones": 40,
                               "novel_fraction": 0.15, "fragment_rate": 0.1,
                               "error_model": "ont", "ont_profile": "r9.4.1",
                               "read_accuracy_mean": 0.92, "read_accuracy_sd": "scaled:0.5"},
                 sweeps=[
                     {"name": "accuracy", "sim": {"read_accuracy_mean": ACCURACIES}, "pipelines": "all"},
                     {"name": "chemistry_r10", "set": {"ont_profile": "r10.4.1", "read_accuracy_mean": None,
                                                       "read_accuracy_sd": None}, "pipelines": "all"},
                     {"name": "insertion_rate", "sim": {"insertion_rate": [0.01, 0.05, 0.1, 0.3, 0.5]},
                      "pipelines": "default"},
                     {"name": "depth", "sim": {"reads": [100, 300, 1000, 3000, 10000]}, "pipelines": "default"},
                     {"name": "novel_fraction", "sim": {"novel_fraction": [0.0, 0.15, 0.3, 0.5]},
                      "pipelines": "default"},
                     {"name": "fragmentation", "sim": {"fragment_rate": [0.0, 0.1, 0.3, 0.5]},
                      "pipelines": "default"},
                     {"name": "novel_len", "sim": {"novel_len": [[300, 600], [600, 1200], [1200, 3000]]},
                      "set": {"novel_fraction": 0.5}, "pipelines": "all"},
                 ]),
    # Smoke test: 2 accuracies x 1 seed x 300 reads, both checkouts, 2 E-values.
    "smoke": dict(copy.deepcopy(_COMMON), name="smoke", seeds=[1],
                  sim_defaults={"barcodes": 1, "reads": 300, "error_model": "ont", "ont_profile": "r9.4.1",
                                "read_accuracy_mean": 0.92, "read_accuracy_sd": "scaled:0.5"},
                  sweeps=[{"name": "accuracy", "sim": {"read_accuracy_mean": [0.97, 0.87]},
                           "pipelines": "all"}]),
}
PRESETS["smoke"]["pipelines"]["dev"]["settings"] = {"evalue": [1e-10, 1e-50], "blast_task": ["megablast"]}

KEY_METRICS = [
    "detection_finder.recall", "detection_finder.precision", "detection_finder.f1",
    "detection_annot.recall", "detection_annot.f1",
    "position.error.pct_within_0bp", "position.error.pct_within_5bp",
    "classification.is_insertion.accuracy", "classification.novel_insertion.accuracy",
    "derived.true_ins.pct_Contamination", "derived.true_ins.pct_Unannotated",
    "sites.event_recall", "sites.site_precision",
    "tsd.length_accuracy", "tsd.consistent_accuracy", "runtime.total_s",
]

_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


# ------------------------------------------------------------ identity
def _canon(x):
    """Canonical form for hashing: sorted keys, integral floats as ints."""
    if isinstance(x, dict):
        return {str(k): _canon(v) for k, v in sorted(x.items())}
    if isinstance(x, (list, tuple)):
        return [_canon(v) for v in x]
    if isinstance(x, float) and x.is_integer():
        return int(x)
    return x


def stable_hash(obj, n=12):
    s = json.dumps(_canon(obj), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(s.encode()).hexdigest()[:n]


def load_simulator_profiles():
    spec = importlib.util.spec_from_file_location("_sim_for_grid", os.path.join(HERE, "simulate_traposome.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.ONT_PROFILES


def resolve_sim_params(params, ont_profiles=None):
    """Drop nulls and resolve read_accuracy_sd = "scaled:F"."""
    p = {k: v for k, v in params.items() if v is not None}
    sd = p.get("read_accuracy_sd")
    if isinstance(sd, str) and sd.startswith("scaled:"):
        mean = p.get("read_accuracy_mean")
        if mean is None:
            ont_profiles = ont_profiles or load_simulator_profiles()
            mean = ont_profiles[p.get("ont_profile", "r9.4.1")]["read_accuracy_mean"]
        p["read_accuracy_sd"] = round(float(sd.split(":", 1)[1]) * (1 - float(mean)), 6)
    return p


def nominal_accuracy(sim, ont_profiles=None):
    if sim.get("error_model") != "ont":
        return None
    if sim.get("read_accuracy_mean") is not None:
        return float(sim["read_accuracy_mean"])
    ont_profiles = ont_profiles or load_simulator_profiles()
    return ont_profiles[sim.get("ont_profile", "r9.4.1")]["read_accuracy_mean"]


def variant_label(name, settings):
    if not settings:
        return name
    parts = [name]
    for k in sorted(settings):
        v = settings[k]
        parts.append(f"{k}={v:g}" if isinstance(v, float) else f"{k}={v}")
    return " ".join(parts)


def pipeline_variants(spec, which):
    out = []
    for name, p in spec["pipelines"].items():
        base = {"name": name, "profile": p.get("profile", name), "checkout": p["checkout"], "tag": p.get("tag", "")}
        if not p.get("settings"):
            out.append(dict(base, settings={}))
            continue
        if which == "default":
            out.append(dict(base, settings=dict(p.get("default_settings", {}))))
            continue
        keys = sorted(p["settings"])
        for combo in itertools.product(*(p["settings"][k] for k in keys)):
            out.append(dict(base, settings=dict(zip(keys, combo))))
    return out


def expand_spec(spec, ont_profiles=None):
    """All cells of the grid, de-duplicated; each lists the sweeps it belongs to."""
    cells = {}
    for sw in spec["sweeps"]:
        axes = sw.get("sim", {})
        keys = list(axes)
        combos = list(itertools.product(*(axes[k] for k in keys))) if keys else [()]
        for combo, seed in itertools.product(combos, spec["seeds"]):
            sim = dict(spec["sim_defaults"])
            sim.update(sw.get("set", {}))
            sim.update(dict(zip(keys, combo)))
            sim["seed"] = seed
            sim = resolve_sim_params(sim, ont_profiles)
            sim_full = {"refs": spec["refs"], **sim}
            sim_id = stable_hash(sim_full)
            for var in pipeline_variants(spec, sw.get("pipelines", "default")):
                ident = {"sim": sim_full, "pipeline": {k: var[k] for k in ("name", "profile", "checkout",
                                                                         "settings", "tag")},
                         "run": spec.get("run", {})}
                cid = stable_hash(ident)
                if cid in cells:
                    if sw["name"] not in cells[cid]["sweeps"]:
                        cells[cid]["sweeps"].append(sw["name"])
                    continue
                cells[cid] = {"cell_id": cid, "sim_id": sim_id, "sim": sim_full, "pipeline": var,
                              "variant": variant_label(var["name"], var["settings"]),
                              "run": spec.get("run", {}), "sweeps": [sw["name"]],
                              "accuracy_nominal": nominal_accuracy(sim, ont_profiles)}
    return list(cells.values())


# ------------------------------------------------------------ environment
def _first_line(cmd, env):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=env)
        txt = (p.stdout or "") + (p.stderr or "")
        lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
        return lines[0] if lines else None
    except Exception as e:  # tool missing / broken
        return f"unavailable ({type(e).__name__})"


def tool_versions(python, env=None):
    env = env or os.environ
    v = {"platform": platform.platform(), "python": _first_line([python, "--version"], env),
         "minimap2": _first_line(["minimap2", "--version"], env),
         "samtools": _first_line(["samtools", "--version"], env),
         "blastn": _first_line(["blastn", "-version"], env),
         "makeblastdb": _first_line(["makeblastdb", "-version"], env),
         "mafft": _first_line(["mafft", "--version"], env),
         "biopython": _first_line([python, "-c", "import Bio; print(Bio.__version__)"], env)}
    return v


def git_info(path):
    def g(*a):
        p = subprocess.run(["git", "-C", path, *a], capture_output=True, text=True)
        return p.stdout if p.returncode == 0 else ""
    commit = g("rev-parse", "--short", "HEAD").strip() or None
    # compiled files are tracked in v0.1.0 and change whenever it runs; they are not code changes
    excl = ["--", ".", ":(exclude)*.pyc", ":(exclude)*/__pycache__/*"]
    status = [ln for ln in g("status", "--porcelain", "--untracked-files=no", *excl).splitlines() if ln.strip()]
    diff = g("diff", "HEAD", *excl)
    return {"path": path, "commit": commit, "dirty": bool(status), "n_modified": len(status),
            "dirty_hash": hashlib.sha1(diff.encode()).hexdigest()[:10] if status else ""}


def ensure_baseline(path, commit=BASELINE_COMMIT):
    if os.path.exists(os.path.join(path, "trapo-seq.py")):
        return path
    log(f"creating baseline worktree {path} @ {commit}")
    subprocess.run(["git", "-C", REPO, "worktree", "prune"], check=False)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    subprocess.run(["git", "-C", REPO, "worktree", "add", "--detach", path, commit], check=True)
    return path


def checkout_path(key, ctx):
    if key == "dev":
        return REPO
    if key == "baseline":
        return ctx["baseline_checkout"]
    return os.path.abspath(key)


# ------------------------------------------------------------ simulation
def sim_command(python, sim_full, out_dir):
    refs = sim_full["refs"]
    rp = lambda p: p if os.path.isabs(p) else os.path.join(REPO, p)  # noqa: E731
    cmd = [python, os.path.join(HERE, "simulate_traposome.py"), "-p", rp(refs["plasmid"]),
           "-g", rp(refs["genome"]), "-i", rp(refs["is_fasta"]), "-o", out_dir]
    for k, v in sim_full.items():
        if k == "refs" or v is None:
            continue
        cmd.append(f"--{k}")
        if isinstance(v, (list, tuple)):
            cmd += [str(x) for x in v]
        else:
            cmd.append(str(v))
    return cmd


def sim_fastqs(sim_dir):
    out = []
    if not os.path.isdir(sim_dir):
        return out
    for d in sorted(os.listdir(sim_dir)):
        f = os.path.join(sim_dir, d, f"{d}_raw.fastq")
        if os.path.exists(f):
            out.append(f)
    return out


def sim_is_ready(sim_dir):
    return os.path.exists(os.path.join(sim_dir, "sim_stats.json")) and bool(sim_fastqs(sim_dir))


def truth_stats(truth_tsv):
    cats, acc = {}, []
    with open(truth_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cats[r["category"]] = cats.get(r["category"], 0) + 1
            a = r.get("read_accuracy")
            if a not in (None, "", "N/A"):
                acc.append(float(a))
    st = {"n_reads": sum(cats.values()), "categories": cats}
    if acc:
        st["read_accuracy_mean"] = statistics.fmean(acc)
        st["read_accuracy_median"] = statistics.median(acc)
    return st


def run_simulation(sim_id, sim_full, ctx):
    sim_dir = os.path.join(ctx["out"], "sims", sim_id)
    if sim_is_ready(sim_dir):
        return sim_dir
    if os.path.exists(sim_dir):  # FASTQ removed earlier, or interrupted: regenerate (deterministic)
        shutil.rmtree(sim_dir)
    os.makedirs(sim_dir)
    t0 = time.time()
    cmd = sim_command(ctx["python"], sim_full, sim_dir)
    with open(os.path.join(sim_dir, "simulate.log"), "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=ctx["env"])
    if p.returncode != 0:
        raise RuntimeError(f"simulation {sim_id} failed (see {sim_dir}/simulate.log)")
    stats = truth_stats(os.path.join(sim_dir, "truth.tsv"))
    stats["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(sim_dir, "sim_params.json"), "w") as fh:
        json.dump({"sim_id": sim_id, "params": sim_full, "cmd": cmd}, fh, indent=2)
    with open(os.path.join(sim_dir, "sim_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)
    log(f"sim {sim_id} done in {stats['seconds']} s")
    return sim_dir


# ------------------------------------------------------------ cells
def cell_dir(out, cell_id):
    return os.path.join(out, "cells", cell_id)


def cell_is_done(out, cell_id):
    return os.path.exists(os.path.join(cell_dir(out, cell_id), "metrics.json"))


def blast_steps_json(profile, settings):
    """{"blast_annot_batch": template + options} for run_benchmark --steps_json."""
    if not settings:
        return None
    sys.path.insert(0, HERE)
    import run_benchmark  # noqa: E402  (template source of truth)
    tmpl = run_benchmark.PROFILES[profile]["blast_annot_batch"]
    for k in sorted(settings):
        v = settings[k]
        tmpl += f" --{k} {v:g}" if isinstance(v, float) else f" --{k} {v}"
    return {"blast_annot_batch": tmpl}


def run_cell(cell, ctx):
    """Run run_benchmark.py for one cell. Returns 'done' or raises."""
    out = ctx["out"]
    final = cell_dir(out, cell["cell_id"])
    partial, failed = final + ".partial", final + ".failed"
    for d in (partial, failed):
        if os.path.exists(d):
            shutil.rmtree(d)
    os.makedirs(partial)
    work = os.path.join(ctx["work_root"], cell["cell_id"])
    pipe = cell["pipeline"]
    checkout = checkout_path(pipe["checkout"], ctx)
    sim_dir = os.path.join(out, "sims", cell["sim_id"])
    cmd = [ctx["python"], os.path.join(HERE, "run_benchmark.py"), "--checkout", checkout,
           "--sim_dir", sim_dir, "--work_dir", work, "--out_dir", partial, "--profile", pipe["profile"],
           "--threads", str(ctx["threads"]), "--python", ctx["python"],
           "--min_len", str(cell["run"].get("min_len", 2000)), "--tol", str(cell["run"].get("tol", 20))]
    steps = blast_steps_json(pipe["profile"], pipe["settings"])
    if steps:
        sj = os.path.join(partial, "steps_override.json")
        with open(sj, "w") as fh:
            json.dump(steps, fh, indent=2)
        cmd += ["--steps_json", sj]
    t0 = time.time()
    with open(os.path.join(partial, "runner.log"), "w") as fh:
        fh.write("$ " + " ".join(cmd) + "\n\n")
        fh.flush()
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=ctx["env"])
    rec = {"cell_id": cell["cell_id"], "sim_id": cell["sim_id"], "variant": cell["variant"],
           "sweeps": cell["sweeps"], "params": {"sim": cell["sim"], "pipeline": pipe, "run": cell["run"]},
           "checkout": ctx["checkouts"].get(pipe["checkout"]), "tools": ctx["tools"],
           "threads": ctx["threads"], "returncode": p.returncode,
           "wall_seconds": round(time.time() - t0, 1),
           "finished": datetime.datetime.now().isoformat(timespec="seconds")}
    with open(os.path.join(partial, "cell.json"), "w") as fh:
        json.dump(rec, fh, indent=2)
    ok = p.returncode == 0 and os.path.exists(os.path.join(partial, "metrics.json"))
    if not ok:
        os.rename(partial, failed)
        raise RuntimeError(f"cell {cell['cell_id']} ({cell['variant']}) failed rc={p.returncode}; "
                           f"see {failed}/runner.log (work dir kept: {work})")
    if not ctx["keep"] and os.path.exists(work):
        shutil.rmtree(work, ignore_errors=True)
    if os.path.exists(final):
        shutil.rmtree(final)
    os.rename(partial, final)
    return "done"


def execute(cells, ctx, run_cell_fn=run_cell, run_sim_fn=run_simulation):
    """Run every unfinished cell. Returns {cell_id: status}."""
    status, todo = {}, []
    for c in cells:
        if cell_is_done(ctx["out"], c["cell_id"]):
            status[c["cell_id"]] = "skipped (done)"
        else:
            todo.append(c)
    if ctx.get("max_cells"):
        for c in todo[ctx["max_cells"]:]:
            status[c["cell_id"]] = "pending"
        todo = todo[:ctx["max_cells"]]
    log(f"{len(cells)} cells: {len(cells) - len(todo)} done/deferred, {len(todo)} to run "
        f"(jobs={ctx['jobs']}, threads={ctx['threads']})")
    if not todo:
        return status
    sims = {}
    for c in todo:
        sims.setdefault(c["sim_id"], c["sim"])
    sim_ok = {}
    with cf.ThreadPoolExecutor(max_workers=ctx["jobs"]) as ex:
        futs = {ex.submit(run_sim_fn, sid, s, ctx): sid for sid, s in sims.items()}
        for f in cf.as_completed(futs):
            try:
                f.result()
                sim_ok[futs[f]] = True
            except Exception as e:
                sim_ok[futs[f]] = False
                log(f"ERROR {e}")
    runnable = [c for c in todo if sim_ok.get(c["sim_id"])]
    for c in todo:
        if not sim_ok.get(c["sim_id"]):
            status[c["cell_id"]] = "failed (simulation)"
    n = len(runnable)
    with cf.ThreadPoolExecutor(max_workers=ctx["jobs"]) as ex:
        futs = {ex.submit(run_cell_fn, c, ctx): c for c in runnable}
        for i, f in enumerate(cf.as_completed(futs), 1):
            c = futs[f]
            try:
                status[c["cell_id"]] = f.result()
                log(f"[{i}/{n}] {c['cell_id']} {c['variant']} acc={c['accuracy_nominal']} seed={c['sim']['seed']} done")
            except Exception as e:
                status[c["cell_id"]] = "failed"
                log(f"[{i}/{n}] ERROR {e}")
    if not ctx["keep"]:
        # Simulated FASTQ is regenerable from the parameters; drop it once all cells of a sim are finished.
        for sid in sims:
            if all(cell_is_done(ctx["out"], c["cell_id"]) for c in cells if c["sim_id"] == sid):
                for fq in sim_fastqs(os.path.join(ctx["out"], "sims", sid)):
                    os.remove(fq)
    return status


# ------------------------------------------------------------ aggregation
def flatten(d, prefix=""):
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            yield from flatten(v, key + ".")
        else:
            yield key, v


def to_number(v):
    if isinstance(v, bool):
        return float(v)
    if isinstance(v, (int, float)):
        return float(v) if math.isfinite(v) else None
    if isinstance(v, str):
        try:
            x = float(v)
            return x if math.isfinite(x) else None
        except ValueError:
            return None
    return None


def read_summary_overall(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return out
    idx = rows[0].index("overall") if "overall" in rows[0] else 1
    for r in rows[1:]:
        if len(r) > idx:
            out[r[0]] = to_number(r[idx])
    return out


def derived_metrics(per_read_tsv):
    """Annotation outcome of true insertion reads called by the insert finder."""
    if not os.path.exists(per_read_tsv):
        return {}
    called = {c: [] for c in INSERTION_CATS}
    with open(per_read_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("category") in called and (to_number(r.get("n_finder_rows")) or 0) > 0:
                called[r["category"]].append(r)

    def outcome(r):
        if (to_number(r.get("n_annot_rows")) or 0) == 0 or r.get("annot_note") in ("Unannotated", "", "N/A"):
            return "Unannotated"
        return r.get("annot_note")

    out = {}
    allr = called["is_insertion"] + called["novel_insertion"]
    for scope, rs in (("true_ins", allr), ("is_insertion", called["is_insertion"]),
                      ("novel_insertion", called["novel_insertion"])):
        out[f"derived.{scope}.n_called"] = float(len(rs))
        for lab in ("IS_DB", "Genome", "Contamination", "Unannotated"):
            out[f"derived.{scope}.pct_{lab}"] = (100.0 * sum(outcome(r) == lab for r in rs) / len(rs)) if rs else None
    rs = called["is_insertion"]
    out["derived.is_insertion.name_accuracy_of_called"] = (
        sum(r.get("class_ok") == "True" for r in rs) / len(rs)) if rs else None
    return out


def runtime_metrics(run_log_json):
    if not os.path.exists(run_log_json):
        return {}
    with open(run_log_json) as fh:
        rl = json.load(fh)
    t = {s: 0.0 for s in STEP_NAMES}
    failed = []
    for st in rl.get("steps", []):
        name = st.get("step", "").split(":")[0]
        if name in t:
            t[name] += st.get("seconds") or 0.0
        if st.get("returncode") or st.get("traceback_in_log"):
            failed.append(st.get("step"))
    out = {f"runtime.{k}_s": v for k, v in t.items()}
    out["runtime.total_s"] = sum(t.values())
    out["run.failed_steps"] = float(len(failed))
    return out


def param_columns(cell):
    row = {"cell_id": cell["cell_id"], "sim_id": cell["sim_id"], "variant": cell["variant"],
           "sweeps": ",".join(cell["sweeps"]), "pipe.name": cell["pipeline"]["name"],
           "pipe.profile": cell["pipeline"]["profile"], "pipe.tag": cell["pipeline"].get("tag", "")}
    for k in ("evalue", "blast_task"):
        row[f"pipe.{k}"] = cell["pipeline"]["settings"].get(k, "default")
    for k, v in cell["sim"].items():
        if k == "refs":
            continue
        row[f"sim.{k}"] = "-".join(str(x) for x in v) if isinstance(v, (list, tuple)) else v
    row["sim.accuracy_nominal"] = cell.get("accuracy_nominal")
    return row


def collect_rows(out, cells):
    rows = []
    for c in cells:
        d = cell_dir(out, c["cell_id"])
        mj = os.path.join(d, "metrics.json")
        if not os.path.exists(mj):
            continue
        with open(mj) as fh:
            m = json.load(fh)
        row = param_columns(c)
        cj = os.path.join(d, "cell.json")
        if os.path.exists(cj):
            with open(cj) as fh:
                info = json.load(fh)
            co = info.get("checkout") or {}
            row.update({"pipe.commit": co.get("commit"), "pipe.dirty": co.get("dirty"),
                        "pipe.dirty_hash": co.get("dirty_hash"), "threads": info.get("threads")})
        ss = os.path.join(out, "sims", c["sim_id"], "sim_stats.json")
        metrics = {}
        if os.path.exists(ss):
            with open(ss) as fh:
                st = json.load(fh)
            metrics["sim.read_accuracy_observed"] = st.get("read_accuracy_mean")
            metrics["runtime.simulate_s"] = st.get("seconds")
        metrics.update(read_summary_overall(os.path.join(d, "summary.tsv")))
        for k, v in flatten({k: x for k, x in m.get("overall", {}).items() if k != "files"}):
            if v is not None and not isinstance(v, (int, float)):
                continue  # lists of IDs, confusion strings, ...
            x = to_number(v)
            if x is not None or k not in metrics:
                metrics[k] = x
        metrics.update(derived_metrics(os.path.join(d, "per_read.tsv")))
        metrics.update(runtime_metrics(os.path.join(d, "run_log.json")))
        row["_metrics"] = metrics
        rows.append(row)
    return rows


def fmt(v):
    if v is None:
        return NA
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_tidy(rows, path):
    pcols = []
    for r in rows:
        for k in r:
            if k != "_metrics" and k not in pcols:
                pcols.append(k)
    mcols = []
    for r in rows:
        for k, v in r["_metrics"].items():
            if k not in mcols:
                mcols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(pcols + mcols)
        for r in rows:
            w.writerow([fmt(r.get(k)) for k in pcols] + [fmt(r["_metrics"].get(k)) for k in mcols])
    return pcols, mcols


def summarise(rows, mcols, group_exclude=("cell_id", "sim_id", "sim.seed", "pipe.commit", "pipe.dirty",
                                           "pipe.dirty_hash", "threads", "sweeps")):
    groups = {}
    for r in rows:
        key = tuple((k, fmt(r[k])) for k in r if k != "_metrics" and k not in group_exclude)
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in groups.items():
        g = dict(key)
        g["sweeps"] = ",".join(sorted({s for r in rs for s in r.get("sweeps", "").split(",") if s}))
        g["n_seeds"] = len(rs)
        commits = {fmt(r.get("pipe.commit")) + ("+dirty" if r.get("pipe.dirty") else "") for r in rs}
        g["commits"] = ",".join(sorted(commits))
        stats = {}
        for m in mcols:
            xs = [r["_metrics"].get(m) for r in rs]
            xs = [x for x in xs if isinstance(x, float)]
            stats[m] = (statistics.fmean(xs) if xs else None, statistics.stdev(xs) if len(xs) > 1 else None, len(xs))
        out.append((g, stats))
    return out


def write_summaries(summary, mcols, out):
    gcols = []
    for g, _ in summary:
        for k in g:
            if k not in gcols:
                gcols.append(k)
    with open(os.path.join(out, "summary_long.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(gcols + ["metric", "mean", "sd", "n"])
        for g, st in summary:
            for m in mcols:
                mean, sd, n = st[m]
                w.writerow([fmt(g.get(k)) for k in gcols] + [m, fmt(mean), fmt(sd), n])
    keys = [m for m in KEY_METRICS if m in mcols]
    with open(os.path.join(out, "summary_key.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(gcols + keys)
        for g, st in summary:
            cells = []
            for m in keys:
                mean, sd, n = st[m]
                cells.append(NA if mean is None else (f"{mean:.4g}" if sd is None else f"{mean:.4g} ± {sd:.2g}"))
            w.writerow([fmt(g.get(k)) for k in gcols] + cells)


def aggregate(out, cells):
    rows = collect_rows(out, cells)
    if not rows:
        log("no finished cells to aggregate")
        return None
    pcols, mcols = write_tidy(rows, os.path.join(out, "grid_results.tsv"))
    summary = summarise(rows, mcols)
    write_summaries(summary, mcols, out)
    mixed = {}
    for r in rows:
        mixed.setdefault(r["variant"], set()).add((r.get("pipe.commit"), r.get("pipe.dirty_hash")))
    for v, s in mixed.items():
        if len(s) > 1:
            log(f"WARNING: variant '{v}' mixes code versions {sorted(map(str, s))}")
    log(f"aggregated {len(rows)} cells, {len(mcols)} metrics -> {out}/grid_results.tsv, summary_*.tsv")
    return os.path.join(out, "grid_results.tsv")


# ------------------------------------------------------------ manifest
def write_manifest(out, spec, cells, status, ctx):
    path = os.path.join(out, "manifest.json")
    rows = []
    for c in cells:
        st = "done" if cell_is_done(out, c["cell_id"]) else status.get(c["cell_id"], "pending")
        rows.append({"cell_id": c["cell_id"], "sim_id": c["sim_id"], "variant": c["variant"],
                     "sweeps": c["sweeps"], "status": st, "sim": c["sim"], "pipeline": c["pipeline"],
                     "run": c["run"]})
    man = {"updated": datetime.datetime.now().isoformat(timespec="seconds"), "spec": spec,
           "tools": ctx.get("tools"), "checkouts": ctx.get("checkouts"), "python": ctx.get("python"),
           "threads": ctx.get("threads"), "jobs": ctx.get("jobs"), "work_root": ctx.get("work_root"),
           "n_cells": len(cells), "n_done": sum(r["status"] == "done" for r in rows), "cells": rows}
    with open(path, "w") as fh:
        json.dump(man, fh, indent=2)
    return path


# ------------------------------------------------------------ main
def load_spec(args):
    if args.spec:
        with open(args.spec) as fh:
            spec = json.load(fh)
    else:
        spec = copy.deepcopy(PRESETS[args.preset])
    if args.seeds:
        spec["seeds"] = args.seeds
    return spec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preset", choices=sorted(PRESETS), default="full", help="Built-in grid (default full)")
    ap.add_argument("--spec", help="JSON grid specification (overrides --preset)")
    ap.add_argument("--dump_spec", action="store_true", help="Print the grid specification as JSON and exit")
    ap.add_argument("--out", default=DEFAULT_OUT, help="Output root (default publication/benchmark/grid)")
    ap.add_argument("--work_root", default=os.path.join(tempfile.gettempdir(), "trapo-seq-grid-work"),
                    help="Scratch for pipeline work dirs (no spaces if v0.1.0 is in the grid)")
    ap.add_argument("--baseline_checkout", default=os.path.join(tempfile.gettempdir(), "trapo-seq-grid",
                                                                "baseline-v0.1.0"),
                    help=f"v0.1.0 checkout; created with `git worktree add` at {BASELINE_COMMIT} if missing")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--jobs", type=int, default=2, help="Cells run in parallel (default 2)")
    ap.add_argument("--threads", type=int, default=4, help="Threads per pipeline run (default 4)")
    ap.add_argument("--seeds", type=int, nargs="+", help="Override the spec's seeds")
    ap.add_argument("--keep", action="store_true", help="Keep work dirs (BAM/FASTQ) and simulated FASTQ")
    ap.add_argument("--max_cells", type=int, default=0, help="Run at most N unfinished cells this time")
    ap.add_argument("--dry_run", action="store_true", help="List cells and their status, run nothing")
    ap.add_argument("--aggregate_only", action="store_true", help="Only aggregate finished cells (and plot)")
    ap.add_argument("--no_plot", action="store_true", help="Do not call plot_grid.py")
    args = ap.parse_args(argv)

    spec = load_spec(args)
    if args.dump_spec:
        print(json.dumps(spec, indent=2))
        return 0
    out = os.path.abspath(args.out)
    cells = expand_spec(spec)
    os.makedirs(out, exist_ok=True)

    if args.dry_run:
        n_done = 0
        for c in cells:
            done = cell_is_done(out, c["cell_id"])
            n_done += done
            sim = c["sim"]
            print(f"{c['cell_id']}  {'done   ' if done else 'pending'}  sim={c['sim_id']}  seed={sim['seed']}  "
                  f"acc={c['accuracy_nominal']}  reads={sim.get('reads')}  {c['variant']}  [{','.join(c['sweeps'])}]")
        print(f"\n{len(cells)} cells ({len({c['sim_id'] for c in cells})} simulations), {n_done} done")
        return 0

    env = dict(os.environ, MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    ctx = {"out": out, "python": args.python, "jobs": max(1, args.jobs), "threads": args.threads,
           "keep": args.keep, "max_cells": args.max_cells, "env": env,
           "work_root": os.path.abspath(args.work_root), "baseline_checkout": os.path.abspath(args.baseline_checkout)}
    status = {}
    if not args.aggregate_only:
        keys = {c["pipeline"]["checkout"] for c in cells}
        if any(c["pipeline"]["profile"] == "v0.1.0" for c in cells) and " " in ctx["work_root"]:
            sys.exit("--work_root must not contain spaces when v0.1.0 is in the grid (AUDIT B18)")
        if "baseline" in keys:
            ensure_baseline(ctx["baseline_checkout"])
            if " " in ctx["baseline_checkout"]:
                sys.exit("--baseline_checkout must not contain spaces (AUDIT B18)")
        ctx["checkouts"] = {k: git_info(checkout_path(k, ctx)) for k in keys}
        ctx["tools"] = tool_versions(args.python, env)
        for k, info in ctx["checkouts"].items():
            log(f"checkout {k}: {info['path']} @ {info['commit']}{' (DIRTY: %d files)' % info['n_modified'] if info['dirty'] else ''}")
        os.makedirs(ctx["work_root"], exist_ok=True)
        write_manifest(out, spec, cells, {}, ctx)
        status = execute(cells, ctx)
        write_manifest(out, spec, cells, status, ctx)
        n_fail = sum(s.startswith("failed") for s in status.values())
        if n_fail:
            log(f"{n_fail} cells failed; re-run the same command to retry them")
    results = aggregate(out, cells)
    if results and not args.no_plot:
        subprocess.run([args.python, os.path.join(HERE, "plot_grid.py"), "--results", results,
                        "--out_dir", os.path.join(out, "figures")], env=dict(os.environ, MPLBACKEND="Agg"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
