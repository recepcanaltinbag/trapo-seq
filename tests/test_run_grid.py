"""Unit tests for simulate/run_grid.py: cell identity, resume/skip, aggregation.

Fast: no simulator or pipeline runs (runner functions are replaced by fakes).
"""
import copy
import csv
import json
import os

import pytest

import run_grid as rg

from conftest import REPO, TEST_DATA, requires_tools

PROFILES = {"r9.4.1": {"read_accuracy_mean": 0.92}, "r10.4.1": {"read_accuracy_mean": 0.98}}


def tiny_spec():
    return {
        "name": "t", "refs": {"plasmid": "p.fa", "genome": "g.fa", "is_fasta": "i.fa"},
        "seeds": [1, 2],
        "sim_defaults": {"reads": 100, "error_model": "ont", "ont_profile": "r9.4.1",
                         "read_accuracy_mean": 0.92, "read_accuracy_sd": "scaled:0.5"},
        "pipelines": {"v0.1.0": {"checkout": "baseline", "profile": "v0.1.0"},
                      "dev": {"checkout": "dev", "profile": "dev",
                              "settings": {"evalue": [1e-5, 1e-10], "blast_task": ["megablast"]},
                              "default_settings": {"evalue": 1e-10, "blast_task": "megablast"}}},
        "sweeps": [{"name": "accuracy", "sim": {"read_accuracy_mean": [0.97, 0.92]}, "pipelines": "all"},
                   {"name": "depth", "sim": {"reads": [100, 300]}, "pipelines": "default"},
                   {"name": "r10", "set": {"ont_profile": "r10.4.1", "read_accuracy_mean": None,
                                           "read_accuracy_sd": None}, "pipelines": "default"}],
        "run": {"min_len": 2000, "tol": 20},
    }


# ---------------------------------------------------------------- naming
def test_hash_is_order_independent_and_canonical():
    a = {"x": 1, "y": [1.0, 2], "z": {"b": 0.3, "a": "s"}}
    b = {"z": {"a": "s", "b": 0.3}, "y": [1, 2], "x": 1.0}
    assert rg.stable_hash(a) == rg.stable_hash(b)
    assert len(rg.stable_hash(a)) == 12
    assert rg.stable_hash({"evalue": 1e-5}) != rg.stable_hash({"evalue": 1e-10})


def test_scaled_sd_and_nulls():
    p = rg.resolve_sim_params({"read_accuracy_mean": 0.9, "read_accuracy_sd": "scaled:0.5", "x": None}, PROFILES)
    assert p == {"read_accuracy_mean": 0.9, "read_accuracy_sd": 0.05}
    p = rg.resolve_sim_params({"ont_profile": "r10.4.1", "read_accuracy_sd": "scaled:0.5"}, PROFILES)
    assert p["read_accuracy_sd"] == pytest.approx(0.01)


def test_expand_spec_counts_dedup_and_sweeps():
    cells = rg.expand_spec(tiny_spec(), PROFILES)
    # accuracy: 2 acc x 2 seeds x (v0.1.0 + 2 dev) = 12; depth: 2 x 2 x 2 = 8, of which
    # reads=100 at acc 0.92 duplicates the accuracy sweep's centre (4 cells); r10: 2 seeds x 2 = 4
    assert len(cells) == 12 + 8 - 4 + 4
    ids = [c["cell_id"] for c in cells]
    assert len(set(ids)) == len(ids)
    shared = [c for c in cells if c["sweeps"] == ["accuracy", "depth"]]
    assert len(shared) == 4 and all(c["sim"]["read_accuracy_mean"] == 0.92 for c in shared)
    # every pipeline variant of one dataset shares the simulation
    by_sim = {}
    for c in cells:
        by_sim.setdefault(c["sim_id"], set()).add(c["variant"])
    assert max(len(v) for v in by_sim.values()) == 3
    r10 = [c for c in cells if "r10" in c["sweeps"]]
    assert all("read_accuracy_mean" not in c["sim"] and c["accuracy_nominal"] == 0.98 for c in r10)
    # deterministic across calls and independent of thread/out settings
    assert ids == [c["cell_id"] for c in rg.expand_spec(tiny_spec(), PROFILES)]


def test_cell_id_changes_with_settings_tag_and_seed():
    spec = tiny_spec()
    ids = {c["cell_id"] for c in rg.expand_spec(spec, PROFILES)}
    spec2 = copy.deepcopy(spec)
    spec2["pipelines"]["dev"]["tag"] = "after-fix"
    ids2 = {c["cell_id"] for c in rg.expand_spec(spec2, PROFILES)}
    assert len(ids & ids2) == 8  # only the v0.1.0 cells are unchanged
    spec3 = copy.deepcopy(spec)
    spec3["seeds"] = [3]
    assert not ids & {c["cell_id"] for c in rg.expand_spec(spec3, PROFILES)}


def test_variant_label_and_blast_steps():
    assert rg.variant_label("v0.1.0", {}) == "v0.1.0"
    assert rg.variant_label("dev", {"evalue": 1e-20, "blast_task": "dc-megablast"}) == \
        "dev blast_task=dc-megablast evalue=1e-20"
    steps = rg.blast_steps_json("dev", {"evalue": 1e-20, "blast_task": "dc-megablast"})
    assert steps["blast_annot_batch"].endswith("--blast_task dc-megablast --evalue 1e-20")
    assert rg.blast_steps_json("v0.1.0", {}) is None


def test_sim_command_lists_and_refs(tmp_path):
    cmd = rg.sim_command("py", {"refs": {"plasmid": "/a/p.fa", "genome": "g.fa", "is_fasta": "i.fa"},
                                "novel_len": [300, 600], "seed": 2, "skip": None}, str(tmp_path))
    i = cmd.index("--novel_len")
    assert cmd[i + 1:i + 3] == ["300", "600"]
    assert cmd[cmd.index("-p") + 1] == "/a/p.fa"
    assert os.path.isabs(cmd[cmd.index("-g") + 1])
    assert "--skip" not in cmd


# ---------------------------------------------------------------- fakes
def fake_metrics(acc, is_contam_pct=0.0):
    return {"meta": {}, "per_barcode": {}, "overall": {
        "detection_finder": {"precision": 1.0, "recall": acc, "f1": acc, "available": True},
        "detection_annot": {"precision": 1.0, "recall": acc - 0.01, "f1": acc - 0.01},
        "position": {"error": {"median": 0.0, "pct_within_0bp": 100 * acc, "pct_within_5bp": 99.0}},
        "classification": {"is_insertion": {"accuracy": acc, "pct_note_Contamination": is_contam_pct,
                                            "top_confusions": ["a->b:1"]},
                           "contamination_note_by_truth": {"deletion": 1}},
        "tsd": {"length_accuracy": 0.5, "consistent_accuracy": 0.7},
        "files": {"x": "y"}}}


PER_READ_COLS = ["read_id", "category", "n_finder_rows", "n_annot_rows", "annot_note", "class_ok"]


def write_fake_cell(out, cell, acc, notes):
    d = rg.cell_dir(out, cell["cell_id"])
    os.makedirs(d)
    with open(os.path.join(d, "metrics.json"), "w") as fh:
        json.dump(fake_metrics(acc), fh)
    with open(os.path.join(d, "summary.tsv"), "w") as fh:
        fh.write(f"metric\toverall\tbarcode01\ndetection_finder.f1\t{acc:.4f}\t{acc:.4f}\n"
                 "sites.event_recall\tN/A\tN/A\n")
    with open(os.path.join(d, "per_read.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(PER_READ_COLS)
        for i, (cat, nf, na, note, ok) in enumerate(notes):
            w.writerow([f"r{i}", cat, nf, na, note, ok])
    with open(os.path.join(d, "run_log.json"), "w") as fh:
        json.dump({"steps": [{"step": "filter:barcode01", "seconds": 1.0, "returncode": 0},
                             {"step": "map", "seconds": 2.0, "returncode": 0},
                             {"step": "blast_annot_batch", "seconds": 3.5, "returncode": 1}]}, fh)
    with open(os.path.join(d, "cell.json"), "w") as fh:
        json.dump({"checkout": {"commit": "abc1234", "dirty": False, "dirty_hash": ""}, "threads": 4}, fh)


NOTES = [("is_insertion", 1, 1, "IS_DB", "True"), ("is_insertion", 1, 1, "Contamination", "False"),
         ("is_insertion", 1, 0, "N/A", "N/A"), ("novel_insertion", 1, 1, "Genome", "True"),
         ("is_insertion", 0, 0, "N/A", "N/A"),  # not called: excluded from derived.*
         ("empty_plasmid", 0, 0, "N/A", "N/A")]


def ctx_for(out):
    return {"out": str(out), "jobs": 2, "threads": 1, "keep": False, "max_cells": 0,
            "work_root": os.path.join(str(out), "work")}


# ---------------------------------------------------------------- resume
def test_resume_skips_finished_cells(tmp_path):
    cells = rg.expand_spec(tiny_spec(), PROFILES)
    done = cells[:5]
    for c in done:
        write_fake_cell(str(tmp_path), c, 0.9, NOTES)
    ran, sims = [], []

    def fake_sim(sid, sim, ctx):
        sims.append(sid)

    def fake_cell(cell, ctx):
        ran.append(cell["cell_id"])
        write_fake_cell(ctx["out"], cell, 0.8, NOTES)
        return "done"

    status = rg.execute(cells, ctx_for(tmp_path), fake_cell, fake_sim)
    assert sorted(ran) == sorted(c["cell_id"] for c in cells[5:])
    assert all(status[c["cell_id"]].startswith("skipped") for c in done)
    assert len(sims) == len({c["sim_id"] for c in cells[5:]})  # one simulation per dataset
    ran.clear()
    rg.execute(cells, ctx_for(tmp_path), fake_cell, fake_sim)
    assert ran == []  # second invocation: nothing left to do


def test_failed_cell_is_reported_and_retried(tmp_path):
    cells = rg.expand_spec(tiny_spec(), PROFILES)[:3]
    calls = []

    def flaky(cell, ctx):
        calls.append(cell["cell_id"])
        if cell is cells[0] and calls.count(cell["cell_id"]) == 1:
            raise RuntimeError("boom")
        write_fake_cell(ctx["out"], cell, 0.9, NOTES)
        return "done"

    st = rg.execute(cells, ctx_for(tmp_path), flaky, lambda *a: None)
    assert st[cells[0]["cell_id"]] == "failed"
    st = rg.execute(cells, ctx_for(tmp_path), flaky, lambda *a: None)
    assert st[cells[0]["cell_id"]] == "done"
    assert calls.count(cells[0]["cell_id"]) == 2


def test_max_cells_limits_work(tmp_path):
    cells = rg.expand_spec(tiny_spec(), PROFILES)
    ran = []

    def fake_cell(cell, ctx):
        ran.append(1)
        write_fake_cell(ctx["out"], cell, 0.9, NOTES)
        return "done"

    ctx = ctx_for(tmp_path)
    ctx["max_cells"] = 3
    st = rg.execute(cells, ctx, fake_cell, lambda *a: None)
    assert len(ran) == 3 and list(st.values()).count("pending") == len(cells) - 3


# ---------------------------------------------------------------- aggregation
def test_derived_metrics(tmp_path):
    c = rg.expand_spec(tiny_spec(), PROFILES)[0]
    write_fake_cell(str(tmp_path), c, 0.9, NOTES)
    d = rg.derived_metrics(os.path.join(rg.cell_dir(str(tmp_path), c["cell_id"]), "per_read.tsv"))
    assert d["derived.true_ins.n_called"] == 4
    assert d["derived.true_ins.pct_Contamination"] == pytest.approx(25.0)
    assert d["derived.true_ins.pct_Unannotated"] == pytest.approx(25.0)
    assert d["derived.true_ins.pct_IS_DB"] == pytest.approx(25.0)
    assert d["derived.is_insertion.name_accuracy_of_called"] == pytest.approx(1 / 3)
    assert d["derived.novel_insertion.pct_Genome"] == pytest.approx(100.0)


def test_aggregate_tidy_and_summary(tmp_path):
    out = str(tmp_path)
    cells = [c for c in rg.expand_spec(tiny_spec(), PROFILES) if c["sweeps"] == ["accuracy"]]
    accs = {}
    for i, c in enumerate(cells):
        a = 0.8 + 0.1 * (c["sim"]["seed"] - 1)  # seed 1 -> 0.8, seed 2 -> 0.9
        accs[c["cell_id"]] = a
        write_fake_cell(out, c, a, NOTES)
    path = rg.aggregate(out, cells)
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == len(cells)
    r = rows[0]
    for col in ("cell_id", "variant", "pipe.evalue", "sim.seed", "sim.accuracy_nominal",
                "detection_finder.f1", "position.error.pct_within_0bp", "classification.is_insertion.accuracy",
                "derived.true_ins.pct_Contamination", "runtime.map_s", "runtime.total_s", "run.failed_steps"):
        assert col in r, col
    assert "classification.is_insertion.top_confusions" not in r  # non-numeric leaves dropped
    assert "files.x" not in r
    assert float(r["detection_finder.f1"]) == pytest.approx(accs[r["cell_id"]], abs=1e-4)
    assert r["sites.event_recall"] == "NA"
    assert float(r["runtime.filter_s"]) == 1.0 and float(r["runtime.total_s"]) == 6.5
    assert float(r["run.failed_steps"]) == 1
    # summary over seeds: mean 0.85, sd of (0.8, 0.9)
    with open(os.path.join(out, "summary_long.tsv")) as fh:
        long = [x for x in csv.DictReader(fh, delimiter="\t") if x["metric"] == "detection_finder.f1"]
    assert len(long) == len(cells) // 2
    assert all(x["n"] == "2" for x in long)
    assert float(long[0]["mean"]) == pytest.approx(0.85, abs=1e-4)
    assert float(long[0]["sd"]) == pytest.approx(0.0707107, abs=1e-4)
    with open(os.path.join(out, "summary_key.tsv")) as fh:
        key = list(csv.DictReader(fh, delimiter="\t"))
    assert "±" in key[0]["detection_finder.f1"]
    assert "sim.seed" not in key[0] and key[0]["n_seeds"] == "2"


def test_aggregate_ignores_unfinished(tmp_path):
    cells = rg.expand_spec(tiny_spec(), PROFILES)[:2]
    write_fake_cell(str(tmp_path), cells[0], 0.9, NOTES)
    os.makedirs(rg.cell_dir(str(tmp_path), cells[1]["cell_id"]) + ".partial")
    path = rg.aggregate(str(tmp_path), cells)
    with open(path) as fh:
        assert len(list(csv.DictReader(fh, delimiter="\t"))) == 1


def test_manifest_records_status(tmp_path):
    cells = rg.expand_spec(tiny_spec(), PROFILES)[:2]
    write_fake_cell(str(tmp_path), cells[0], 0.9, NOTES)
    p = rg.write_manifest(str(tmp_path), tiny_spec(), cells, {cells[1]["cell_id"]: "failed"},
                          {"tools": {"minimap2": "2.31"}, "checkouts": {}})
    with open(p) as fh:
        m = json.load(fh)
    assert [c["status"] for c in m["cells"]] == ["done", "failed"]
    assert m["n_done"] == 1 and m["tools"]["minimap2"] == "2.31"


# ---------------------------------------------------------------- grouping (fakes)
def test_group_key_shares_upstream_only_across_annotation_settings():
    cells = rg.expand_spec(tiny_spec(), PROFILES)
    groups = {}
    for c in cells:
        groups.setdefault(rg.group_key(c), []).append(c)
    for cs in groups.values():
        assert len({(c["sim_id"], c["pipeline"]["checkout"], c["pipeline"]["profile"]) for c in cs}) == 1
    # per dataset: one v0.1.0 group + one dev group
    assert len(groups) == 2 * len({c["sim_id"] for c in cells})
    assert max(len(cs) for cs in groups.values()) == 2


def _grouping_run(tmp_path, make_snapshot=True, group=True):
    cells = [c for c in rg.expand_spec(tiny_spec(), PROFILES) if c["sweeps"] == ["accuracy"]]
    ctx = ctx_for(tmp_path)
    ctx["group"] = group
    calls, finished = [], set()

    def fake_cell(cell, ctx):
        up = dict((ctx.get("upstream") or {}).get(cell["cell_id"]) or {"role": "single"})
        calls.append((cell["cell_id"], up))
        if up["role"] == "follower":
            lead = [cid for cid, u in calls if u.get("role") == "lead" and u["group"] == up["group"]]
            assert lead and lead[0] in finished              # followers start after their lead
            assert os.path.exists(os.path.join(up["snapshot"], ".run_benchmark_stamp.json"))
        if up["role"] == "lead" and make_snapshot:
            os.makedirs(up["snapshot"])
            open(os.path.join(up["snapshot"], ".run_benchmark_stamp.json"), "w").close()
        write_fake_cell(ctx["out"], cell, 0.9, NOTES)
        finished.add(cell["cell_id"])
        return "done"

    st = rg.execute(cells, ctx, fake_cell, lambda *a: None)
    assert all(v == "done" for v in st.values()) and len(calls) == len(cells)
    return cells, calls, ctx


def test_grouping_lead_then_followers_and_cleanup(tmp_path):
    cells, calls, ctx = _grouping_run(tmp_path)
    roles = [u["role"] for _, u in calls]
    groups = {}
    for c in cells:
        groups.setdefault(rg.group_key(c), []).append(c)
    multi = [cs for cs in groups.values() if len(cs) > 1]
    assert multi and roles.count("lead") == len(multi)
    assert roles.count("follower") == sum(len(cs) - 1 for cs in multi)
    assert roles.count("single") == len(groups) - len(multi)   # v0.1.0, lone dev cells
    assert not [d for d in os.listdir(ctx["work_root"]) if d.startswith("upstream-")]  # snapshots removed


def test_grouping_without_snapshot_runs_from_scratch(tmp_path):
    _, calls, _ = _grouping_run(tmp_path, make_snapshot=False)
    assert "follower" not in [u["role"] for _, u in calls]


def test_no_group_runs_every_cell_alone(tmp_path):
    _, calls, _ = _grouping_run(tmp_path, group=False)
    assert {u["role"] for _, u in calls} == {"single"}


def test_runtime_metrics_cpu_columns(tmp_path):
    p = tmp_path / "run_log.json"
    p.write_text(json.dumps({"reused_steps": ["filter", "map"], "steps": [
        {"step": "filter:barcode01", "seconds": 1.0, "cpu_seconds": 0.5, "returncode": 0, "reused": True},
        {"step": "filter:barcode02", "seconds": 1.0, "cpu_seconds": 0.25, "returncode": 0, "reused": True},
        {"step": "map", "seconds": 2.0, "cpu_seconds": 6.0, "returncode": 0, "reused": True},
        {"step": "blast_annot_batch", "seconds": 3.0, "cpu_seconds": 2.0, "returncode": 0}]}))
    m = rg.runtime_metrics(str(p))
    assert m["runtime.filter_cpu_s"] == 0.75 and m["runtime.map_cpu_s"] == 6.0
    assert m["runtime.total_cpu_s"] == 8.75 and m["runtime.total_s"] == 7.0
    assert m["run.upstream_reused"] == 1.0
    p.write_text(json.dumps({"steps": [{"step": "map", "seconds": 2.0, "returncode": 0}]}))
    m = rg.runtime_metrics(str(p))   # old run logs: no CPU columns
    assert "runtime.total_cpu_s" not in m and m["run.upstream_reused"] == 0.0


def test_full_preset_is_variation_sweeps():
    cells = rg.expand_spec(copy.deepcopy(rg.PRESETS["full"]), PROFILES)
    div = [c for c in cells if "is_divergence" in c["sweeps"]]
    trunc = [c for c in cells if "is_truncation" in c["sweeps"]]
    n_var = 1 + len(rg.EVALUES) * len(rg.BLAST_TASKS)   # v0.1.0 + every dev annotation setting
    seeds = len(rg.PRESETS["full"]["seeds"])
    assert len(div) == len(rg.IS_DIVERGENCES) * seeds * n_var
    assert len(trunc) == len(rg.IS_TRUNCATIONS) * seeds * n_var
    # divergence 0 is the accuracy sweep's centre point, not a separate dataset
    zero = [c for c in div if "is_divergence" not in c["sim"]]
    assert len(zero) == seeds * n_var and all("accuracy" in c["sweeps"] for c in zero)
    assert {c["sim"].get("is_divergence") for c in div} == {None, 0.05, 0.10, 0.15, 0.20}
    assert all(c["sim"]["is_truncate_frac"] == 0.5 for c in trunc)
    cmd = rg.sim_command("py", trunc[0]["sim"], "/x")
    i = cmd.index("--is_truncate")
    assert cmd[i + 1:i + 3] in (["300", "600"], ["600", "1200"])


# ---------------------------------------------------------------- grouping (real pipeline)
@requires_tools
def test_grouped_and_ungrouped_runs_give_identical_metrics(tmp_path):
    spec = {"name": "grouping", "seeds": [1],
            "refs": {k: os.path.join(TEST_DATA, f) for k, f in
                     (("plasmid", "plasmid.fasta"), ("genome", "genome.fasta"), ("is_fasta", "IS_curated.fasta"))},
            "sim_defaults": {"barcodes": 1, "reads": 200, "error_model": "ont", "ont_profile": "r9.4.1",
                             "is_divergence": 0.15},
            "pipelines": {"dev": {"checkout": "dev", "profile": "dev",
                                  "settings": {"evalue": [1e-10], "blast_task": ["megablast", "dc-megablast"]},
                                  "default_settings": {"evalue": 1e-10, "blast_task": "megablast"}}},
            "sweeps": [{"name": "s", "pipelines": "all"}], "run": {"min_len": 2000, "tol": 20}}
    sp = tmp_path / "spec.json"
    sp.write_text(json.dumps(spec))
    outs = {}
    for mode, extra in (("grouped", []), ("ungrouped", ["--no_group"])):
        out = tmp_path / mode
        rg.main(["--spec", str(sp), "--out", str(out), "--work_root", str(tmp_path / f"work_{mode}"),
                 "--jobs", "2", "--threads", "2", "--no_plot"] + extra)
        outs[mode] = out
    cells = rg.expand_spec(spec)
    assert len(cells) == 2
    roles = []
    for c in cells:
        res = {}
        for mode, out in outs.items():
            d = rg.cell_dir(str(out), c["cell_id"])
            with open(os.path.join(d, "metrics.json")) as fh:
                m = json.load(fh)["overall"]
            m.pop("files", None)
            with open(os.path.join(d, "per_read.tsv")) as fh:
                res[mode] = (m, fh.read())
            with open(os.path.join(d, "cell.json")) as fh:
                info = json.load(fh)
            if mode == "grouped":
                roles.append(info["upstream"]["role"])
                if info["upstream"]["role"] == "follower":
                    assert info["upstream"]["reused"] is True
            else:
                assert info["upstream"]["role"] == "single"
        assert res["grouped"] == res["ungrouped"], c["variant"]
    assert sorted(roles) == ["follower", "lead"]
    with open(outs["grouped"] / "grid_results.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert all(float(r["runtime.total_cpu_s"]) > 0 for r in rows)
    assert sorted(r["run.upstream_reused"] for r in rows) == ["0", "1"]
    assert not [d for d in os.listdir(tmp_path / "work_grouped") if d.startswith("upstream-")]
