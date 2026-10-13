"""simulate/run_benchmark.py: CPU time per step and reuse of upstream work
(--reuse_work / --from_step / --snapshot_after) giving the same results as a
run from scratch."""
import json
import os
import subprocess
import sys

import pytest

import run_benchmark as rb

from conftest import REPO, requires_tools

RB = os.path.join(REPO, "simulate", "run_benchmark.py")
DC = "{py} {cli} blast_annot_batch -d {data} -g {genome} --is_fasta {is_fasta} --temp {temp} " \
     "--threads {threads} --blast_task dc-megablast --evalue 1e-20"


def bench(sim, work, out, *extra, check=True):
    p = subprocess.run([sys.executable, RB, "--checkout", REPO, "--sim_dir", str(sim), "--work_dir", str(work),
                        "--out_dir", str(out), "--threads", "2", *map(str, extra)],
                       capture_output=True, text=True)
    if check:
        assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    return p


def results(out):
    with open(os.path.join(out, "metrics.json")) as fh:
        m = json.load(fh)["overall"]
    m.pop("files", None)
    tables = {f: open(os.path.join(out, f)).read() for f in ("per_read.tsv", "sites.tsv", "summary.tsv")}
    return m, tables


def run_log(out):
    with open(os.path.join(out, "run_log.json")) as fh:
        return json.load(fh)


# ------------------------------------------------------------------ unit
def test_manifest_and_check_reuse(tmp_path):
    w = tmp_path / "w"
    (w / "data" / "bc").mkdir(parents=True)
    (w / "data" / "bc" / "a.txt").write_text("a")
    ident = {"sim_dir": "s"}
    tmpl = {s: s for s in rb.STEP_ORDER}
    stamp = {"version": rb.STAMP_VERSION, "identity": ident, "templates": tmpl, "barcodes": ["bc"],
             "done": [], "records": {}, "manifests": {"prepare": rb.manifest(str(w))}}
    (w / "data" / "bc" / "b.txt").write_text("b")
    stamp["done"] = ["filter"]
    stamp["records"]["filter"] = [{"step": "filter:bc", "seconds": 1.0, "cpu_seconds": 0.5}]
    stamp["manifests"]["filter"] = rb.manifest(str(w))
    (w / "data" / "out").mkdir()
    (w / "data" / "out" / "c.txt").write_text("c")    # written by a later step
    stamp["done"].append("map")
    stamp["manifests"]["map"] = rb.manifest(str(w))
    rb.write_stamp(str(w), stamp)
    assert "data/bc/" in stamp["manifests"]["prepare"] and rb.STAMP not in stamp["manifests"]["map"]

    with pytest.raises(rb.ReuseError, match="different run"):
        rb.check_reuse(str(w), {"sim_dir": "other"}, tmpl, "map")
    with pytest.raises(rb.ReuseError, match="different command template"):
        rb.check_reuse(str(w), ident, dict(tmpl, filter="x"), "map")
    with pytest.raises(rb.ReuseError, match="never completed"):
        rb.check_reuse(str(w), ident, tmpl, "blast_annot_batch")   # insert_finder_batch missing
    # resume after the last completed step; nothing is removed
    _, reused, start = rb.check_reuse(str(w), ident, tmpl, None)
    assert reused == ["filter", "map"] and start == "insert_finder_batch"
    # from map: later files and directories are removed
    st, reused, start = rb.check_reuse(str(w), ident, tmpl, "map")
    assert reused == ["filter"] and start == "map" and st["done"] == ["filter"]
    assert not (w / "data" / "out").exists() and (w / "data" / "bc" / "b.txt").exists()
    # a changed upstream file is refused
    (w / "data" / "bc" / "b.txt").write_text("changed")
    with pytest.raises(rb.ReuseError, match="missing or changed"):
        rb.check_reuse(str(w), ident, tmpl, "map")


def test_copy_work_preserves_manifest(tmp_path):
    w = tmp_path / "w"
    (w / "d").mkdir(parents=True)
    (w / "d" / "f").write_text("x" * 1000)
    rb.copy_work(str(w), str(tmp_path / "c"))
    assert rb.manifest(str(w)) == rb.manifest(str(tmp_path / "c"))


# ------------------------------------------------------------------ end to end
@requires_tools
def test_from_step_equals_scratch_and_records_cpu(sim_dir, tmp_path):
    dc = tmp_path / "dc.json"
    dc.write_text(json.dumps({"blast_annot_batch": DC}))
    # reference: dc-megablast from scratch
    bench(sim_dir, tmp_path / "w_ref", tmp_path / "e_ref", "--steps_json", dc)
    # lead: default settings, snapshot after insert_finder_batch
    snap = tmp_path / "snap"
    bench(sim_dir, tmp_path / "w_lead", tmp_path / "e_lead", "--snapshot_after", "insert_finder_batch", snap)
    assert (snap / rb.STAMP).exists() and not (snap / "data" / "insertions").exists()
    # follower: copy of the snapshot, only the annotation steps with dc-megablast
    rb.copy_work(str(snap), str(tmp_path / "w_fol"))
    bench(sim_dir, tmp_path / "w_fol", tmp_path / "e_fol", "--steps_json", dc, "--from_step", "blast_annot_batch")
    assert results(tmp_path / "e_fol") == results(tmp_path / "e_ref")
    assert results(tmp_path / "e_lead") != results(tmp_path / "e_ref")   # the setting matters

    log = run_log(tmp_path / "e_fol")
    assert log["reused_steps"] == ["filter", "map", "insert_finder_batch"]
    assert log["from_step"] == "blast_annot_batch"
    steps = {s["step"].split(":")[0]: s for s in log["steps"]}
    assert set(steps) == set(rb.STEP_ORDER)
    assert steps["map"]["reused"] and "reused" not in steps["blast_annot_batch"]
    lead = {s["step"]: s for s in run_log(tmp_path / "e_lead")["steps"]}
    assert steps["map"]["seconds"] == lead["map"]["seconds"]
    assert all(s["cpu_seconds"] >= 0 for s in log["steps"])
    assert steps["map"]["cpu_seconds"] > 0 and steps["blast_annot_batch"]["cpu_seconds"] > 0

    # rerunning the annotation in the lead's own work dir: later outputs are
    # discarded first, so the result is again identical to a fresh run
    bench(sim_dir, tmp_path / "w_lead", tmp_path / "e_lead2", "--steps_json", dc, "--from_step",
          "blast_annot_batch")
    assert results(tmp_path / "e_lead2") == results(tmp_path / "e_ref")
    # --reuse_work alone resumes after the last step: nothing runs, evaluation is the same
    bench(sim_dir, tmp_path / "w_lead", tmp_path / "e_lead3", "--steps_json", dc, "--reuse_work")
    assert results(tmp_path / "e_lead3") == results(tmp_path / "e_ref")
    assert run_log(tmp_path / "e_lead3")["reused_steps"] == rb.STEP_ORDER


@requires_tools
def test_reuse_refused_on_mismatch(sim_dir, tmp_path):
    snap = tmp_path / "snap"
    bench(sim_dir, tmp_path / "w", tmp_path / "e", "--steps", "filter", "map", "insert_finder_batch",
          "--snapshot_after", "insert_finder_batch", snap)

    def refused(work, *extra, sim=sim_dir):
        p = bench(sim, work, tmp_path / "e_x", "--from_step", "blast_annot_batch", *extra, check=False)
        assert p.returncode == rb.REUSE_MISMATCH_RC, p.stdout[-2000:] + p.stderr[-2000:]
        assert "refused" in p.stderr
        return p

    refused(tmp_path / "nowhere")                                      # no stamp
    rb.copy_work(str(snap), str(tmp_path / "w1"))
    other = tmp_path / "other sim"
    other.mkdir()
    for f in os.listdir(sim_dir):
        if os.path.isfile(os.path.join(sim_dir, f)):
            (other / f).write_bytes(open(os.path.join(sim_dir, f), "rb").read())
    with open(other / "truth.tsv", "a") as fh:
        fh.write("\n")
    refused(tmp_path / "w1", sim=other)                                # different simulation
    mp = tmp_path / "map.json"
    mp.write_text(json.dumps({"map": rb.PROFILES["dev"]["map"] + " "}))
    refused(tmp_path / "w1", "--steps_json", mp)                       # different upstream command
    refused(tmp_path / "w1", "--min_len", "1500")                      # different run option
    bam = next((tmp_path / "w1" / "data").glob("barcode01/03_sorted_mapped_to_plasmid.bam"))
    with open(bam, "ab") as fh:
        fh.write(b"x")
    refused(tmp_path / "w1")                                           # changed upstream output
    # the refused work dir was left alone
    assert bam.exists() and (tmp_path / "w1" / rb.STAMP).exists()
