"""Unit tests for simulate/evaluate.py on hand-made truth / prediction tables."""
import json
import os

import pytest

import evaluate as ev

TRUTH_COLS = ["barcode", "read_id", "category", "clone_id", "element", "element_family", "ref_pos",
              "orientation", "tsd_len", "tsd_seq", "insert_len", "genome_contig", "genome_start",
              "genome_end", "read_strand", "fragmented", "read_len"]
CLONE_COLS = ["barcode", "clone_id", "category", "element", "element_family", "ref_pos", "orientation",
              "tsd_len", "tsd_seq", "insert_len", "genome_contig", "genome_start", "genome_end", "weight", "n_reads"]
ANNOT_HEADER = ("Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\tSubject Start\t"
                "Subject End\tNote\tExplained\tref_pos\tis_Reverse\n")

BC = "barcode01"
# read_id, category, clone, element, ref_pos, tsd_len, tsd_seq, genome_start, genome_end, fragmented
READS = [
    ("barcode01_r000001", "is_insertion", "c0", "IS4_ISVsa5", 1000, 11, "ACGTACGTACG", "", "", False),
    ("barcode01_r000002", "is_insertion", "c0", "IS4_ISVsa5", 1000, 11, "ACGTACGTACG", "", "", False),
    ("barcode01_r000003", "is_insertion", "c1", "IS5_IS5", 3000, 4, "GGAT", "", "", False),
    ("barcode01_r000004", "novel_insertion", "c2", "novel", 5000, 5, "AAAAC", 100000, 102000, False),
    ("barcode01_r000005", "is_insertion", "c1", "IS5_IS5", 3000, 4, "GGAT", "", "", True),  # missed, fragmented
    ("barcode01_r000006", "empty_plasmid", "", "", "", "", "", "", "", False),
    ("barcode01_r000007", "empty_plasmid", "", "", "", "", "", "", "", False),  # false call
    ("barcode01_r000008", "contamination", "", "", "", "", "", 5, 5000, False),
    ("barcode01_r000009", "deletion", "", "", 200, "", "", "", "", False),
]
CLONES = [
    ("c0", "is_insertion", "IS4_ISVsa5", 1000, 11, 2),
    ("c1", "is_insertion", "IS5_IS5", 3000, 4, 2),
    ("c2", "novel_insertion", "novel", 5000, 5, 1),
    ("c3", "is_insertion", "IS5_IS5", 6000, 4, 0),  # no reads -> not a true event
]


def write_tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


@pytest.fixture
def sim(tmp_path):
    sim = tmp_path / "sim dir"
    sim.mkdir()
    rows = []
    for rid, cat, cid, el, pos, tsd, seq, gs, ge, frag in READS:
        rows.append(dict(barcode=BC, read_id=rid, category=cat, clone_id=f"{BC}_{cid}" if cid else "",
                         element=el, ref_pos=pos, tsd_len=tsd, tsd_seq=seq, genome_contig="chr" if gs != "" else "",
                         genome_start=gs, genome_end=ge, fragmented=frag, read_len=7000))
    write_tsv(sim / "truth.tsv", TRUTH_COLS, rows)
    write_tsv(sim / "clones.tsv", CLONE_COLS,
              [dict(barcode=BC, clone_id=f"{BC}_{c}", category=cat, element=el, ref_pos=p, tsd_len=t, n_reads=n)
               for c, cat, el, p, t, n in CLONES])
    (sim / "simulation.json").write_text(json.dumps({"plasmid_len": 7000}))
    return sim


@pytest.fixture
def results(tmp_path):
    res = tmp_path / "results"
    bdir = res / BC
    bdir.mkdir(parents=True)
    # insert finder: read IDs rewritten with '-', duplicate rows for r1 (IN + SC),
    # r2 placed at the right end of the TSD, r7 is a false call on an empty plasmid
    (bdir / "insertion_from_bam.tab").write_text(
        "barcode01-r000001\t100\t1400\t1300\t1000\t20\tIN\tFalse\n"
        "barcode01-r000001\t100\t700\t600\t1500\t20\tSC\tFalse\n"
        "barcode01-r000002\t200\t1500\t1300\t1011\t20\tIN\tTrue\n"
        "barcode01-r000003\t300\t1500\t1200\t3030\t20\tIN\tFalse\n"
        "barcode01-r000004\t300\t2300\t2000\t5002\t20\tSC\tFalse\n"
        "barcode01-r000007\t300\t1000\t700\t6500\t20\tSC\tFalse\n")
    (bdir / "best_alignment.tab").write_text(
        ANNOT_HEADER +
        "barcode01-r000001_100_1400\tIS4_ISVsa5\t97\t1200\t0.0\t1\t1300\t1\t1300\tIS_DB\t99\t1000\tFalse\n"
        "barcode01-r000001_100_1400\tchr\t97\t800\t0.0\t1\t1300\t1\t1300\tGenome\t99\t1000\tFalse\n"
        "barcode01-r000002_200_1500\tIS4_IS186B\t97\t1200\t0.0\t1\t1300\t1\t1300\tIS_DB\t99\t1011\tTrue\n"
        "barcode01-r000003_300_1500\tN/A\tN/A\tN/A\tN/A\t0\t0\t0\t0\tContamination\t0\t3030\tFalse\n"
        "barcode01-r000004_300_2300\tchr\t97\t1900\t0.0\t1\t2000\t101950\t100001\tGenome\t99\t5002\tFalse\n"
        "barcode01-r000007_300_1000\tN/A\tN/A\tN/A\tN/A\t0\t0\t0\t0\tContamination\t0\t6500\tFalse\n")
    ins = res / "insertions"
    ins.mkdir()
    (ins / f"{BC}_insertions.csv").write_text(
        "Read ID,Subject ID,Query Start CSV,Query End CSV,Start Subject,End Subject,Insertion Point,Repeat Length,Overlap Sequence\n"
        "barcode01-r000001,IS4_ISVsa5,100,1400,1001,1011,1001,11,ACGTACGTACG\n"
        "barcode01-r000002,IS4_IS186B,200,1500,1001,1012,1001,12,ACGTACGTACGT\n"
        "barcode01-r000003,N/A,300,1500,3031,-1,3031,0,N/A\n"
        "barcode01-r000004,chr,300,2300,5001,5005,5001,5,AAAAC\n")
    return res


def run(sim, res, **kw):
    report, per_read, sites = ev.evaluate(str(sim), str(res), **kw)
    return report, {r["read_id"]: r for r in per_read}, sites


def test_helpers():
    assert ev.circ_dist(10, 6990, 7000) == 20
    assert ev.circ_dist(10, 30) == 20
    assert ev.interval_dist(1005, 1000, 1011) == 0
    assert ev.interval_dist(990, 1000, 1011) == 10
    assert ev.interval_dist(1015, 1000, 1011) == 4
    assert ev.split_query_id("barcode01_r000001_100_1400") == ("barcode01_r000001", 100, 1400)
    assert ev.split_query_id("weird") == ("weird", None, None)
    assert ev.f1(1.0, 0.5) == pytest.approx(2 / 3)
    assert ev.f1(None, 0.5) is None
    assert ev.ratio(1, 0) is None


def test_repeat_extension():
    plasmid = "TTTTGACGTTTT"   # TSD = plasmid[4:7] = "GAC"
    # insert ends with 'T' == plasmid[3] -> one-base left extension; starts with 'C' != plasmid[7]
    left, right, ext = ev.repeat_extension(plasmid, 4, 3, "CCCCT")
    assert (left, right) == (1, 0)
    assert ext == "TGAC"
    left, right, ext = ev.repeat_extension(plasmid, 4, 3, "AAAAA")
    assert (left, right, ext) == (0, 0, "GAC")
    # insert starting with 'G' == plasmid[7] extends the repeat to the right
    assert ev.repeat_extension(plasmid, 4, 3, "GAAAA") == (0, 1, "GACG")


def test_read_level_detection(sim, results):
    report, per_read, _ = run(sim, results, min_read_len=2000)
    d = report["per_barcode"][BC]["detection_finder"]
    assert (d["tp"], d["fp"], d["fn"]) == (4, 1, 1)
    assert d["precision"] == pytest.approx(0.8)
    assert d["recall"] == pytest.approx(0.8)
    assert d["recall_unfragmented_eligible"] == pytest.approx(1.0)   # the missed read is fragmented
    assert d["false_calls_empty_plasmid"] == 1
    assert d["false_call_rate_empty_plasmid"] == pytest.approx(0.5)
    assert d["false_calls_contamination"] == 0
    # annotation stage: Contamination rows are not biological calls
    a = report["per_barcode"][BC]["detection_annot"]
    assert (a["tp"], a["fp"], a["fn"]) == (3, 0, 2)
    # overall equals the single barcode here
    assert report["overall"]["detection_finder"]["tp"] == 4


def test_dedupe_and_position(sim, results):
    report, per_read, _ = run(sim, results)
    r1 = per_read["barcode01_r000001"]
    assert r1["n_finder_rows"] == 2 and r1["pred_ref_pos"] == 1000   # longest insertion wins
    assert r1["pos_err"] == 0
    assert per_read["barcode01_r000002"]["pos_err"] == 0               # right end of TSD window
    assert per_read["barcode01_r000002"]["pos_err_raw"] == 11
    assert per_read["barcode01_r000003"]["pos_err"] == 26              # 3030 vs [3000, 3004]
    pos = report["per_barcode"][BC]["position"]["error"]
    assert pos["n"] == 4 and pos["median"] == 0
    assert pos["pct_within_20bp"] == pytest.approx(75.0)
    m = report["per_barcode"][BC]["multiplicity_finder"]
    assert m["pct_reads_multiple_rows"] == pytest.approx(20.0)
    assert m["unmatched_read_ids"] == 0
    d = report["per_barcode"][BC]["detection_finder"]
    assert d["tp_correct_position"] == 3


def test_classification(sim, results):
    report, per_read, _ = run(sim, results)
    assert per_read["barcode01_r000001"]["annot_note"] == "IS_DB"      # IS_DB beats Genome
    assert per_read["barcode01_r000001"]["class_ok"] is True
    assert per_read["barcode01_r000002"]["class_ok"] is False          # wrong IS
    assert per_read["barcode01_r000003"]["class_ok"] is False          # labelled Contamination
    assert per_read["barcode01_r000004"]["class_ok"] is True           # Genome hit overlaps truth (minus strand)
    c = report["per_barcode"][BC]["classification"]
    assert c["is_insertion"]["n_annotated"] == 3
    assert c["is_insertion"]["accuracy"] == pytest.approx(1 / 3)
    assert c["is_insertion"]["pct_note_Contamination"] == pytest.approx(100 / 3)
    assert c["is_insertion"]["top_confusions"] == ["IS4_ISVsa5->IS4_IS186B:1"]
    assert c["novel_insertion"]["accuracy"] == 1.0
    assert c["contamination_note_by_truth"] == {"is_insertion": 1, "empty_plasmid": 1}
    assert c["contamination_note_true_insertions"] == 1
    assert c["contamination_reads_not_called"] == 1.0


def test_tsd(sim, results):
    report, per_read, _ = run(sim, results)
    t = report["per_barcode"][BC]["tsd"]
    assert t["reads_with_dr_rows"] == 4
    assert t["length_accuracy"] == pytest.approx(0.5)    # r1 and r4 exact; r2 too long; r3 none
    assert t["sequence_accuracy"] == pytest.approx(0.5)
    assert t["pct_tsd_reported_nonzero"] == pytest.approx(75.0)
    assert t["length_accuracy_when_reported"] == pytest.approx(2 / 3)
    assert t["insertion_point_error"]["median"] == 0


def test_sites(sim, results):
    report, _, sites = run(sim, results, tol=20)
    s = report["per_barcode"][BC]["sites"]
    assert s["true_events"] == 3                       # c3 has no reads
    # sites are split by label: IS4_ISVsa5@1000, IS4_IS186B@1011, Contamination@3030,
    # Genome@5002, Contamination@6500
    assert s["predicted_sites"] == 5
    assert s["matched"] == 2                           # 1011 duplicates c0; 3030 is 26 bp from c1
    assert s["spurious_sites"] == 3
    assert s["event_recall"] == pytest.approx(2 / 3)
    assert s["missed_events"] == ["barcode01_c1"]
    s30 = run(sim, results, tol=30)[0]["per_barcode"][BC]["sites"]
    assert s30["matched"] == 3 and s30["spurious_sites"] == 2
    # same label, positions chained within --site_link -> a single site
    p = results / BC / "best_alignment.tab"
    p.write_text(p.read_text().replace("IS4_IS186B", "IS4_ISVsa5"))
    assert run(sim, results, site_link=20)[0]["per_barcode"][BC]["sites"]["predicted_sites"] == 4
    assert run(sim, results, site_link=5)[0]["per_barcode"][BC]["sites"]["predicted_sites"] == 5
    lab = {x["clone_id"]: x["label_ok"] for x in sites if x["clone_id"]}
    assert lab["barcode01_c2"] is True


def test_missing_files(sim, tmp_path):
    empty = tmp_path / "nothing"
    (empty / BC).mkdir(parents=True)
    report, per_read, sites = ev.evaluate(str(sim), str(empty))
    b = report["per_barcode"][BC]
    assert b["detection_finder"] == {"available": False}
    assert b["tsd"] == {"available": False}
    assert b["sites"] == {"available": False}
    assert report["overall"]["classification"] == {"available": False}
    out = tmp_path / "out"
    ev.write_outputs(report, per_read, sites, str(out))
    text = (out / "metrics.tsv").read_text()
    assert "available\tFalse" in text
    assert "\n".join(ev.summary_lines(report)).count("N/A") > 5


def test_cli_writes_outputs(sim, results, tmp_path):
    out = tmp_path / "eval out"
    ev.main(["-s", str(sim), "-r", str(results), "-o", str(out), "-q"])
    for f in ("metrics.json", "metrics.tsv", "per_read.tsv", "sites.tsv", "summary.tsv"):
        assert os.path.exists(out / f)
    data = json.loads((out / "metrics.json").read_text())
    assert data["overall"]["detection_finder"]["tp"] == 4


def test_headered_finder_table(sim, results):
    # A future insert finder may write a header row and extra columns.
    p = results / BC / "insertion_from_bam.tab"
    p.write_text("read_id\tquery_start\tquery_end\tinsertion_length\tref_pos\tquality\tinsert_type\tis_reverse\textra\n"
                 "barcode01_r000001\t100\t1400\t1300\t1003\t20\tIN\tFalse\tx\n")
    report, per_read, _ = run(sim, results)
    assert per_read["barcode01_r000001"]["pred_ref_pos"] == 1003
    assert report["per_barcode"][BC]["detection_finder"]["tp"] == 1
