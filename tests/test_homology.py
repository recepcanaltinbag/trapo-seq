"""Plasmid-genome shared segments (GitHub issue #5): detection, corrected
mapping statistics, flagging of homology artefacts, annotation/stats/evaluation
outputs and the simulator's --homology_len."""
import csv
import json
import os
import random
import shutil

import pandas as pd
import pytest

from conftest import TEST_DATA, requires_tools
from simulate_traposome import main as simulate, revcomp
from src import h_is_stats
from src import plasmid_genome_homology as hom
from src.g_blast_annot import HOMOLOGY_PREFIX, OUTPUT_COLUMNS, load_homology_flags

PLASMID = os.path.join(TEST_DATA, "plasmid.fasta")
GENOME = os.path.join(TEST_DATA, "genome.fasta")
IS_FASTA = os.path.join(TEST_DATA, "IS_curated.fasta")
SIM_ARGS = ["-p", PLASMID, "-g", GENOME, "-i", IS_FASTA]


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _reads(fastq):
    with open(fastq) as fh:
        lines = fh.read().split("\n")
    return {lines[i][1:]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 1, 4)}


def _seq(path):
    return "".join(hom.read_fasta(path)[0][1]).upper()


# ------------------------------------------------------------- units

def paf(qname, qstart, qend, tname, tstart, tend, matches, block, strand="+"):
    return "\t".join(map(str, [qname, 99999, qstart, qend, strand, tname, 10 ** 6, tstart, tend,
                               matches, block, 60]))


def test_parse_paf_hits_filters_dedupes_and_handles_origin():
    lengths = {"p": 1000}
    lines = [
        paf("p", 100, 400, "g", 5000, 5300, 300, 300),         # kept
        paf("p", 1100, 1400, "g", 5000, 5300, 300, 300),       # same copy found in the appended start
        paf("p", 900, 1150, "g", 9000, 9250, 250, 250),        # crosses the origin
        paf("p", 0, 150, "g", 9100, 9250, 150, 150),           # piece of the crossing hit
        paf("p", 500, 560, "g", 100, 160, 60, 60),             # too short
        paf("p", 600, 800, "g", 200, 400, 140, 200),           # identity 0.7
        paf("p", 100, 400, "g", 70000, 70300, 290, 300, "-"),  # second genome copy
    ]
    segs = hom.parse_paf_hits(lines, lengths, min_len=100, min_identity=0.8)
    got = [(s["plasmid_start"], s["plasmid_end"], s["genome_start"], s["strand"]) for s in segs]
    assert got == [(100, 400, 5000, "+"), (100, 400, 70000, "-"), (900, 1150, 9000, "+")]
    assert segs[1]["identity"] == pytest.approx(0.9667, abs=1e-4)


def test_segments_at_is_circular():
    segs = [{"plasmid": "p", "plasmid_len": 1000, "plasmid_start": 900, "plasmid_end": 1150}]
    assert hom.segments_at("p", 950, segs) == [0]
    assert hom.segments_at("p", 100, segs) == [0]          # past the origin
    assert hom.segments_at("p", 190, segs) == [0]          # 40 bp after the end (tolerance 50)
    assert hom.segments_at("p", 300, segs) == []
    assert hom.segments_at("p", 860, segs) == [0]          # 40 bp before the start
    assert hom.segments_at("q", 950, segs) == []


def test_junctions_and_genome_contiguity():
    # read: plasmid 0-1000 | insert 1000-3000 | plasmid 3000-4000
    pieces = [(0, 1000, 10, 1010, "p"), (3000, 4000, 1010, 2010, "p")]
    j = hom.insertion_junctions(1000, 3000, pieces)
    assert [(x[0], x[1], x[2]) for x in j] == [("left", 1000, 1010), ("right", 3000, 1010)]
    # an end clip has one flank
    assert len(hom.insertion_junctions(3000, 4500, [(0, 3000, 0, 3000, "p")])) == 1
    genome = [(500, 4000)]
    assert hom.genome_contiguous("left", 1000, 1000, 2000, genome)
    assert hom.genome_contiguous("right", 3000, 1000, 2000, genome)
    assert not hom.genome_contiguous("left", 1000, 1000, 2000, [(1000, 3000)])   # insert only
    assert not hom.genome_contiguous("left", 1000, 1000, 2000, [(0, 1000)])      # flank only


def test_classify_genome_hits():
    plasmid = [(0, 3000)]
    assert hom.classify_genome_hits(plasmid, [(1000, 2000)]) == "homology_only"     # inside the plasmid part
    assert hom.classify_genome_hits(plasmid, [(3000, 4500)]) == "genome"             # insert only
    assert hom.classify_genome_hits(plasmid, [(2000, 5000)]) == "genomic"            # runs across
    assert hom.classify_genome_hits(plasmid, [(2950, 4500)]) == "genome"             # 50 bp overlap only


def test_unaligned_stretches():
    assert hom.unaligned_stretches([(0, 1000), (1050, 2000)], 2500) == [(2000, 2500)]
    assert hom.unaligned_stretches([(300, 2000)], 2000) == [(0, 300)]


def test_messages_warn_only_for_long_segments():
    seg = {"plasmid": "p", "plasmid_start": 0, "plasmid_end": 230, "length": 230, "identity": 1.0}
    warnings, infos = hom.homology_messages([seg], warn_len=500)
    assert warnings == [] and "shorter than 500" in infos[0]
    warnings, _ = hom.homology_messages([dict(seg, plasmid_end=1500, length=1500)], warn_len=500)
    assert "1500 bp" in warnings[0]
    assert hom.homology_messages([], 500) == ([], ["no plasmid-genome shared segments found."])


def test_homology_table_round_trip(tmp_path):
    seg = {"plasmid": "p", "plasmid_len": 1000, "plasmid_start": 900, "plasmid_end": 1150,
           "genome_contig": "g", "genome_start": 5, "genome_end": 255, "strand": "-",
           "identity": 0.98, "length": 250}
    path = tmp_path / "t.tsv"
    hom.write_homology([seg], str(path))
    assert hom.read_homology(str(path)) == [seg]
    hom.write_homology([], str(path))
    assert hom.read_homology(str(path)) == []


def test_load_homology_flags(tmp_path):
    tab = tmp_path / "insertion_from_bam.tab"
    tab.write_text("r1\t0\t900\t900\t10\t20\tSC\tFalse\n")
    assert load_homology_flags(str(tab)) == set()                       # no flags file
    flags = tmp_path / "insertion_from_bam_homology_flags.tsv"
    assert hom.flags_path_for(str(tab)) == str(flags)
    hom.write_flags([{"Query ID": "r1_0_900", "read_id": "r1", "query_start": 0, "query_end": 900,
                      "insert_type": "SC", "junctions": 1, "junctions_in_shared_segment": 1,
                      "junctions_genome_contiguous": 1, "shared_segments": "0-1500",
                      "homology_flag": True},
                     {"Query ID": "r2_5_900", "read_id": "r2", "query_start": 5, "query_end": 900,
                      "insert_type": "SC", "junctions": 1, "junctions_in_shared_segment": 1,
                      "junctions_genome_contiguous": 0, "shared_segments": "0-1500",
                      "homology_flag": False}], str(flags))
    assert load_homology_flags(str(tab)) == {"r1_0_900"}
    assert load_homology_flags(str(tab), None) == set()
    os.utime(tab, (os.path.getmtime(flags) + 10,) * 2)                    # stale flags are ignored
    assert load_homology_flags(str(tab)) == set()
    with pytest.raises(FileNotFoundError):
        load_homology_flags(str(tab), str(tmp_path / "missing.tsv"))


def test_stats_report_homology_flagged(tmp_path):
    cols = OUTPUT_COLUMNS
    rows = [["q1", "IS5_IS5", 99, 900, 0, 1, 900, 1, 900, "IS_DB", 99, 10, False, "-", "-"],
            ["q2", "chr", 99, 900, 0, 1, 900, 1, 900, "Contamination", 99, 10, False,
             HOMOLOGY_PREFIX + "Genome", "-"]]
    os.makedirs(tmp_path / "bc1")
    pd.DataFrame(rows, columns=cols).to_csv(tmp_path / "bc1" / "best_alignment.tab", sep="\t", index=False)
    s = h_is_stats.summarize_condition(str(tmp_path / "bc1" / "best_alignment.tab"))
    assert s["homology_flagged"] == 1 and s["is_count"] == 1
    text = h_is_stats.format_condition("bc1", s)
    assert "#Homology_flagged=1" in text and "#Contamination=1" in text
    s["homology_flagged"] = 0
    assert "Homology" not in h_is_stats.format_condition("bc1", s)


def test_evaluate_homology_metrics():
    import evaluate as ev

    def r(cat, notes, unflagged, flagged, overlap):
        return {"category": cat, "annot_notes_all": notes, "annot_notes_unflagged_all": unflagged,
                "homology_flagged": flagged, "homology_overlap": overlap}
    rows = [r("contamination", "Contamination", "Genome", True, 800),        # caught
            r("contamination", "Genome", "Genome", False, 500),              # missed
            r("contamination", "Genome", "Genome", False, 0),                # not homology
            r("contamination", "", "", False, 0),
            r("is_insertion", "IS_DB", "IS_DB", False, 0),
            r("empty_plasmid", "", "", False, 0)]
    d = ev.homology_metrics(rows, True, 1, {"Transposition Rate": "0.50",
                                             "Transposition Rate (uncorrected)": "1.00"})
    assert d["false_calls_contamination_before"] == 3 and d["false_calls_contamination_after"] == 2
    assert d["false_calls_homology_before"] == 2 and d["false_calls_homology_after"] == 1
    assert d["flag_caught_homology"] == 1 and d["flag_sensitivity"] == 0.5
    assert d["flagged_true_insertions"] == 0 and d["recall_annot_after"] == 1.0
    assert d["transposition_rate_true"] == 0.5 and d["transposition_rate_uncorrected"] == 1.0
    assert ev.homology_metrics(rows, False, 1) == {"available": False}


# ------------------------------------------------------------ simulator

def test_simulator_homology_outputs(tmp_path):
    out_a, out_b = tmp_path / "a", tmp_path / "b"
    common = SIM_ARGS + ["--barcodes", "1", "--reads", "300", "--contamination_rate", "0.2", "--seed", "3"]
    simulate(common + ["-o", str(out_a)])
    simulate(common + ["-o", str(out_b), "--homology_len", "1500", "1500", "--n_homology", "2",
                       "--homology_contam_frac", "0.5"])
    plasmid = _seq(PLASMID)
    genome = dict(hom.read_fasta(str(out_b / "genome_with_homology.fasta")))
    original = dict(hom.read_fasta(GENOME))
    segs = _tsv(out_b / "homology.tsv")
    assert len(segs) == 2
    changed = 0
    for sg in segs:
        ps, pe, gs, ge = (int(sg[k]) for k in ("plasmid_start", "plasmid_end", "genome_start", "genome_end"))
        assert pe - ps == ge - gs == int(sg["length"]) == 1500
        copy = genome[sg["genome_contig"]][gs:ge]
        seg = (plasmid + plasmid)[ps:pe]
        assert copy == (seg if sg["genome_strand"] == "+" else revcomp(seg))
        changed += ge - gs
    g, o = next(iter(genome.values())), next(iter(original.values()))
    assert len(g) == len(o)                                          # replaced, not inserted
    assert sum(a != b for a, b in zip(g, o)) <= changed
    meta = json.loads((out_b / "simulation.json").read_text())
    assert meta["homology_len"] == [1500, 1500] and len(meta["homology_segments"]) == 2
    assert "homology_len" not in json.loads((out_a / "simulation.json").read_text())
    assert not (out_a / "homology.tsv").exists()

    # contamination reads: about half overlap a segment copy; all other reads are
    # identical to the run without shared segments (same main random stream)
    truth_a, truth_b = _tsv(out_a / "truth.tsv"), _tsv(out_b / "truth.tsv")
    assert [t["category"] for t in truth_a] == [t["category"] for t in truth_b]
    contam = [t for t in truth_b if t["category"] == "contamination"]
    hits = [t for t in contam if any(
        int(t["genome_start"]) < int(sg["genome_end"]) - 99 and int(t["genome_end"]) > int(sg["genome_start"]) + 99
        for sg in segs)]
    assert 0.25 * len(contam) <= len(hits) <= 0.8 * len(contam)
    ra = _reads(out_a / "barcode01" / "barcode01_raw.fastq")
    rb = _reads(out_b / "barcode01" / "barcode01_raw.fastq")
    for ta, tb in zip(truth_a, truth_b):
        if ta["category"] != "contamination":
            assert ra[ta["read_id"]] == rb[tb["read_id"]]
        else:
            assert len(ra[ta["read_id"]][0]) == len(rb[tb["read_id"]][0]) or ta["genome_contig"] != tb["genome_contig"]


def test_simulator_homology_validation(tmp_path):
    with pytest.raises(SystemExit):
        simulate(SIM_ARGS + ["-o", str(tmp_path), "--homology_len", "2000", "1000"])
    with pytest.raises(SystemExit):
        simulate(SIM_ARGS + ["-o", str(tmp_path), "--homology_len", "100", "200", "--homology_contam_frac", "2"])


# --------------------------------------------------------- with tools

@requires_tools
def test_detect_homology_finds_planted_segments(tmp_path):
    plasmid = _seq(PLASMID)
    rng = random.Random(1)
    rnd = "".join(rng.choice("ACGT") for _ in range(60000))
    n = len(plasmid)
    inner = plasmid[2000:3200]                          # 1200 bp, forward
    across = plasmid[n - 300:] + plasmid[:500]          # 800 bp crossing the origin, reverse
    genome = rnd[:10000] + inner + rnd[10000:40000] + revcomp(across) + rnd[40000:]
    gpath = tmp_path / "genome with space.fasta"
    gpath.write_text(">chr\n" + genome + "\n")
    segs = hom.detect_homology(PLASMID, str(gpath), min_len=100, min_identity=0.8)
    got = sorted((s["plasmid_start"], s["plasmid_end"], s["genome_start"], s["strand"]) for s in segs)
    assert len(got) == 2, segs
    (a0, a1, ag, astr), (b0, b1, bg, bstr) = got
    assert abs(a0 - 2000) <= 5 and abs(a1 - 3200) <= 5 and abs(ag - 10000) <= 5 and astr == "+"
    assert abs(b0 - (n - 300)) <= 5 and abs(b1 - (n + 500)) <= 5 and bstr == "-"


@requires_tools
def test_homology_end_to_end(tmp_path):
    """Contamination reads through a 2 kb shared segment: called by the insert
    finder, flagged, relabelled Contamination; true insertions untouched;
    08_mapping_stats.txt corrected."""
    from src.a_data_mapping import main_mapping
    from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch
    from src.g_blast_annot_batch import main_annot_batch

    sim = tmp_path / "sim"
    simulate(SIM_ARGS + ["-o", str(sim), "--barcodes", "1", "--reads", "200", "--seed", "4",
                         "--contamination_rate", "0.25", "--homology_len", "2000", "2000",
                         "--homology_contam_frac", "1", "--error_model", "ont",
                         "--ont_profile", "r10.4.1"])
    data = tmp_path / "data dir"
    (data / "barcode01").mkdir(parents=True)
    shutil.copy(sim / "barcode01" / "barcode01_raw.fastq", data / "barcode01" / "barcode01.fastq")
    genome = str(sim / "genome_with_homology.fasta")
    segs = main_mapping(str(data), PLASMID, genome)
    planted = _tsv(sim / "homology.tsv")[0]
    assert any(abs(s["plasmid_start"] - int(planted["plasmid_start"]) % s["plasmid_len"]) <= 10
               and s["length"] >= 1990 for s in segs)
    assert (data / hom.HOMOLOGY_NAME).is_file()

    st = {}
    for line in (data / "barcode01" / "08_mapping_stats.txt").read_text().splitlines():
        k, v = line.rsplit(":", 1)
        st[k.strip()] = v.strip()
    assert float(st["Transposition Rate"]) < float(st["Transposition Rate (uncorrected)"])
    assert int(st["Genomic reads reaching the plasmid through a shared segment (excluded)"]) > 0

    main_insert_finder_from_bam_batch(str(data), 500)
    flags = _tsv(data / "barcode01" / "insertion_from_bam_homology_flags.tsv")
    truth = {t["read_id"]: t for t in _tsv(sim / "truth.tsv")}
    flagged = {f["read_id"] for f in flags if f["homology_flag"] == "True"}
    assert flagged and all(truth[r]["category"] == "contamination" for r in flagged)

    main_annot_batch(str(data), genome, is_fasta=IS_FASTA)
    df = pd.read_csv(data / "barcode01" / "best_alignment.tab", sep="\t")
    assert list(df.columns) == OUTPUT_COLUMNS
    marked = df[df["Homology"].str.startswith(HOMOLOGY_PREFIX)]
    assert len(marked) and set(marked["Note"]) == {"Contamination"}
    called = df[df["Note"].isin(["IS_DB", "Genome", "Unannotated"])]
    reads = {q.rsplit("_", 2)[0] for q in called["Query ID"]}
    assert all(truth[r]["category"] in ("is_insertion", "novel_insertion") for r in reads)
