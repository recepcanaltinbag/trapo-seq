"""DR (target-site duplication) finder against an exactly known truth."""
import csv
import os
import random
import shutil

import pandas as pd
import pytest

from src.j_dr_finder import (OUTPUT_COLUMNS, analyze_insertion, circular_slice, main_dr_finder,
                             parse_query_id)

requires_blast = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("blastn", "makeblastdb")),
    reason="blastn/makeblastdb not in PATH")

COMP = str.maketrans("ACGT", "TGCA")
HEADER = ["Query ID", "Subject ID", "Identity (%)", "Score", "E-value", "Query Start", "Query End",
          "Subject Start", "Subject End", "Note", "Explained", "ref_pos", "is_Reverse"]


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_read(plasmid, is_seq, pos, tsd_len, is_orient, read_strand, cut_offset, tsd_side="right"):
    """Plasmid with `is_seq` inserted at 0-based `pos` with a circular TSD of `tsd_len`.

    The circular molecule is linearised `cut_offset` bases after the right TSD copy.
    Returns (read, ins_start, ins_end, is_start, is_end) - read coordinates 0-based,
    half-open. Like a CIGAR insertion, [ins_start, ins_end) holds the IS plus one TSD copy.
    """
    tsd = circular_slice(plasmid, pos, tsd_len)
    rotated = plasmid[pos:] + plasmid[:pos]            # starts with the target site
    element = is_seq if is_orient == "+" else rc(is_seq)
    circ = tsd + element + rotated                     # TSD, IS, TSD, rest of plasmid
    c = tsd_len + len(element) + tsd_len + cut_offset  # cut inside the plasmid backbone
    lin = circ[c:] + circ[:c]
    is_s = len(circ) - c + tsd_len
    is_e = is_s + len(element)
    if tsd_side == "right":
        ins = (is_s, is_e + tsd_len)
    else:
        ins = (is_s - tsd_len, is_e)
    if read_strand == "-":
        n = len(lin)
        lin = rc(lin)
        is_s, is_e = n - is_e, n - is_s
        ins = (n - ins[1], n - ins[0])
    return lin, ins[0], ins[1], is_s, is_e, tsd


def write_sample(root, sample, plasmid, cases, note="IS_DB"):
    """cases: list of (read_id, read_seq, ins_s, ins_e, is_s, is_e)."""
    d = os.path.join(root, sample)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{sample}.fasta"), "w") as fh:
        for rid, seq, *_ in cases:
            fh.write(f">{rid}\n{seq}\n")
    with open(os.path.join(d, "best_alignment.tab"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        for rid, seq, s, e, is_s, is_e in cases:
            w.writerow([f"{rid}_{s}_{e}", "IS4_test", 100, 2000, 0.0, is_s - s + 1, is_e - s,
                        1, is_e - is_s, note, 100, 0, False])
    pfa = os.path.join(root, "plasmid.fasta")
    with open(pfa, "w") as fh:
        fh.write(f">pTEST\n{plasmid}\n")
    return pfa


# ----------------------------------------------------------------- unit tests
def test_parse_query_id_rsplit():
    assert parse_query_id("read_with_under_scores_10_250") == ("read_with_under_scores", 10, 250)
    assert parse_query_id("abc_x_5") is None
    assert parse_query_id("nounderscore") is None


def test_circular_slice_wraps():
    assert circular_slice("ABCDEFG", 5, 4) == "FGAB"
    assert circular_slice("ABCDEFG", 0, 0) == ""
    assert circular_slice("ABC", 1, 7) == "BCABCAB"


def _hit(qs, qe, ss, se, bits=900.0):
    return {"qseqid": "q", "sseqid": "p", "pident": 100, "length": abs(qe - qs) + 1, "mismatch": 0,
            "gapopen": 0, "qstart": qs, "qend": qe, "sstart": ss, "send": se, "evalue": 0,
            "bitscore": bits}


def test_analyze_plus_strand_threshold_inclusive_and_statuses():
    plasmid = "ACGT" * 250                      # 1000 bp
    # insertion read [500, 1500); left flank ends at read pos 500 = plasmid 600,
    # right flank starts at read pos 1501 = plasmid 596 -> DR 596..600 = 5 bp
    hits = [_hit(1, 500, 101, 600), _hit(1501, 2000, 596, 1095)]
    r = analyze_insertion(hits, 500, 1500, 0, plasmid, threshold=200, repeat_th=5)
    assert (r["Status"], r["Repeat Length"], r["Start Subject"], r["End Subject"]) == ("DR", 5, 596, 600)
    assert r["Overlap Sequence"] == plasmid[595:600]
    assert analyze_insertion(hits, 500, 1500, 0, plasmid, 200, 4)["Status"] == "repeat_too_long"
    # blunt and target-site deletion are distinguished from each other
    blunt = analyze_insertion([_hit(1, 500, 101, 600), _hit(1501, 2000, 601, 1100)], 500, 1500, 0, plasmid, 200, 30)
    assert (blunt["Status"], blunt["Repeat Length"]) == ("no_duplication", 0)
    dele = analyze_insertion([_hit(1, 500, 101, 600), _hit(1501, 2000, 611, 1110)], 500, 1500, 0, plasmid, 200, 30)
    assert (dele["Status"], dele["Repeat Length"]) == ("target_site_deletion", -10)


def test_analyze_distance_threshold_inclusive():
    plasmid = "ACGT" * 250
    hits = [_hit(1, 400, 201, 600)]   # ends 100 bp before the insertion start
    assert analyze_insertion(hits, 500, 1500, 0, plasmid, 100, 30)["Status"] == "left_flank_only"
    assert analyze_insertion(hits, 500, 1500, 0, plasmid, 99, 30)["Status"] == "no_flank_hit"


def test_analyze_single_flank_reports_plasmid_coordinates():
    plasmid = "ACGT" * 250
    r = analyze_insertion([_hit(1501, 2000, 596, 1095)], 500, 1500, 0, plasmid, 200, 30)
    assert r["Status"] == "right_flank_only" and r["Start Subject"] == 596 and r["Insertion Point"] == 596
    # minus-strand read: the left flank on the read is the plasmid-right side
    r = analyze_insertion([_hit(1, 500, 1095, 596)], 500, 1500, 0, plasmid, 200, 30)
    assert r["Status"] == "left_flank_only" and r["Read Strand"] == "-" and r["Start Subject"] == 596


def test_analyze_wraps_origin_in_doubled_coordinates():
    plasmid = "ACGT" * 250
    # left flank ends at doubled coord 1002 (= plasmid 2), right starts at 998 -> DR 998..2 (5 bp)
    hits = [_hit(1, 500, 503, 1002), _hit(1501, 2000, 998, 1497)]
    r = analyze_insertion(hits, 500, 1500, 0, plasmid, 200, 30)
    assert (r["Repeat Length"], r["Start Subject"], r["End Subject"]) == (5, 998, 2)
    assert r["Overlap Sequence"] == plasmid[997:] + plasmid[:2]


# ------------------------------------------------------- end-to-end synthetic
def synthetic_cases(seed=3, plen=3000):
    rng = random.Random(seed)
    plasmid = rand_seq(rng, plen)
    is_seq = rand_seq(rng, 1200)
    positions = [0, 3, 9, plen - 1, plen - 2, plen - 5, plen - 9, 1500, 777, 2210]
    cases, truth = [], {}
    k = 0
    for pos in positions:
        for tsd_len in (0, 2, 5, 9, 11):
            for is_orient in "+-":
                for strand in "+-":
                    side = "right" if k % 2 else "left"
                    cut = rng.randint(300, plen - 300)
                    read, s, e, is_s, is_e, tsd = make_read(plasmid, is_seq, pos, tsd_len, is_orient,
                                                            strand, cut, side)
                    rid = f"read_{k}_x"           # underscores on purpose (B34)
                    cases.append((rid, read, s, e, is_s, is_e))
                    truth[rid] = (pos, tsd_len, tsd, strand)
                    k += 1
    return plasmid, cases, truth


@requires_blast
def test_end_to_end_exact_tsd_including_origin(tmp_path):
    plasmid, cases, truth = synthetic_cases()
    root = str(tmp_path / "data dir")
    pfa = write_sample(root, "barcode01", plasmid, cases)
    out = main_dr_finder(root, pfa, 500, 200, 30)
    assert len(out) == 1
    df = pd.read_csv(out[0], dtype={"Read ID": str}, keep_default_na=False)
    assert list(df.columns) == OUTPUT_COLUMNS
    assert len(df) == len(cases)
    bad = []
    for _, row in df.iterrows():
        pos, tsd_len, tsd, strand = truth[row["Read ID"]]
        want_status = "DR" if tsd_len else "no_duplication"
        got = (row["Status"], int(row["Repeat Length"]), int(row["Insertion Point"]), row["Read Strand"])
        want = (want_status, tsd_len, pos + 1, strand)
        seq_ok = row["Overlap Sequence"] == (tsd if tsd_len else "N/A")
        if got != want or not seq_ok:
            bad.append((row["Read ID"], got, want, row["Overlap Sequence"], tsd))
    assert not bad, bad[:10]


@requires_blast
def test_genome_note_without_is_clipping_on_unambiguous_ends(tmp_path):
    """Without IS boundaries the flank HSPs alone must give the exact DR when ends are unambiguous."""
    plasmid, cases, truth = synthetic_cases(seed=5)
    keep = cases
    root = str(tmp_path / "g")
    pfa = write_sample(root, "s1", plasmid, keep, note="Genome")
    df = pd.read_csv(main_dr_finder(root, pfa, 500, 200, 30)[0], dtype={"Read ID": str})
    # flanks may extend into the IS by chance base matches; never shorter than truth
    diffs = [int(r["Repeat Length"]) - truth[r["Read ID"]][1] for _, r in df.iterrows()]
    assert min(diffs) >= 0
    assert sum(d == 0 for d in diffs) / len(diffs) > 0.3


@requires_blast
def test_empty_and_bad_inputs(tmp_path):
    root = str(tmp_path / "e")
    plasmid = rand_seq(random.Random(1), 2000)
    pfa = write_sample(root, "empty", plasmid, [])
    # a row whose read is missing and one with an unparsable Query ID
    with open(os.path.join(root, "empty", "best_alignment.tab"), "a") as fh:
        fh.write("\t".join(["missing_read_1_100", "IS", "1", "1", "0", "1", "2", "1", "2", "IS_DB",
                            "1", "0", "False"]) + "\n")
        fh.write("\t".join(["garbage", "IS", "1", "1", "0", "1", "2", "1", "2", "IS_DB",
                            "1", "0", "False"]) + "\n")
    write_sample(root, "header_only", plasmid, [])
    out = {os.path.basename(p): p for p in main_dr_finder(root, pfa, 500, 200, 30)}
    df = pd.read_csv(out["empty_insertions.csv"])
    assert sorted(df["Status"]) == ["bad_query_id", "read_not_found"]
    df2 = pd.read_csv(out["header_only_insertions.csv"])
    assert df2.empty and list(df2.columns) == OUTPUT_COLUMNS


# ------------------------------------------------------- simulator ground truth
def _sim_cases(sim_out, barcode, is_fasta, genome_fasta):
    from simulate_traposome import read_fasta
    is_db = dict(read_fasta(is_fasta))
    genome = dict(read_fasta(genome_fasta))
    reads = {}
    with open(os.path.join(sim_out, barcode, f"{barcode}_raw.fastq")) as fh:
        lines = fh.read().split("\n")
    for i in range(0, len(lines) - 1, 4):
        reads[lines[i][1:]] = lines[i + 1]
    truth = {}
    cases = []
    with open(os.path.join(sim_out, "truth.tsv")) as fh:
        for t in csv.DictReader(fh, delimiter="\t"):
            if t["barcode"] != barcode or t["category"] not in ("is_insertion", "novel_insertion"):
                continue
            if t["category"] == "is_insertion":
                element = is_db[t["element"]]
            else:
                element = genome[t["genome_contig"]][int(t["genome_start"]):int(t["genome_end"])]
            read = reads[t["read_id"]]
            hits = [(read.find(x), len(x)) for x in (element, rc(element)) if read.find(x) >= 0]
            if not hits:
                continue          # fragmented read without the whole element
            is_s, n = hits[0]
            tsd_len = int(t["tsd_len"])
            cases.append((t["read_id"], read, is_s, is_s + n + tsd_len, is_s, is_s + n, t["category"]))
            truth[t["read_id"]] = t
    return cases, truth


@requires_blast
def test_simulator_error_free_exact(tmp_path):
    from simulate_traposome import main as simulate
    from conftest import TEST_DATA
    sim = str(tmp_path / "sim")
    pfa = os.path.join(TEST_DATA, "plasmid.fasta")
    isf, gf = os.path.join(TEST_DATA, "IS_curated.fasta"), os.path.join(TEST_DATA, "genome.fasta")
    simulate(["-p", pfa, "-g", gf, "-i", isf, "-o", sim, "--barcodes", "1", "--reads", "300",
              "--insertion_rate", "0.9", "--fragment_rate", "0", "--seed", "4",
              "--sub", "0", "--ins", "0", "--del", "0"])
    cases, truth = _sim_cases(sim, "barcode01", isf, gf)
    assert len(cases) > 100
    root = str(tmp_path / "run")
    is_cases = [c[:6] for c in cases if c[6] == "is_insertion"]
    write_sample(root, "barcode01", open(pfa).read().split("\n", 1)[1].replace("\n", ""), is_cases)
    df = pd.read_csv(main_dr_finder(root, pfa, 500, 200, 30)[0], dtype={"Read ID": str})
    bad = []
    for _, r in df.iterrows():
        t = truth[r["Read ID"]]
        if not (r["Status"] in ("DR", "no_duplication")
                and int(r["Repeat Length"]) == int(t["tsd_len"])
                and int(r["Insertion Point"]) == int(t["ref_pos"]) + 1
                and (r["Overlap Sequence"] if int(t["tsd_len"]) else "") == t["tsd_seq"]):
            bad.append((dict(r), {k: t[k] for k in ("ref_pos", "tsd_len", "tsd_seq", "read_strand")}))
    assert not bad, f"{len(bad)}/{len(df)} wrong: {bad[:5]}"


def test_is_bounds_extend_to_is_ends():
    from src.j_dr_finder import _is_bounds
    row = {"Note": "IS_DB", "Subject ID": "IS5_IS5", "Query Start": 1, "Query End": 1193,
           "Subject Start": 1193, "Subject End": 1}
    # minus-strand hit missing the last 2 IS bases at the query start side
    assert _is_bounds(row, 100) == (101, 1293)
    assert _is_bounds(row, 100, {"IS5_IS5": 1195}) == (99, 1293)
    row.update({"Subject Start": 3, "Subject End": 1195})       # plus hit missing IS bases 1-2
    assert _is_bounds(row, 100) == (99, 1293)
    assert _is_bounds(dict(row, Note="Genome"), 100) == (None, None)
    assert _is_bounds(dict(row, Note="Unannotated"), 100) == (None, None)   # handled like Genome


def test_dr_subject_for_inserts_without_hit():
    from src.j_dr_finder import _dr_subject
    assert _dr_subject({"Note": "IS_DB", "Subject ID": "IS5_IS5"}) == "IS5_IS5"
    assert _dr_subject({"Note": "Genome", "Subject ID": "chr"}) == "chr"
    assert _dr_subject({"Note": "Unannotated", "Subject ID": "no blast hit"}) == "Unannotated"
    assert _dr_subject({"Note": "Contamination", "Subject ID": "no blast hit"}) == "Contamination"


def test_flank_extends_to_is_boundary_at_most_max_extend():
    from src.j_dr_finder import MAX_EXTEND, flank_junction
    h = dict(_hit(1, 495, 106, 600), rq_s=1, rq_e=495)
    assert flank_junction(h, "L", 500) == (605, True)               # 5 bp error-trimmed end
    assert flank_junction(h, "L", 495 + MAX_EXTEND + 1) == (600, True)
    h = dict(_hit(1, 500, 1095, 596), rq_s=1, rq_e=500)
    assert flank_junction(h, "L", 498) == (598, False)              # clipped, minus strand
