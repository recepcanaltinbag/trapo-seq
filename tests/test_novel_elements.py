"""Novel mobile element characterisation (src/l_novel_elements.py)."""
import contextlib
import csv
import io
import os
import random
import shutil

import pytest

from src.l_novel_elements import (boundary_candidates, characterise, cluster_intervals,
                                  find_dr, find_tir, longest_orf, main_novel_elements,
                                  parse_query_id, revcomp, robust_end, scan_known_tirs)

from conftest import TEST_DATA, requires_tools


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, positions):
    s = list(seq)
    for i in positions:
        s[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]]
    return "".join(s)


def element_in_context(seed, tir_len=25, mismatches=(), dr="CATCCA", length=1500):
    """flank + DR + [TIR ... revcomp(TIR)] + DR + flank; returns (seq, start, end, tir)."""
    rng = random.Random(seed)
    tir = rand_seq(rng, tir_len)
    right = mutate(revcomp(tir), mismatches)
    body = tir + rand_seq(rng, length - 2 * tir_len) + right
    left_flank, right_flank = rand_seq(rng, 300), rand_seq(rng, 300)
    seq = left_flank + dr + body + dr + right_flank
    start = len(left_flank) + len(dr)
    return seq, start, start + len(body), tir


# ------------------------------------------------------------------ find_tir

def test_find_tir_exact():
    seq, s, e, tir = element_in_context(1, tir_len=25)
    t = find_tir(seq, s, e)
    assert t is not None
    assert (t["start"], t["end"]) == (s, e)
    assert abs(t["left_len"] - 25) <= 2 and abs(t["right_len"] - 25) <= 2
    assert t["left_seq"].startswith(tir[:20])
    assert t["identity"] >= 0.9


def test_find_tir_with_mismatches():
    # 3 of 30 positions mismatched (10 %), outer 2 bp intact
    seq, s, e, _ = element_in_context(2, tir_len=30, mismatches=(6, 14, 22))
    t = find_tir(seq, s, e)
    assert t is not None
    assert (t["start"], t["end"]) == (s, e)
    assert 26 <= t["left_len"] <= 33
    assert 0.8 <= t["identity"] < 1.0


@pytest.mark.parametrize("ds,de", [(-4, 3), (5, -5), (2, 4), (-6, -2)])
def test_find_tir_from_offset_estimates(ds, de):
    seq, s, e, _ = element_in_context(3, tir_len=22, mismatches=(9,))
    t = find_tir(seq, s + ds, e + de)
    assert t is not None
    assert (t["start"], t["end"]) == (s, e)


def test_find_tir_on_element_alone():
    _, s, e, _ = element_in_context(4, tir_len=28)
    seq = element_in_context(4, tir_len=28)[0][s:e]
    t = find_tir(seq)
    assert t is not None and (t["start"], t["end"]) == (0, len(seq))


def test_find_tir_none_on_random_sequence():
    hits = 0
    for seed in range(40):
        rng = random.Random(1000 + seed)
        seq = rand_seq(rng, 2000)
        hits += find_tir(seq, 300, 1700) is not None
    assert hits <= 1          # calibrated false-positive rate ~1-2 %
    rng = random.Random(7)
    assert find_tir(rand_seq(rng, 2000), 300, 1700) is None


def test_find_tir_too_short_or_divergent():
    # an 8 bp perfect TIR is below min_len / min_score
    seq, s, e, _ = element_in_context(5, tir_len=8)
    assert find_tir(seq, s, e) is None
    # a 20 bp TIR with 7 mismatches (65 % identity) is rejected
    seq, s, e, _ = element_in_context(6, tir_len=20, mismatches=(3, 5, 8, 10, 13, 15, 17))
    assert find_tir(seq, s, e) is None


# -------------------------------------------------------------------- find_dr

def test_find_dr_exact_and_longest():
    assert find_dr("TTTTGATCAG", "GATCAGCCCC")[:3] == (6, 0, "GATCAG")
    # the longest matching k wins; shorter suffix/prefix pairs are different registers
    k, mm, a, b = find_dr("ACGTACGTAC", "GTACCCCCCC")
    assert (k, mm, a, b) == (4, 0, "GTAC", "GTAC")


def test_find_dr_none():
    assert find_dr("AAAAAAAAAA", "CCCCCCCCCC") == (0, 0, "", "")
    assert find_dr("", "ACGT") == (0, 0, "", "")


def test_find_dr_mismatch_only_for_long_repeats():
    # 9 bp with one mismatch: accepted (k >= 8)
    assert find_dr("CCCC" + "ATCAGTCAG", "ATCTGTCAG" + "TTTT", max_len=9) == \
        (9, 1, "ATCAGTCAG", "ATCTGTCAG")
    # 5 bp with one mismatch: not accepted
    assert find_dr("CCCCCCCC" + "GATCA", "GATTA" + "TTTTTTT", max_len=5) == (0, 0, "", "")


def test_find_dr_respects_min_len():
    assert find_dr("TTTTTA", "AGGGGG") == (0, 0, "", "")
    assert find_dr("TTTTTA", "AGGGGG", min_len=1)[0] == 1


# ------------------------------------------------------------------- longest ORF

def test_longest_orf_both_strands():
    orf = "ATG" + "GCT" * 120 + "TAA"
    seq = "CCCC" + orf + "GGGG"
    o = longest_orf(seq)
    assert o == {"aa": 121, "start": 4, "end": 4 + len(orf), "strand": "+"}
    rc = revcomp(seq)
    o = longest_orf(rc)
    assert o["strand"] == "-" and o["aa"] == 121
    assert revcomp(rc[o["start"]:o["end"]]) == orf
    assert longest_orf("CCCCCCCC") is None


# --------------------------------------------------------- parsing / clustering

def test_parse_query_id_with_underscores():
    assert parse_query_id("barcode01_r000001_100_1400") == ("barcode01_r000001", 100, 1400)
    assert parse_query_id("a-b_c_5_10") == ("a-b_c", 5, 10)


def test_cluster_intervals_and_partial_absorption():
    items = [("c", 1000, 3000, "a"), ("c", 1002, 2999, "b"), ("c", 998, 3001, "c"),
             ("c", 1000, 1600, "partial"),           # fragmented read, inside the element
             ("c", 8000, 9500, "other"), ("d", 1000, 3000, "other contig")]
    clusters = cluster_intervals(items)
    assert len(clusters) == 3
    main = [c for c in clusters if c["contig"] == "c" and c["seed"][0] < 5000][0]
    assert sorted(m[2] for m in main["members"]) == ["a", "b", "c"]
    assert [m[2] for m in main["partial"]] == ["partial"]


def test_robust_end_and_candidates():
    vals = [100] * 6 + [101] * 3 + [40, 400]
    assert robust_end(vals)[0] == 100
    cands = boundary_candidates(vals, 100)
    assert cands == {100: 6, 101: 3}


def test_characterise_full_element():
    rng = random.Random(9)
    orf = "ATG" + "".join(rng.choice(["GCT", "AAA", "CTG", "GGC", "TCC"]) for _ in range(330)) + "TGA"
    tir = rand_seq(rng, 26)
    body = tir + rand_seq(rng, 100) + "TAA" + orf + rand_seq(rng, 100) + mutate(revcomp(tir), (10,))
    flank_l, flank_r, dr = rand_seq(rng, 400), rand_seq(rng, 400), "CTTAGA"
    seq = flank_l + dr + body + dr + flank_r
    s = len(flank_l) + len(dr)
    e = s + len(body)
    row, element = characterise(seq, s + 2, e - 3)       # estimates a few bp off
    assert (row["start"], row["end"]) == (s, e)
    assert element == body
    assert row["tir_found"] and abs(row["tir_left_len"] - 26) <= 2
    assert (row["dr_len"], row["dr_left_seq"]) == (6, dr)
    assert row["longest_orf_aa"] >= 331     # other frames of the low-complexity ORF may be open too
    assert row["verdict"] == "mobile_element_candidate"


def test_characterise_genome_segment():
    rng = random.Random(10)
    seq = rand_seq(rng, 3000)
    row, _ = characterise(seq, 500, 2500)
    assert not row["tir_found"]
    assert row["verdict"] in ("genome_segment", "partial_evidence")


def test_scan_known_tirs():
    found = {name: t for name, _, t in scan_known_tirs(os.path.join(TEST_DATA, "IS_curated.fasta"))}
    # values observed for test-data/IS_curated.fasta (see module report)
    assert found["IS4_IS186B"]["left_len"] == 23 and found["IS4_IS186B"]["identity"] == 1.0
    assert found["Tn3_ISEc63"]["left_len"] == 35
    assert found["IS5_IS5"]["left_len"] == 16
    assert found["ISNCY_ISPsncy"] is None
    for t in found.values():
        if t:
            assert t["start"] == 0


# ----------------------------------------------------------------- inputs

ANNOT_HEADER = ("Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\t"
                "Subject Start\tSubject End\tNote\tExplained\tref_pos\tis_Reverse\n")


@requires_tools
def test_main_on_handmade_tables(tmp_path):
    """Explicit `samples` list, read IDs with underscores, paths with spaces."""
    rng = random.Random(11)
    genome = rand_seq(rng, 20000)
    tir = rand_seq(rng, 24)
    element = tir + rand_seq(rng, 1400) + revcomp(tir)
    dr = genome[8000:8005]
    genome = genome[:8005] + element + genome[8000:]
    s, e = 8005, 8005 + len(element)
    plasmid = rand_seq(rng, 3000)
    d = tmp_path / "my data"
    (d / "s 1").mkdir(parents=True)
    gpath = tmp_path / "gen ome.fasta"
    gpath.write_text(f">chr_1 test\n{genome}\n")
    reads, rows = [], []
    for i in range(4):
        pos = 1000 + 300 * (i % 2)
        ins = element if i % 2 else revcomp(element)
        r = plasmid[:pos] + ins + plasmid[pos:]
        rid = f"read_{i}_x"
        reads.append(f">{rid}\n{r}\n")
        rows.append(f"{rid}_{pos}_{pos + len(ins)}\tchr_1\t99\t1\t0\t1\t1\t1\t1\tGenome\t99\t{pos}\tFalse\n")
    rows.append("read_9_9_0_10\tIS\t99\t1\t0\t1\t1\t1\t1\tIS_DB\t99\t5\tFalse\n")
    (d / "s 1" / "s 1.fasta").write_text("".join(reads))
    (d / "s 1" / "best_alignment.tab").write_text(ANNOT_HEADER + "".join(rows))
    out = tmp_path / "out dir"
    res = main_novel_elements(genome_fasta=str(gpath), out_dir=str(out),
                              samples=[("s 1", str(d / "s 1" / "best_alignment.tab"),
                                        str(d / "s 1" / "s 1.fasta"))])
    assert len(res) == 1
    r = res[0]
    assert (r["contig"], r["start"], r["end"]) == ("chr_1", s, e)
    assert r["n_reads"] == 4 and r["samples"] == "s 1:4" and r["n_plasmid_sites"] == 2
    assert r["dr_len"] == 5 and r["dr_left_seq"] == dr and r["tir_found"]
    assert r["source"] == "genome" and r["plasmid_dr_len"] == "NA"     # no DR-finder output
    for name in ("novel_elements.tsv", "novel_elements.fasta", "novel_elements.log"):
        assert (out / name).is_file()
    assert sorted(os.listdir(out)) == ["novel_elements.fasta", "novel_elements.log", "novel_elements.tsv"]
    # data_dir discovery gives the same result
    res2 = main_novel_elements(str(d), str(gpath))
    assert [(x["start"], x["end"]) for x in res2] == [(s, e)]
    assert (d / "novel_elements" / "novel_elements.tsv").is_file()


# ------------------------------------------------------------------ end to end

@requires_tools
def test_end_to_end_recovers_simulated_elements(tmp_path, capsys):
    from simulate_traposome import main as simulate
    from src.a_data_mapping import main_mapping
    from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch
    from src.g_blast_annot_batch import main_annot_batch

    data = tmp_path / "e2e data"
    simulate(["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
              "-g", os.path.join(TEST_DATA, "genome.fasta"),
              "-i", os.path.join(TEST_DATA, "IS_curated.fasta"),
              "-o", str(data), "--barcodes", "1", "--reads", "400", "--novel_fraction", "0.5",
              "--novel_mode", "element", "--seed", "10"])
    genome = str(data / "genome_with_elements.fasta")
    shutil.copy(data / "barcode01" / "barcode01_raw.fastq", data / "barcode01" / "barcode01.fastq")
    with contextlib.redirect_stdout(io.StringIO()):
        main_mapping(str(data), os.path.join(TEST_DATA, "plasmid.fasta"), genome)
        main_insert_finder_from_bam_batch(str(data), 500)
        main_annot_batch(str(data), genome, os.path.join(TEST_DATA, "IS_curated.fasta"))
    rows = main_novel_elements(str(data), genome, is_fasta=os.path.join(TEST_DATA, "IS_curated.fasta"))

    with open(data / "elements.tsv") as fh:
        elements = list(csv.DictReader(fh, delimiter="\t"))
    matched = set()
    report = []
    for el in elements:
        s, e = int(el["genome_insert_pos"]), int(el["genome_end"])
        hits = [r for r in rows if r["contig"] == el["genome_contig"]
                and min(r["end"], e) - max(r["start"], s) > 0.5 * (e - s)]
        assert len(hits) == 1, f"{el['element_id']}: {len(hits)} clusters"
        r = hits[0]
        matched.add(r["cluster_id"])
        report.append((el["element_id"], r["start"] - s, r["end"] - e,
                       r["tir_left_len"] - int(el["tir_len"]), r["dr_len"], int(el["genome_dr_len"])))
        assert abs(r["start"] - s) <= 3 and abs(r["end"] - e) <= 3
        assert abs(r["tir_left_len"] - int(el["tir_len"])) <= 4
        assert r["dr_len"] == int(el["genome_dr_len"]) and r["dr_left_seq"] == el["genome_dr_seq"]
        assert r["verdict"] == "mobile_element_candidate"
        orf_aa = (int(el["orf_end"]) - int(el["orf_start"])) // 3 - 1
        assert r["longest_orf_aa"] == orf_aa
    # Other clusters may only be host copies of curated IS elements whose
    # inserts were labelled Genome; they must be flagged by known_hit.
    others = [r for r in rows if r["cluster_id"] not in matched]
    assert all(r["source"] == "genome" for r in rows)       # the genome contains the elements
    assert all(r["known_hit"] for r in others)
    with capsys.disabled():
        print("\n[novel elements] (element, start err, end err, TIR len err, DR found, DR true):", report,
              "| other clusters:", [(r["cluster_id"], r["known_hit"]) for r in others])
