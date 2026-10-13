"""Unit tests for src/e_insert_finder_from_bam.py and its batch wrapper.

Tiny BAMs with hand-written CIGARs are built with pysam; every expected
coordinate is derived by hand in the comments.

Output tuple: (query_start, query_end, length, ref_pos, insert_type, is_reverse),
query coordinates 0-based half-open in the original read orientation.
"""
import os
import random

import pysam
import pytest

from conftest import REPO  # noqa: F401  (puts the repo on sys.path)
from src.e_insert_finder_from_bam import (
    OUTPUT_NAME,
    find_insertions,
    main_insert_finder_from_bam,
)
from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch

REF_LEN = 7000
FLAG_REVERSE, FLAG_SECONDARY, FLAG_SUPPLEMENTARY = 16, 256, 2048


def query_len(cigar):
    """Length of SEQ for a CIGAR string (M/I/S/=/X; hard clips excluded)."""
    rec = pysam.AlignedSegment()
    rec.cigarstring = cigar
    return sum(n for op, n in rec.cigartuples if op in (0, 1, 4, 7, 8))


def make_bam(path, records):
    """records: dicts with name, cigar, pos, flag (default 0), qual (default 30),
    seq (default random; None for '*')."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "plasmid", "LN": REF_LEN}]}
    rng = random.Random(1)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in sorted(records, key=lambda r: r["pos"]):
            a = pysam.AlignedSegment()
            a.query_name = r["name"]
            a.flag = r.get("flag", 0)
            a.reference_id = 0
            a.reference_start = r["pos"]
            a.mapping_quality = 60
            a.cigarstring = r["cigar"]
            if r.get("seq", "random") is not None:
                n = query_len(r["cigar"])
                a.query_sequence = "".join(rng.choice("ACGT") for _ in range(n))
                a.query_qualities = pysam.qualitystring_to_array(chr(33 + r.get("qual", 30)) * n)
            out.write(a)
    return str(path)


def rows(bam, threshold=500):
    return [(i.query_start, i.query_end, i.length, i.ref_pos, i.insert_type, i.is_reverse)
            for i in find_insertions(bam, threshold)]


# --------------------------------------------------------------- CIGAR (IN)

def test_cigar_insertion_forward(tmp_path):
    # q: 500 M, then 700 I at q=500..1200; ref after the 500 M = 100 + 500 = 600
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "500M700I800M", "pos": 100}])
    assert rows(bam) == [(500, 1200, 700, 600, "IN", False)]


def test_cigar_insertion_reverse(tmp_path):
    # read length 2000; alignment-orientation q=500..1200 -> original 2000-1200 .. 2000-500
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "500M700I800M", "pos": 100,
                                         "flag": FLAG_REVERSE}])
    assert rows(bam) == [(800, 1500, 700, 600, "IN", True)]


def test_insertion_at_threshold_not_reported(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "500M500I800M", "pos": 100}])
    assert rows(bam, threshold=500) == []
    assert rows(bam, threshold=499) == [(500, 1000, 500, 600, "IN", False)]


def test_eqx_ops_are_matches(tmp_path):
    # B05: 300= 1X 199= equals 500M; deletions/N still advance the reference.
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "300=1X199=700I10D800=", "pos": 100}])
    assert rows(bam) == [(500, 1200, 700, 600, "IN", False)]


def test_hard_clipped_reverse_supplementary_insertion(tmp_path):
    # supplementary, reverse: 1000H 500M 600I 900M (read length 3000).
    # aligned part: aln-orientation q 1000..3000 -> original 0..2000.
    # insertion: aln-orientation q 1500..2100 -> original 3000-2100 .. 3000-1500 = 900..1500,
    # ref = 100 + 500 = 600. The primary covers original q 2000..3000.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "2000S1000M", "pos": 5000},
        {"name": "r", "cigar": "1000H500M600I900M", "pos": 100,
         "flag": FLAG_REVERSE | FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == [(900, 1500, 600, 600, "IN", True)]


def test_insertions_from_several_records_are_all_kept(tmp_path):
    # B02: the IN of the primary must not be overwritten by the supplementary's.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "500M700I800M2000S", "pos": 100},                      # q 0..2000
        {"name": "r", "cigar": "2000H500M600I900M", "pos": 3000, "flag": FLAG_SUPPLEMENTARY},  # q 2000..4000
    ])
    assert sorted(rows(bam)) == [(500, 1200, 700, 600, "IN", False),
                                 (2500, 3100, 600, 3500, "IN", False)]


# -------------------------------------------------------- uncovered (SC)

def test_soft_clip_forward(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "1000M800S", "pos": 2000}])
    # gap q 1000..1800, flanked on the left by ref_end = 3000
    assert rows(bam) == [(1000, 1800, 800, 3000, "SC", False)]


def test_soft_clip_reverse(tmp_path):
    # reverse: the trailing clip is at the read start in original orientation;
    # the read-start end of the alignment corresponds to ref_end = 3000.
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "1000M800S", "pos": 2000,
                                         "flag": FLAG_REVERSE}])
    assert rows(bam) == [(0, 800, 800, 3000, "SC", True)]


def test_short_soft_clip_ignored(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "1000M500S", "pos": 2000}])
    assert rows(bam) == []


def test_split_read_gap_between_primary_and_supplementary(tmp_path):
    # primary q 0..1000 at ref 100..1100, supplementary q 1700..2700 at ref 1090..2090
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "1000M1700S", "pos": 100},
        {"name": "r", "cigar": "1700H1000M", "pos": 1090, "flag": FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == [(1000, 1700, 700, 1100, "SC", False)]


def test_strand_comes_from_flanking_alignment(tmp_path):
    # B04: read length 2800.
    #   primary   fwd  1000M1800S at ref 5000 -> q 0..1000
    #   supp      rev  800H1000M1000H at ref 200 -> aln q 800..1800 -> original 1000..2000
    # gap q 2000..2800 is flanked only by the reverse supplementary; its
    # read-orientation right end corresponds to ref_start = 200.
    # (The legacy code used the strand of the last record in BAM order: False.)
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "1000M1800S", "pos": 5000},
        {"name": "r", "cigar": "800H1000M1000H", "pos": 200,
         "flag": FLAG_REVERSE | FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == [(2000, 2800, 800, 200, "SC", True)]


def test_flank_at_contig_edge_is_avoided(tmp_path):
    # Circular plasmid, insertion shortly after the origin:
    #   primary fwd 3000M2000S at ref 4000 -> ref_end 7000 = contig end (q 0..3000)
    #   supp    fwd 4000H1000M at ref 490 (q 4000..5000)
    # The real site is the right flank start (490), not the contig end.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "3000M2000S", "pos": 4000},
        {"name": "r", "cigar": "4000H1000M", "pos": 490, "flag": FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == [(3000, 4000, 1000, 490, "SC", False)]


def test_insertion_exactly_at_origin_reported_as_zero(tmp_path):
    # B09: with the circular mapping an alignment crossing the origin is split
    # there. Both flanks of a gap at the origin sit on a contig edge; the
    # contig end (7000 == REF_LEN) is reported as position 0.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "3000M2000S", "pos": 4000},
        {"name": "r", "cigar": "4000H1000M", "pos": 0, "flag": FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == [(3000, 4000, 1000, 0, "SC", False)]


def test_adjacent_pieces_split_at_origin_give_no_gap(tmp_path):
    # A read crossing the origin, split into two records with no read bases
    # between them, contains no insertion.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "3000M1000S", "pos": 4000},
        {"name": "r", "cigar": "3000H1000M", "pos": 0, "flag": FLAG_SUPPLEMENTARY},
    ])
    assert rows(bam) == []


# -------------------------------------------------------------- secondary

def test_secondary_insertion_not_duplicated_and_quality_from_primary(tmp_path):
    # B03: secondary records (often SEQ '*') must not duplicate rows or zero the quality.
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "500M700I800M", "pos": 100, "qual": 25},
        {"name": "r", "cigar": "500M700I800M", "pos": 3000, "flag": FLAG_SECONDARY, "seq": None},
    ])
    ins = find_insertions(bam, 500)
    assert [(i.query_start, i.query_end, i.ref_pos, i.quality) for i in ins] == [(500, 1200, 600, 25)]


def test_secondary_does_not_hide_gap(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "1000M1000S", "pos": 100},
        {"name": "r", "cigar": "1000S1000M", "pos": 3000, "flag": FLAG_SECONDARY},
    ])
    assert rows(bam) == [(1000, 2000, 1000, 1100, "SC", False)]


def test_quality_from_primary_not_supplementary(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r", "cigar": "1000M1700S", "pos": 100, "qual": 20},
        {"name": "r", "cigar": "1700H1000M", "pos": 1090, "flag": FLAG_SUPPLEMENTARY, "qual": 40},
    ])
    assert [i.quality for i in find_insertions(bam, 500)] == [20]


def test_unmapped_and_fully_aligned_reads_give_nothing(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "3000M", "pos": 100}])
    assert rows(bam) == []


# ----------------------------------------------------------------- output

def test_output_file_format(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [
        {"name": "r1", "cigar": "500M700I800M", "pos": 100, "qual": 30},
        {"name": "r2", "cigar": "1000M800S", "pos": 2000, "flag": FLAG_REVERSE, "qual": 12},
    ])
    out = tmp_path / "dir with space" / "out.tab"
    out.parent.mkdir()
    main_insert_finder_from_bam(bam, str(out), 500)
    assert out.read_text().splitlines() == [
        "r1\t500\t1200\t700\t600\t30\tIN\tFalse",
        "r2\t0\t800\t800\t3000\t12\tSC\tTrue",
    ]
    assert os.listdir(out.parent) == ["out.tab"]       # no temp files left


def test_batch_matches_single(tmp_path):
    data = tmp_path / "data dir"
    for name in ("s1", "s2", "empty"):
        (data / name).mkdir(parents=True)
    recs = [{"name": "r1", "cigar": "500M700I800M", "pos": 100},
            {"name": "r2", "cigar": "1000M800S", "pos": 2000}]
    make_bam(data / "s1" / "03_sorted_mapped_to_plasmid.bam", recs)
    make_bam(data / "s2" / "03_sorted_mapped_to_plasmid.bam", recs[:1])
    result = main_insert_finder_from_bam_batch(str(data), 500)
    assert sorted(result) == ["s1", "s2"]
    assert not (data / "empty" / OUTPUT_NAME).exists()
    single = tmp_path / "single.tab"
    main_insert_finder_from_bam(str(data / "s1" / "03_sorted_mapped_to_plasmid.bam"), str(single), 500)
    assert (data / "s1" / OUTPUT_NAME).read_text() == single.read_text()
    assert len((data / "s2" / OUTPUT_NAME).read_text().splitlines()) == 1


def test_batch_module_reuses_single_module_logic():
    import src.e_insert_finder_from_bam as single
    import src.e_insert_finder_from_bam_batch as batch
    assert batch.main_insert_finder_from_bam is single.main_insert_finder_from_bam
    for dup in ("get_read_info", "get_middle_inserts", "find_large_insertions", "find_zero_sequences"):
        assert not hasattr(batch, dup)


def test_no_stdout_redirection(tmp_path, capsys):
    import sys
    before = sys.stdout
    bam = make_bam(tmp_path / "a.bam", [{"name": "r", "cigar": "1000M800S", "pos": 2000}])
    main_insert_finder_from_bam(bam, str(tmp_path / "o.tab"), 500)
    assert sys.stdout is before
    assert capsys.readouterr().out == ""
