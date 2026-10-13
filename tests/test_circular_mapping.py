"""Unit tests for src/circular_mapping.py (AUDIT B09).

Coordinates: the rotated reference starts at original position s = L // 2,
so rotated position b is original (b + s) % L and the original origin lies at
rotated position cut = L - s.
"""
import random

import pysam
import pytest

from conftest import REPO  # noqa: F401  (puts the repo on sys.path)
from src import circular_mapping as cm
from src.e_insert_finder_from_bam import segment_of

L = 1000                # s = 500, cut = 500
NAME = "plasmid"
FLAG_REVERSE, FLAG_SECONDARY, FLAG_SUPPLEMENTARY = 16, 256, 2048
rng = random.Random(3)
PLASMID = "".join(rng.choice("ACGT") for _ in range(L))
ROTATED = PLASMID[L // 2:] + PLASMID[:L // 2]


def header(length=L):
    return pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"},
                                            "SQ": [{"SN": NAME, "LN": length}]})


def record(hdr, cigar, pos, flag=0, seq=None, name="r"):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = cigar
    if seq is not None:
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.set_tag("AS", 100)
    a.set_tag("tp", "P", value_type="A")
    return a


def test_rotate_fasta(tmp_path):
    src = tmp_path / "p.fasta"
    src.write_text(f">{NAME} some description\n{PLASMID[:400]}\n{PLASMID[400:]}\n>small\nACGTAC\n")
    out = tmp_path / "r.fasta"
    cm.write_rotated_fasta(str(src), str(out))
    recs = cm.read_fasta(str(out))
    assert recs == [(NAME, ROTATED), ("small", "TACACG")]


def test_record_not_crossing_origin_is_shifted():
    hdr = header()
    # rotated 100..300 -> original 600..800; rotated 600..700 -> original 100..200
    a = cm.convert_record(record(hdr, "50S200M", 100, seq=ROTATED[50:300]), hdr, L, PLASMID)
    b = cm.convert_record(record(hdr, "100M", 600, seq=ROTATED[600:700]), hdr, L, PLASMID)
    assert [(r.reference_start, r.cigarstring, r.get_tag("AS")) for r in a] == [(600, "50S200M", 100)]
    assert [(r.reference_start, r.cigarstring) for r in b] == [(100, "100M")]


def test_forward_record_split_at_origin():
    # Read = 30 bp clip + rotated 400..650 (= original 900..1000 + 0..150) + 20 bp clip.
    hdr = header()
    seq = "A" * 30 + ROTATED[400:650] + "C" * 20
    out = cm.convert_record(record(hdr, "30S250M20S", 400, seq=seq), hdr, L, PLASMID)
    assert [(r.reference_start, r.reference_end, r.cigarstring, r.is_supplementary) for r in out] == [
        (900, 1000, "30H100M170H", True),     # smaller piece -> supplementary, hard clipped
        (0, 150, "130S150M20S", False),       # larger piece keeps primary flag and full SEQ
    ]
    assert out[1].query_sequence == seq
    assert out[0].query_sequence == seq[30:130]
    assert all(r.get_tag("NM") == 0 for r in out)
    assert not any(r.has_tag("AS") for r in out)
    # pieces are adjacent on the read
    segs = sorted((segment_of(r).q_start, segment_of(r).q_end) for r in out)
    assert segs == [(30, 130), (130, 280)]
    # every aligned base matches the original reference
    for r in out:
        for q, p in r.get_aligned_pairs(matches_only=True):
            assert r.query_sequence[q] == PLASMID[p]


def test_reverse_record_split_keeps_read_orientation():
    hdr = header()
    seq = ROTATED[450:560]                       # alignment orientation
    out = cm.convert_record(record(hdr, "110M", 450, flag=FLAG_REVERSE, seq=seq), hdr, L, PLASMID)
    assert [(r.reference_start, r.cigarstring, r.is_reverse) for r in out] == [
        (950, "50M60H", True), (0, "50S60M", True)]
    # original read orientation: rotated 500..560 (original 0..60) comes first
    assert sorted((segment_of(r).q_start, segment_of(r).q_end) for r in out) == [(0, 60), (60, 110)]


def test_insertion_and_deletion_at_the_cut():
    hdr = header()
    # Insertion exactly at the cut -> becomes clipping on both pieces (an
    # uncovered read interval remains, as for an inserted element).
    out = cm.convert_record(record(hdr, "100M600I100M", 400), hdr, L)
    assert [(r.reference_start, r.cigarstring) for r in out] == [(900, "100M700H"), (0, "700H100M")]
    # Deletion straddling the cut: the deleted bases are dropped.
    out = cm.convert_record(record(hdr, "90M20D100M", 400), hdr, L)
    assert [(r.reference_start, r.cigarstring) for r in out] == [(900, "90M100H"), (10, "90H100M")]


def test_supplementary_hard_clipped_split():
    hdr = header()
    seq = ROTATED[480:540]
    rec = record(hdr, "1000H60M7H", 480, flag=FLAG_SUPPLEMENTARY, seq=seq)
    out = cm.convert_record(rec, hdr, L, PLASMID)
    assert [(r.reference_start, r.cigarstring, r.is_supplementary, r.query_sequence) for r in out] == [
        (980, "1000H20M47H", True, seq[:20]), (0, "1020H40M7H", True, seq[20:])]


def test_choose_prefers_linear_unless_rotated_covers_more():
    hdr = header()
    lin = [record(hdr, "400M300S", 0)]                                    # ends at 400, not an edge
    assert cm.choose(lin, [record(hdr, "405M295S", 0)]) == "linear"       # within MIN_GAIN
    assert cm.choose(lin, [record(hdr, "420M280S", 0)]) == "linear"       # small gain, no contig edge
    assert cm.choose(lin, [record(hdr, "450M250S", 0)]) == "rotated"      # > LARGE_GAIN
    # a long CIGAR insertion is not aligned sequence
    assert cm.choose(lin, [record(hdr, "200M200I200M100S", 0)]) == "linear"
    assert cm.choose([], lin) == "rotated"
    assert cm.choose(lin, []) == "linear"


def test_choose_small_gain_used_when_linear_alignment_stops_at_origin():
    hdr = header()
    # Linear: 300 bp ending at the contig end (1000) in the middle of the read;
    # the 20 bp after the origin are not aligned. Rotated aligns them.
    lin = [record(hdr, "300M400S", 700)]
    assert cm.edge_inside_read(lin)
    assert cm.choose(lin, [record(hdr, "320M380S", 180)]) == "rotated"
    # The same alignment ending at the read end is not an origin crossing.
    assert not cm.edge_inside_read([record(hdr, "10S300M", 700)])
    # reverse strand: the contig end is at the read start of the alignment
    assert cm.edge_inside_read([record(hdr, "300M400S", 700, flag=FLAG_REVERSE)])
    assert not cm.edge_inside_read([record(hdr, "400S300M", 700, flag=FLAG_REVERSE)])
    # A plain origin crossing (the next alignment starts where this one stops
    # at the contig end) is not a cut alignment.
    assert not cm.edge_inside_read([record(hdr, "300M400S", 700),
                                    record(hdr, "300H400M", 0, flag=FLAG_SUPPLEMENTARY)])
    assert cm.choose([record(hdr, "300M400S", 700), record(hdr, "300H400M", 0, flag=FLAG_SUPPLEMENTARY)],
                     [record(hdr, "720M", 180)]) == "linear"


def write_bam(path, hdr, records):
    with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
        for r in records:
            out.write(r)


def unmapped(hdr, name, seq):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.flag = 4
    a.reference_id = -1
    a.reference_start = -1
    a.query_sequence = seq
    return a


def test_merge_picks_per_read_and_converts(tmp_path):
    hdr = header()
    fasta = tmp_path / "p.fasta"
    fasta.write_text(f">{NAME}\n{PLASMID}\n")
    # r1: crosses the origin; linear mapping has only the long piece (the
    # 100 bp after the origin are unaligned), rotated mapping aligns all.
    r1 = PLASMID[600:] + PLASMID[:100]
    # r2: away from the origin; both mappings agree -> linear record kept.
    r2 = PLASMID[200:600]
    lin = [record(hdr, "400M100S", 600, seq=r1, name="r1"),
           record(hdr, "400M", 200, seq=r2, name="r2"),
           unmapped(hdr, "r3", "ACGT" * 50)]
    rot = [record(hdr, "500M", 100, seq=r1, name="r1"),
           record(hdr, "400M", 700, seq=r2, name="r2"),
           unmapped(hdr, "r3", "ACGT" * 50)]
    write_bam(tmp_path / "lin.bam", hdr, lin)
    write_bam(tmp_path / "rot.bam", hdr, rot)
    counts = cm.merge(str(fasta), str(tmp_path / "lin.bam"), str(tmp_path / "rot.bam"),
                      str(tmp_path / "out.bam"))
    assert counts == {"reads": 3, "mapped": 2, "rotated_used": 1, "split_records": 1}
    with pysam.AlignmentFile(str(tmp_path / "out.bam")) as bam:
        assert bam.header.get_reference_name(0) == NAME and bam.lengths == (L,)
        got = [(r.query_name, r.reference_start, r.cigarstring, r.is_supplementary) for r in bam]
    assert got == [("r1", 600, "400M100S", False), ("r1", 0, "400H100M", True),
                   ("r2", 200, "400M", False)]


def test_merge_rejects_out_of_step_inputs(tmp_path):
    hdr = header()
    fasta = tmp_path / "p.fasta"
    fasta.write_text(f">{NAME}\n{PLASMID}\n")
    write_bam(tmp_path / "lin.bam", hdr, [unmapped(hdr, "a", "ACGT")])
    write_bam(tmp_path / "rot.bam", hdr, [unmapped(hdr, "b", "ACGT")])
    with pytest.raises(ValueError, match="read order"):
        cm.merge(str(fasta), str(tmp_path / "lin.bam"), str(tmp_path / "rot.bam"), str(tmp_path / "o.bam"))


def test_merge_rejects_wrong_plasmid_length(tmp_path):
    hdr = header(L + 1)
    fasta = tmp_path / "p.fasta"
    fasta.write_text(f">{NAME}\n{PLASMID}\n")
    write_bam(tmp_path / "lin.bam", hdr, [])
    write_bam(tmp_path / "rot.bam", hdr, [])
    with pytest.raises(ValueError, match="length"):
        cm.merge(str(fasta), str(tmp_path / "lin.bam"), str(tmp_path / "rot.bam"), str(tmp_path / "o.bam"))
