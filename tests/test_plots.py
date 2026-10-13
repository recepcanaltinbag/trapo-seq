"""Tests for the plotting / BAM-statistics modules (a, b, c, c_01, d, f)."""
import os

import numpy as np
import pysam
import pytest

from src import (a_read_histograms, b_filtering_based_on_len, c_01_polar_plot, c_map_dist,
                 d_kde_mobile, f_in_del_plot)

FLAG_REVERSE, FLAG_UNMAPPED, FLAG_SECONDARY, FLAG_SUPPLEMENTARY = 16, 4, 256, 2048


def make_bam(path, contigs, records):
    """contigs: [(name, length)]; records: [(qname, flag, contig, pos, cigar)].
    Writes a coordinate-sorted, indexed BAM."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in contigs]}
    tid = {n: i for i, (n, _) in enumerate(contigs)}
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for qname, flag, contig, pos, cigar in records:
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.flag = flag
            if flag & FLAG_UNMAPPED:
                a.reference_id = -1
                a.reference_start = -1
                a.query_sequence = "A" * 50
            else:
                a.reference_id = tid[contig]
                a.reference_start = pos
                a.cigarstring = cigar
                a.mapping_quality = 60
                a.query_sequence = "A" * a.infer_query_length()
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
            out.write(a)
    pysam.sort("-o", str(path), unsorted)
    os.remove(unsorted)
    pysam.index(str(path))
    return str(path)


def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def segment(cigar, flag=0):
    a = pysam.AlignedSegment()
    a.flag = flag
    a.cigarstring = cigar
    return a


# ------------------------------------------------------------------ B48 histogram
def write_fastq(path, lengths):
    with open(path, "w") as f:
        for i, n in enumerate(lengths):
            f.write(f"@r{i}\n{'A' * n}\n+\n{'I' * n}\n")
    return str(path)


def test_histogram_zero_width_bins(tmp_path):
    fq = write_fastq(tmp_path / "a.fastq", [100] * 10)
    out = a_read_histograms.plot_histogram(fq, str(tmp_path / "hist"))
    assert out.endswith("hist.pdf") and nonempty(out)


def test_histogram_empty_fastq(tmp_path):
    fq = write_fastq(tmp_path / "a.fastq", [])
    out = a_read_histograms.plot_histogram(fq, str(tmp_path / "hist.png"))
    assert out == str(tmp_path / "hist.png") and nonempty(out)


def test_histogram_keeps_extension(tmp_path):
    fq = write_fastq(tmp_path / "a.fastq", [50, 100, 2000, 3000])
    out = a_read_histograms.plot_histogram(fq, str(tmp_path / "x.pdf"))
    assert out == str(tmp_path / "x.pdf") and nonempty(out)
    assert not os.path.exists(str(tmp_path / "x.pdf.pdf"))


def test_filter_by_length(tmp_path):
    fq = write_fastq(tmp_path / "a.fastq", [50, 100, 150])
    out = tmp_path / "f.fastq"
    assert b_filtering_based_on_len.filter_fastq_by_length(fq, str(out), 100) == 2
    assert a_read_histograms.calculate_read_lengths(str(out)) == [100, 150]


# ------------------------------------------------------------------ B44 read coordinates
@pytest.mark.parametrize("cigar, flag, expected", [
    ("200M300S", 0, (0, 200, 500)),
    ("200M300S", FLAG_REVERSE, (300, 500, 500)),
    ("300H200M", FLAG_SUPPLEMENTARY, (300, 500, 500)),
    # hard-clipped supplementary on the reverse strand: read bases 0..200 in read orientation
    ("300H200M", FLAG_SUPPLEMENTARY | FLAG_REVERSE, (0, 200, 500)),
    ("100H50M10I40=300H", FLAG_SUPPLEMENTARY, (100, 200, 500)),
    ("100H50M10I40X300H", FLAG_SUPPLEMENTARY | FLAG_REVERSE, (300, 400, 500)),
])
def test_read_query_interval(cigar, flag, expected):
    assert d_kde_mobile.read_query_interval(segment(cigar, flag)) == expected


def kde_dataset(tmp_path):
    d = tmp_path / "bc 1"
    d.mkdir()
    # r1 (1000 bp): read[300,600) plasmid; read[0,300) genome fwd; read[600,900) genome rev, hard-clipped
    plasmid = make_bam(d / "03_sorted_mapped_to_plasmid.bam", [("pl", 5000)], [
        ("r1", FLAG_SUPPLEMENTARY, "pl", 100, "300H300M400H"),
        ("r2", 0, "pl", 1000, "3000M"),
    ])
    make_bam(d / "07_sorted_aligned_to_other_genome.bam", [("chrA", 20000), ("chrB", 5000)], [
        ("r1", 0, "chrA", 1000, "300M700S"),
        ("r1", FLAG_SUPPLEMENTARY | FLAG_REVERSE, "chrB", 2000, "100H300M600H"),
        ("r1", FLAG_SECONDARY, "chrA", 5000, "300M700H"),
    ])
    return plasmid, d


def test_kde_hard_clipped_reverse_supplementary_kept(tmp_path):
    plasmid, d = kde_dataset(tmp_path)
    reads = d_kde_mobile.read_unique_alignment_genome(str(d / "07_sorted_aligned_to_other_genome.bam"), "g")
    assert sorted((s, e) for _, s, e in reads["r1"]) == [(0, 300), (600, 900)]

    d_kde_mobile.main_kde_mobile(plasmid, "kde", 100)
    rows = dict(l.split("\t", 1) for l in open(d / "kde.tab").read().splitlines()[1:])
    assert rows["r1"] == "300\t600"
    assert rows["r2"] == "3000\t0"
    assert nonempty(d / "kde.pdf") and nonempty(d / "kde_scatter.pdf")


def test_kde_no_read_passes_threshold(tmp_path):  # B46
    plasmid, d = kde_dataset(tmp_path)
    d_kde_mobile.main_kde_mobile(plasmid, "kde", 10 ** 9)
    assert nonempty(d / "kde.pdf") and nonempty(d / "kde_scatter.pdf")


def test_kde_plot_uses_filtered_reads(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    dist = {f"r{i}": {"plasmid": int(p), "genome": int(g)}
            for i, (p, g) in enumerate(zip(rng.integers(2500, 6000, 40), rng.integers(0, 3000, 40)))}
    dist["low"] = {"plasmid": 10, "genome": 99999}
    seen = {}
    from matplotlib.axes import Axes
    orig = Axes.scatter
    monkeypatch.setattr(Axes, "scatter", lambda self, x, y, **kw: seen.setdefault("x", list(x)) and orig(self, x, y, **kw))
    assert d_kde_mobile.plot_read_distribution(dist, str(tmp_path / "k"), 2000) == 40
    assert 10 not in seen["x"]
    assert nonempty(tmp_path / "k.pdf") and nonempty(tmp_path / "k_scatter.pdf")


# ------------------------------------------------------------------ B45 multi-contig coverage
def test_compute_coverage_per_contig():
    reads = {"a": [("chrA", 0, 10)], "b": [("chrB", 5, 8), ("chrB", 6, 20)], "c": [("chrX", 0, 5)]}
    cov = c_map_dist.compute_coverage(reads, {"chrA": 10, "chrB": 15})
    assert cov["chrA"].tolist() == [1] * 10
    assert cov["chrB"].tolist() == [0] * 5 + [1, 2, 2] + [1] * 7
    assert set(cov) == {"chrA", "chrB"}


def test_map_dist_multi_contig(tmp_path):
    d = tmp_path / "bc 1"
    d.mkdir()
    plasmid = make_bam(d / "03_sorted_mapped_to_plasmid.bam", [("pl", 3000)], [
        ("r1", 0, "pl", 0, "500M"), ("r2", 0, "pl", 100, "500M"), ("r3", 0, "pl", 200, "500M"),
    ])
    genome = make_bam(d / "07_sorted_aligned_to_other_genome.bam", [("chrA", 1000), ("chrB", 500)], [
        ("r1", 0, "chrB", 100, "100M"),
        ("r1", FLAG_SECONDARY, "chrA", 0, "100M"),
        ("r4", 0, "chrA", 900, "100M"),
    ])
    reads = c_map_dist.read_alignment(genome)
    cov = c_map_dist.compute_coverage(reads, c_map_dist.get_contig_lengths(genome))
    assert cov["chrB"][100:200].tolist() == [1] * 100 and cov["chrB"].sum() == 100
    assert cov["chrA"][:900].sum() == 0 and cov["chrA"][900:].tolist() == [1] * 100

    fq = write_fastq(tmp_path / "r.fastq", [500] * 5)
    c_map_dist.main_map_dist(plasmid, "dist", fastq_file=fq)
    assert nonempty(d / "dist.pdf")
    stats = (d / "dist_mapping_stats.txt").read_text()
    assert "Total reads: 5" in stats and "Both: 1" in stats and "Plasmid only: 2" in stats
    assert "Transposition Rate (both / plasmid-mapped): 0.33" in stats


def test_summarize_mappings_no_plasmid_reads(tmp_path):
    res = c_map_dist.summarize_mappings({}, {"g": [("c", 0, 1)]}, str(tmp_path / "s.txt"))
    assert res["rate"] == "NA"


def test_map_dist_downsample_keeps_peak():
    y = np.zeros(100000, dtype=int)
    y[54321] = 7
    xs, ys = c_map_dist._downsample(np.arange(len(y)), y)
    assert len(ys) <= c_map_dist.MAX_PLOT_POINTS + 1 and ys.max() == 7


# ------------------------------------------------------------------ B05-style CIGAR in f_in_del_plot
def test_in_del_cigar_eq_x_n(tmp_path):
    bam = make_bam(tmp_path / "p.bam", [("pl", 1000)], [
        # ref: 100 +10(=) +2(X) +5(N) = 117 -> D(3) at 117, ref -> 120; I(4) at 120
        ("r1", 0, "pl", 100, "10=2X5N3D4I10M"),
        ("r2", FLAG_SECONDARY, "pl", 100, "10=2X5N3D4I10M"),
        ("u", FLAG_UNMAPPED, None, None, None),
    ])
    dels, wdels = f_in_del_plot.extract_deletions_from_bam(bam, 3)
    assert dict(dels) == {117: 1}
    assert dict(wdels) == {117: 3, 118: 3, 119: 3}
    ins, wins = f_in_del_plot.extract_insertions_from_bam(bam, 2)
    assert dict(ins) == {120: 1} and dict(wins) == {120: 4}


def test_in_del_plot_smoke(tmp_path):
    bam = make_bam(tmp_path / "p.bam", [("pl", 1000)], [("r1", 0, "pl", 100, "50M10D50M5I50M")])
    tab = tmp_path / "insertion_from_bam.tab"
    tab.write_text("r9\t10\t1010\t1000\t300\t20\tSC\tFalse\n")
    f_in_del_plot.main_in_del_plot(str(tab), bam, "in_del_plot")
    for name in ("non-weighted", "weighted"):
        assert nonempty(tmp_path / f"in_del_plot_{name}.pdf")
        assert nonempty(tmp_path / f"in_del_plot_{name}_insertion_table.tsv")
    ins = (tmp_path / "in_del_plot_non-weighted_insertion_table.tsv").read_text().splitlines()
    assert ins[1:] == ["210\t1", "300\t1"]  # 100 + 50M + 10D + 50M


def test_in_del_plot_empty_inputs(tmp_path):
    bam = make_bam(tmp_path / "p.bam", [("pl", 1000)], [])
    tab = tmp_path / "insertion_from_bam.tab"
    tab.write_text("")
    f_in_del_plot.main_in_del_plot(str(tab), bam, "o")
    assert nonempty(tmp_path / "o_weighted.pdf")


# ------------------------------------------------------------------ polar plot
@pytest.mark.parametrize("arr", [np.zeros(50), np.arange(30), np.random.default_rng(1).integers(0, 9, 3_000_000)])
def test_polar_plot(tmp_path, arr):
    out = c_01_polar_plot.main_polar(arr, str(tmp_path / "polar.pdf"))
    assert nonempty(out)


def test_polar_plot_empty(tmp_path):
    with pytest.raises(ValueError):
        c_01_polar_plot.main_polar([], str(tmp_path / "polar.pdf"))


def _fastq(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name} extra=desc\n{seq}\n+\n{'I' * len(seq)}\n")


def test_filter_matches_biopython_and_keeps_bytes(tmp_path):
    from Bio import SeqIO
    fq = tmp_path / "in.fastq"
    _fastq(fq, [("a", "ACGT" * 30), ("b", "AC"), ("c", "G" * 100), ("d", "T" * 99)])
    out = tmp_path / "out.fastq"
    assert b_filtering_based_on_len.filter_fastq_by_length(str(fq), str(out), 100) == 2
    ids = [r.id for r in SeqIO.parse(str(out), "fastq")]
    assert ids == ["a", "c"]
    assert "@a extra=desc\n" in out.read_text()  # header copied unchanged


def test_filter_gzip_input(tmp_path):
    import gzip
    fq = tmp_path / "in.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("@a\n" + "A" * 50 + "\n+\n" + "I" * 50 + "\n")
    out = tmp_path / "out.fastq"
    assert b_filtering_based_on_len.filter_fastq_by_length(str(fq), str(out), 10) == 1


def test_filter_multiline_fastq_falls_back(tmp_path):
    fq = tmp_path / "in.fastq"
    fq.write_text("@a\nACGT\nACGT\n+\nIIII\nIIII\n@b\nAC\n+\nII\n")
    out = tmp_path / "out.fastq"
    assert b_filtering_based_on_len.filter_fastq_by_length(str(fq), str(out), 5) == 1
    assert "ACGTACGT" in out.read_text()
    assert not list(tmp_path.glob("*.tmp*"))
