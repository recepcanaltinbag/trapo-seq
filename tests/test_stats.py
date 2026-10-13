"""Tests for IS statistics (h_is_stats), the heatmap (i_heatmap) and in/del count merging."""
import os
from collections import Counter

import pandas as pd
import pytest

from src import h_is_stats, i_heatmap, f_in_del_plot

COLUMNS = ["Query ID", "Subject ID", "Identity (%)", "Score", "E-value", "Query Start",
           "Query End", "Subject Start", "Subject End", "Note", "Explained", "ref_pos", "is_Reverse"]


def write_best(path, rows):
    """rows: (query_id, subject_id, score, note)"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    recs = [[q, s, 99.0, sc, 0.0, 1, 100, 1, 100, note, 1.0, 10, False] for q, s, sc, note in rows]
    pd.DataFrame(recs, columns=COLUMNS).to_csv(path, sep="\t", index=False)


# ---------------------------------------------------------------- B40
@pytest.mark.parametrize("subject, expected", [
    ("ISPs1", ("ISPs1", "ISPs1")),            # 0 underscores
    ("IS5_ISPs1", ("IS5", "ISPs1")),          # 1 underscore
    ("IS3_ISPa_41", ("IS3", "ISPa_41")),      # 2 underscores
    ("_ISX", ("ISX", "ISX")),
])
def test_split_subject_id(subject, expected):
    assert h_is_stats.split_subject_id(subject) == expected


def test_stats_names_with_0_1_2_underscores(tmp_path):
    write_best(str(tmp_path / "bc1" / "best_alignment.tab"), [
        ("q1", "ISPs1", 100, "IS_DB"),
        ("q2", "IS5_ISPs1", 100, "IS_DB"),
        ("q3", "IS3_ISPa_41", 100, "IS_DB"),
    ])
    out = tmp_path / "stats.rcp"
    h_is_stats.main_is_stat(str(tmp_path), str(out))
    text = out.read_text()
    assert "ISPs1, ISPs1, 33.33, 1" in text
    assert "IS5, ISPs1, 33.33, 1" in text
    assert "IS3, ISPa_41, 33.33, 1" in text


# ---------------------------------------------------------------- B41
def test_counts_insertions_not_rows_and_reports_denominator(tmp_path):
    write_best(str(tmp_path / "bc1" / "best_alignment.tab"), [
        ("q1", "IS5_ISPs1", 500, "IS_DB"),
        ("q1", "IS5_ISPs1", 300, "IS_DB"),        # duplicate row of the same insertion
        ("q2", "IS3_ISPen2", 400, "IS_DB"),
        ("q2", "NC_000913", 900, "Genome"),       # same insertion also labelled Genome -> IS_DB wins
        ("q3", "NC_000913", 900, "Genome"),
        ("q4", "unknown", 10, "Contamination"),
        ("q5", "no blast hit", 0, "Unannotated"),
        ("q5", "no blast hit", 0, "Contamination"),  # Unannotated beats Contamination
        ("q6", "no blast hit", 0, "Unannotated"),
    ])
    s = h_is_stats.summarize_condition(str(tmp_path / "bc1" / "best_alignment.tab"))
    assert s["total"] == 6
    assert s["is_count"] == 2
    assert s["note_counts"] == {"IS_DB": 2, "Genome": 1, "Unannotated": 2, "Contamination": 1}
    assert sorted(s["rows"]) == [("IS3", "ISPen2", 50.0, 1), ("IS5", "ISPs1", 50.0, 1)]

    out = tmp_path / "stats.rcp"
    h_is_stats.main_is_stat(str(tmp_path), str(out))
    text = out.read_text()
    assert ">bc1\n?2\n#total_insertions=6\n#IS_DB=2\n#Genome=1\n#Unannotated=2\n#Contamination=1\n" in text
    assert "#percentage_denominator=IS_DB" in text


def test_stats_no_is_hits_and_empty_table(tmp_path):
    write_best(str(tmp_path / "a" / "best_alignment.tab"), [("q1", "chr", 10, "Genome")])
    write_best(str(tmp_path / "b" / "best_alignment.tab"), [])
    os.makedirs(tmp_path / "c")  # folder without a table is skipped
    out = tmp_path / "stats.rcp"
    h_is_stats.main_is_stat(str(tmp_path), str(out))
    text = out.read_text()
    assert ">a\n?0\n#total_insertions=1\n#IS_DB=0\n#Genome=1\n#Unannotated=0\n#Contamination=0\n" in text
    assert ">b\n?0\n#total_insertions=0\n" in text
    assert ">c" not in text


def test_stats_no_folders(tmp_path):
    out = tmp_path / "stats.rcp"
    h_is_stats.main_is_stat(str(tmp_path), str(out))
    assert out.read_text() == ""


# ---------------------------------------------------------------- B42
def write_rcp(path, conditions):
    segs = []
    for cond, rows in conditions.items():
        lines = [f">{cond}", f"?{sum(r[3] for r in rows)}", "#total_insertions=99"]
        lines += [f"{f}, {n}, {p:.2f}, {c}" for f, n, p, c in rows]
        segs.append("\n".join(lines) + "\n")
    path.write_text("---\n".join(segs))


@pytest.mark.parametrize("rows", [
    {"bc1": [("IS5", "ISPs1", 60.0, 6), ("IS3", "ISPen2", 40.0, 4)]},     # all high
    {"bc1": [("IS5", "ISPs1", 10.0, 1)], "bc2": [("IS3", "ISPen2", 5.0, 1)]},  # all low
    {"bc1": [("IS5", "ISPs1", 90.0, 9), ("IS3", "ISPen2", 10.0, 1)]},     # both panels
    {"bc1": []},                                                           # nothing
])
def test_heatmap_panels(tmp_path, rows):
    rcp = tmp_path / "s.rcp"
    write_rcp(rcp, rows)
    pdf, tsv = tmp_path / "h.pdf", tmp_path / "h.tsv"
    assert i_heatmap.main_heatmap(str(rcp), str(pdf), str(tsv)) == 1
    assert pdf.stat().st_size > 0
    assert tsv.exists()


def test_heatmap_reads_stats_output(tmp_path):
    write_best(str(tmp_path / "d" / "bc1" / "best_alignment.tab"),
               [("q1", "IS5_ISPs1", 1, "IS_DB"), ("q2", "IS3_ISPa_41", 1, "IS_DB"), ("q3", "g", 1, "Genome"),
                ("q4", "no blast hit", 0, "Unannotated"), ("q5", "no blast hit", 0, "Contamination")])
    rcp = tmp_path / "s.rcp"
    h_is_stats.main_is_stat(str(tmp_path / "d"), str(rcp))
    df, order = i_heatmap.parse_rcp(str(rcp))
    assert order == ["bc1"]
    assert sorted(df["Species"]) == ["IS3_ISPa_41", "IS5_ISPs1"]   # Unannotated/Contamination not plotted
    i_heatmap.main_heatmap(str(rcp), str(tmp_path / "h.pdf"), str(tmp_path / "h.tsv"))
    assert (tmp_path / "h.pdf").stat().st_size > 0


def test_heatmap_missing_file(tmp_path):
    assert i_heatmap.main_heatmap(str(tmp_path / "nope.rcp"), str(tmp_path / "h.pdf"), str(tmp_path / "h.tsv")) == 0


# ---------------------------------------------------------------- B43
def test_parse_insertion_tab_skips_header_and_in(tmp_path):
    tab = tmp_path / "insertion_from_bam.tab"
    tab.write_text(
        "read_id\tquery_start\tquery_end\tinsertion_length\tref_pos\tquality\tinsert_type\tis_reverse\n"
        "r1\t10\t1010\t1000\t500\t20\tSC\tFalse\n"
        "r2\t10\t1010\t1000\t500\t20\tSC\tTrue\n"
        "r3\t10\t510\t500\t700\t20\tIN\tFalse\n"
    )
    counts, weighted = f_in_del_plot.parse_insertion_tab_file(str(tab))
    assert counts == Counter({500: 2})
    assert weighted == Counter({500: 2000})


def test_counter_merge_adds_sc_and_in(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(f_in_del_plot, "parse_insertion_tab_file",
                        lambda p: (Counter({100: 2, 200: 1}), Counter({100: 2000, 200: 1000})))
    monkeypatch.setattr(f_in_del_plot, "extract_insertions_from_bam",
                        lambda p, t: (Counter({100: 3}), Counter({100: 30})))
    monkeypatch.setattr(f_in_del_plot, "extract_deletions_from_bam",
                        lambda p, t: (Counter(), Counter()))
    monkeypatch.setattr(f_in_del_plot, "plot_insertion_deletion_graphs",
                        lambda ins, dels, name, out: captured.__setitem__(name, ins))
    f_in_del_plot.main_in_del_plot("x.tab", str(tmp_path / "x.bam"), "o")
    assert captured["non-weighted"] == Counter({100: 5, 200: 1})
    assert captured["weighted"] == Counter({100: 2030, 200: 1000})
