"""DR logo module: alignment trimming, MAFFT handling, empty inputs."""
import os
import shutil

import pandas as pd
import pytest

from src import k_dr_logo
from src.k_dr_logo import align_sequences, main_dr_logo, remove_high_gap_columns

requires_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not in PATH")


def test_gap_filter_keeps_columns_aligned():
    aligned = ["-ACGT", "-ACGT", "TACGT", "-ACG-"]
    # column 0 is 75 % gaps, column 4 is 25 %: only column 0 is removed at 60 %
    assert remove_high_gap_columns(aligned, 60) == ["ACGT", "ACGT", "ACGT", "ACG-"]
    # threshold is inclusive
    assert remove_high_gap_columns(aligned, 75) == aligned
    assert remove_high_gap_columns([], 60) == []


def test_equal_length_needs_no_mafft(monkeypatch):
    monkeypatch.setattr(k_dr_logo.shutil, "which", lambda _: None)
    assert align_sequences(["acgta", "ACGTT", "N/A", None, "AAGTA"]) == ["ACGTA", "ACGTT", "AAGTA"]


def test_missing_mafft_clear_error(monkeypatch):
    monkeypatch.setattr(k_dr_logo.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="mafft"):
        align_sequences(["ACGTA", "ACGTAC"])


@requires_mafft
def test_mafft_mixed_lengths_columns_consistent():
    seqs = ["GATTACA", "GATTACA", "GATTACA", "ATTACA", "GATTACAT"]
    aligned = align_sequences(seqs, 60)
    assert len({len(s) for s in aligned}) == 1
    # the dominant motif must be recovered intact in every full-length row
    assert aligned[0] == "GATTACA"
    assert sum(s == "GATTACA" for s in aligned) >= 3


def _plasmid(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(">p\n" + "ACGT" * 500 + "\n")
    return str(p)


def test_empty_input_folder(tmp_path):
    ins = tmp_path / "ins"
    ins.mkdir()
    assert main_dr_logo(_plasmid(tmp_path), str(ins), str(tmp_path / "out")) == []
    # header-only file is also handled
    pd.DataFrame(columns=["Read ID", "Subject ID", "Insertion Point", "Repeat Length",
                          "Overlap Sequence", "Status"]).to_csv(ins / "a_insertions.csv", index=False)
    assert main_dr_logo(_plasmid(tmp_path), str(ins), str(tmp_path / "out")) == []


def test_logo_and_histogram_written(tmp_path):
    ins = tmp_path / "ins dir"
    ins.mkdir()
    pd.DataFrame({
        "Read ID": ["r1", "r2", "r3", "r4", "r5"],
        "Subject ID": ["IS4_x", "IS4_x", "IS4_x", "IS4_x", "IS5/y"],
        "Insertion Point": [10, 500, 501, None, 7],
        "Repeat Length": [4, 4, 4, None, 0],
        "Overlap Sequence": ["ACGT", "ACGA", "TCGT", "N/A", "N/A"],
        "Status": ["DR", "DR", "DR", "no_flank_hit", "no_duplication"],
    }).to_csv(ins / "s1_insertions.csv", index=False)
    out = tmp_path / "out"
    written = main_dr_logo(_plasmid(tmp_path), str(ins), str(out))
    names = sorted(os.path.basename(w) for w in written)
    assert names == ["IS4_x_insertion_points_histogram.pdf", "IS4_x_overlap_sequence_logo.pdf",
                     "IS5_y_insertion_points_histogram.pdf"]
    assert all(os.path.getsize(w) > 0 for w in written)
    import matplotlib.pyplot as plt
    assert plt.get_fignums() == []          # no figures left open
