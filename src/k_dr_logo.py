"""Sequence logos of direct repeats (TSDs) and insertion-point histograms.

INPUT:  ``<input_folder>/*_insertions.csv`` written by ``dr_finder``.
OUTPUT: per Subject ID (IS name) in ``output_folder``:
  ``<IS>_insertion_points_histogram.pdf`` and, when at least two DRs exist,
  ``<IS>_overlap_sequence_logo.pdf``.

DR sequences are reported on the plasmid + strand starting at the first base of
the target site, so equal-length DRs are already aligned and are used as is.
Mixed lengths are aligned with MAFFT; afterwards only columns with at most
``gap_threshold`` percent gaps are kept (columns are never shifted or truncated
from one end, so positions stay aligned).
"""

import logging
import os
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import logomaker  # noqa: E402
import pandas as pd  # noqa: E402
from Bio import SeqIO  # noqa: E402

logger = logging.getLogger(__name__)

COLORS = {"A": "green", "C": "blue", "G": "orange", "T": "red"}


def get_sequence_length(fasta_file):
    """Length of the first sequence in a FASTA file."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        return len(record.seq)
    raise ValueError(f"No sequence found in {fasta_file}")


def clean_filename(name):
    """Make an IS name safe for use in a file name."""
    for ch in ("::", "/", "\\", " ", ":"):
        name = name.replace(ch, "_")
    return name


def run_mafft(input_file, output_file):
    """Align ``input_file`` with MAFFT into ``output_file``; raise RuntimeError on failure."""
    if shutil.which("mafft") is None:
        raise RuntimeError("'mafft' not found in PATH; it is required by dr_logo for DRs of different lengths.")
    with open(output_file, "w") as out:
        res = subprocess.run(["mafft", "--auto", "--nuc", "--quiet", input_file],
                             stdout=out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"MAFFT failed (exit {res.returncode}): {res.stderr.strip()}")


def remove_high_gap_columns(aligned, gap_threshold):
    """Keep only columns whose gap percentage is <= gap_threshold.

    aligned: list of equal-length strings. Returns a list of strings.
    """
    if not aligned:
        return []
    n = len(aligned)
    keep = [i for i in range(len(aligned[0]))
            if 100.0 * sum(s[i] == "-" for s in aligned) / n <= gap_threshold]
    return ["".join(s[i] for i in keep) for s in aligned]


def align_sequences(sequences, gap_threshold=60):
    """Return gap-filtered, column-consistent aligned sequences (list of str)."""
    seqs = [str(s).upper() for s in sequences if isinstance(s, str) and s and s != "N/A"]
    if len(seqs) < 2:
        return seqs
    if len({len(s) for s in seqs}) == 1:
        return seqs                      # anchored at the target-site start: already aligned
    with tempfile.TemporaryDirectory(prefix="trapo_logo_") as tmp:
        fin, fout = os.path.join(tmp, "in.fasta"), os.path.join(tmp, "out.fasta")
        with open(fin, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">s{i}\n{s}\n")
        run_mafft(fin, fout)
        records = {r.id: str(r.seq).upper() for r in SeqIO.parse(fout, "fasta")}
    aligned = [records[f"s{i}"] for i in range(len(seqs))]
    return remove_high_gap_columns(aligned, gap_threshold)


def create_sequence_logo(aligned, output_file, title):
    counts = logomaker.alignment_to_matrix(sequences=aligned, to_type="counts",
                                           characters_to_ignore=".-N")
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * len(counts) + 2), 3))
    try:
        logo = logomaker.Logo(counts, ax=ax, color_scheme=COLORS)
        logo.style_spines(visible=False)
        logo.style_spines(spines=["left", "bottom"], visible=True)
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels([str(i + 1) for i in range(len(counts))])
        ax.set_ylabel("count")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(output_file, dpi=300)
    finally:
        plt.close(fig)


def plot_insertion_histogram(points, plasmid_len, output_file, title):
    fig, ax = plt.subplots(figsize=(10, 6))
    try:
        ax.hist(points, bins=60, range=(1, plasmid_len), histtype="stepfilled",
                edgecolor="orange", linewidth=1.5, color="white")
        ax.set_xlim(0, plasmid_len)
        ax.set_title(title)
        ax.set_xlabel("Insertion point (plasmid bp)")
        ax.set_ylabel("Frequency")
        fig.savefig(output_file)
    finally:
        plt.close(fig)


def _dr_sequences(group):
    if "Status" in group.columns:
        mask = group["Status"] == "DR"
    else:
        mask = pd.to_numeric(group["Repeat Length"], errors="coerce") > 0
    seqs = group.loc[mask, "Overlap Sequence"]
    return [s for s in seqs if isinstance(s, str) and s and s != "N/A"]


def main_dr_logo(plasmid_fasta, input_folder, output_folder, gap_threshold=60):
    """Histogram + DR logo per IS. Returns the list of written files."""
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    plasmid_len = get_sequence_length(plasmid_fasta)
    os.makedirs(output_folder, exist_ok=True)

    files = sorted(f for f in os.listdir(input_folder)
                   if f.endswith("_insertions.csv") and not f.startswith("."))
    frames = []
    for f in files:
        try:
            frames.append(pd.read_csv(os.path.join(input_folder, f),
                                      dtype={"Read ID": str, "Subject ID": str, "Overlap Sequence": str}))
        except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError) as e:
            logger.warning("Could not read %s: %s", f, e)
    frames = [df for df in frames if not df.empty]
    if not frames:
        logger.warning("No insertions found in %s (*_insertions.csv); nothing to plot", input_folder)
        return []
    data = pd.concat(frames, ignore_index=True)

    written = []
    for name, group in data.dropna(subset=["Subject ID"]).groupby("Subject ID"):
        safe = clean_filename(str(name))
        points = pd.to_numeric(group["Insertion Point"], errors="coerce").dropna()
        if not points.empty:
            out = os.path.join(output_folder, f"{safe}_insertion_points_histogram.pdf")
            plot_insertion_histogram(points, plasmid_len, out, f"Insertion points: {name} (n={len(points)})")
            written.append(out)
        seqs = _dr_sequences(group)
        if len(seqs) < 2:
            logger.info("%s: %d DR sequence(s), no logo", name, len(seqs))
            continue
        aligned = align_sequences(seqs, gap_threshold)
        if not aligned or not aligned[0]:
            logger.info("%s: no alignment columns left after gap filtering, no logo", name)
            continue
        out = os.path.join(output_folder, f"{safe}_overlap_sequence_logo.pdf")
        create_sequence_logo(aligned, out, f"Direct repeats: {name} (n={len(seqs)})")
        written.append(out)
    return written
