#-------------------------------------
#This code used for visualization of read lengths

#INPUT: .fastq file (optionally gzipped)
#OUTPUT: Graphical outputs

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import gzip
import logging
import math
import os

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def _open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def calculate_read_lengths(fastq_file):
    """Return a list with the length of every read in a (4-line) FASTQ file."""
    read_lengths = []
    with _open_text(fastq_file) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                read_lengths.append(len(line.strip()))
    return read_lengths


def determine_bin_size(read_lengths):
    """Bin width from the square-root rule; always >= 1 (B48)."""
    if not read_lengths:
        return 1
    bin_count = math.ceil(math.sqrt(len(read_lengths)))
    data_range = max(read_lengths) - min(read_lengths)
    return max(1, math.ceil(data_range / bin_count))


def output_path(output, default_ext=".pdf"):
    """Append `default_ext` only when `output` has no extension (B48: no x.pdf.pdf)."""
    return output if os.path.splitext(output)[1] else output + default_ext


def plot_histogram(fastq_file, output):
    """Plot a read-length histogram of `fastq_file` and save it to `output`.

    `output` is used as-is if it has an extension, otherwise `.pdf` is appended.
    Returns the path of the written figure.
    """
    read_lengths = calculate_read_lengths(fastq_file)
    out = output_path(output)

    fig, ax = plt.subplots(figsize=(10, 6))
    try:
        if read_lengths:
            bin_size = determine_bin_size(read_lengths)
            bins = range(min(read_lengths), max(read_lengths) + bin_size + 1, bin_size)
            ax.hist(read_lengths, bins=bins, color='skyblue', edgecolor='black')
            ax.set_title(f"Read Length Distribution (n = {len(read_lengths)})", fontsize=16)
        else:
            logger.warning("No reads found in %s; writing an empty histogram.", fastq_file)
            ax.text(0.5, 0.5, "No reads", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("Read Length Distribution (n = 0)", fontsize=16)
        ax.set_xlabel("Read Length", fontsize=14)
        ax.set_ylabel("Frequency", fontsize=14)
        fig.savefig(out)
    finally:
        plt.close(fig)
    logger.info("Figure saved as %s", out)
    return out
