#-------------------------------------
#This code used for stats with visualization of mapped reads to the plasmid and the genome
#It can be used to understand which part of the genome is mobile and how many reads cover the plasmid
#Restriction sites and deletions from the plasmid can be seen

#INPUT: BAM FILES of Genome and plasmid maps (the genome BAM is the 07_*.bam next to the plasmid BAM)
#OUTPUT: <output>.pdf (coverage plot) and <output>_mapping_stats.txt

#Coverage is computed per contig (B45). Contigs are drawn one after another in BAM
#header order, separated by dotted lines.

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import logging
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pysam

logger = logging.getLogger(__name__)

MAX_PLOT_POINTS = 20000


def count_reads_in_fastq(fastq_file):
    """Number of reads in a 4-line FASTQ file."""
    with open(fastq_file, 'r') as file:
        return sum(1 for i, _ in enumerate(file) if i % 4 == 0)


def _ensure_index(bam_file):
    if not (os.path.exists(bam_file + ".bai") or os.path.exists(bam_file + ".csi")
            or os.path.exists(os.path.splitext(bam_file)[0] + ".bai")):
        logger.info("Index file for %s not found, creating index...", bam_file)
        pysam.index(bam_file)


def read_alignment(bam_file):
    """Return {read_name: [(contig, ref_start, ref_end), ...]} for primary and
    supplementary alignments (secondary and unmapped records are skipped)."""
    _ensure_index(bam_file)
    reads = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as alignment:
        for read in alignment.fetch():
            if read.is_unmapped or read.is_secondary or read.reference_end is None:
                continue
            reads[read.query_name].append((read.reference_name, read.reference_start, read.reference_end))
    return reads


def get_contig_lengths(bam_file_path):
    """Ordered {contig: length} from the BAM header."""
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        return dict(zip(bam_file.references, bam_file.lengths))


def get_reference_lengths(bam_file_path):
    return list(get_contig_lengths(bam_file_path).values())


def compute_coverage(reads, contig_lengths):
    """Per-contig coverage arrays {contig: np.ndarray}, built with a difference array."""
    diffs = {c: np.zeros(n + 1, dtype=np.int64) for c, n in contig_lengths.items()}
    starts = defaultdict(list)
    ends = defaultdict(list)
    for alignments in reads.values():
        for contig, start, end in alignments:
            if contig in diffs:
                starts[contig].append(start)
                ends[contig].append(end)
    coverage = {}
    for contig, n in contig_lengths.items():
        diff = diffs[contig]
        if starts[contig]:
            s = np.clip(np.asarray(starts[contig]), 0, n)
            e = np.clip(np.asarray(ends[contig]), 0, n)
            np.add.at(diff, s, 1)
            np.add.at(diff, e, -1)
        coverage[contig] = np.cumsum(diff[:-1])
    return coverage


def _downsample(x, y, max_points=MAX_PLOT_POINTS):
    """Max-pool long tracks so the PDF stays small; peaks are preserved."""
    if len(y) <= max_points:
        return x, y
    step = int(np.ceil(len(y) / max_points))
    n = len(y) // step * step
    y_ds = y[:n].reshape(-1, step).max(axis=1)
    x_ds = x[:n:step]
    if n < len(y):
        y_ds = np.append(y_ds, y[n:].max())
        x_ds = np.append(x_ds, x[n])
    return x_ds, y_ds


def _plot_track(ax, coverage, color, label):
    offset = 0
    for i, (contig, cov) in enumerate(coverage.items()):
        x = np.arange(offset, offset + len(cov))
        xs, ys = _downsample(x, cov)
        ax.plot(xs, ys, color=color, label=label if i == 0 else None)
        if i > 0:
            ax.axvline(offset, color="grey", linestyle=":", linewidth=0.8)
        offset += len(cov)
    ax.minorticks_on()


def plot_read_distribution(plasmid_reads, genome_reads, plasmid_lengths, genome_lengths, output):
    """`plasmid_lengths` / `genome_lengths` are ordered {contig: length} dicts."""
    plasmid_cov = compute_coverage(plasmid_reads, plasmid_lengths)
    genome_cov = compute_coverage(genome_reads, genome_lengths)

    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    try:
        _plot_track(ax[0], plasmid_cov, "C0", "Plasmid")
        ax[0].set_title("Read Distribution on Plasmid")
        ax[0].set_xlabel("Plasmid Position" + (" (contigs concatenated)" if len(plasmid_cov) > 1 else ""))
        ax[0].set_ylabel("Read Count")

        _plot_track(ax[1], genome_cov, "red", "Other Genome")
        ax[1].set_title("Read Distribution on Genome")
        ax[1].set_xlabel("Genomic Position" + (" (contigs concatenated in BAM header order)"
                                               if len(genome_cov) > 1 else ""))
        ax[1].set_ylabel("Read Count")

        fig.tight_layout()
        fig.savefig(f"{output}.pdf")
    finally:
        plt.close(fig)
    logger.info("Figure saved as %s.pdf", output)
    return plasmid_cov, genome_cov


def summarize_mappings(plasmid_reads, genome_reads, output_txt, fastq_file=None):
    """Write read-level mapping statistics. The FASTQ total is reported only if given."""
    plasmid_ids = set(plasmid_reads)
    genome_ids = set(genome_reads)
    total_mapped = plasmid_ids | genome_ids
    plasmid_only = plasmid_ids - genome_ids
    both = plasmid_ids & genome_ids
    denom = len(both) + len(plasmid_only)
    rate = f"{len(both) / denom:.2f}" if denom else "NA"

    lines = []
    if fastq_file:
        lines.append(f"Total reads: {count_reads_in_fastq(fastq_file)}")
    lines += [f"Total mapped reads: {len(total_mapped)}",
              f"Plasmid only: {len(plasmid_only)}",
              f"Genome only: {len(genome_ids - plasmid_ids)}",
              f"Both: {len(both)}",
              f"Transposition Rate (both / plasmid-mapped): {rate}"]
    with open(output_txt, 'w') as f:
        f.write("\n".join(lines) + "\n")
    for line in lines:
        logger.info(line)
    logger.info("Mapping statistics written to %s", output_txt)
    return {"total_mapped": len(total_mapped), "plasmid_only": len(plasmid_only),
            "both": len(both), "rate": rate}


def find_file_with_prefix(directory, prefix, enfix):
    """First file (sorted by name) in `directory` that starts with `prefix` and ends with `enfix`."""
    for filename in sorted(os.listdir(directory)):
        if filename.startswith(prefix) and filename.endswith(enfix):
            return os.path.join(directory, filename)
    return None


def main_map_dist(plasmid_bam, output_p, fastq_file=None):
    dir_path = os.path.dirname(os.path.abspath(plasmid_bam))
    output = os.path.join(dir_path, output_p)

    genome_bam = find_file_with_prefix(dir_path, "07_", ".bam")
    if genome_bam is None:
        raise FileNotFoundError(f"No genome BAM (07_*.bam) found next to {plasmid_bam}")
    logger.info("Corresponding genome bam: %s", genome_bam)

    plasmid_lengths = get_contig_lengths(plasmid_bam)
    genome_lengths = get_contig_lengths(genome_bam)
    logger.info("Plasmid length: %d (%d contig(s)); genome length: %d (%d contig(s))",
                sum(plasmid_lengths.values()), len(plasmid_lengths),
                sum(genome_lengths.values()), len(genome_lengths))

    plasmid_reads = read_alignment(plasmid_bam)
    genome_reads = read_alignment(genome_bam)

    plot_read_distribution(plasmid_reads, genome_reads, plasmid_lengths, genome_lengths, output)
    summarize_mappings(plasmid_reads, genome_reads, f"{output}_mapping_stats.txt", fastq_file)
    print(f"Read distribution written to {output}.pdf and {output}_mapping_stats.txt")
