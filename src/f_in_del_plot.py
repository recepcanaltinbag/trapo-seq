#-------------------------------------
#This code used for visualization of insertions and deletions based on reference positions
#Also, This code can deal with fragmented alignments, insertions can be taken from tabular file.

#INPUT: BAM FILE, and TAB FILE of insertions (insertion_from_bam.tab):
#  read_id, query_start, query_end, insertion_length, ref_pos, quality, insert_type, is_reverse
#  (a header line, if present, is skipped)
#OUTPUT: <output>_{non-weighted,weighted}.pdf and
#        <output>_{non-weighted,weighted}_{insertion,deletion}_table.tsv

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import logging
import os
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd
import pysam

logger = logging.getLogger(__name__)

# CIGAR operations
CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_P, CIGAR_EQ, CIGAR_X = range(9)
REF_CONSUMING = {CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X}


def parse_insertion_tab_file(file_path):
    """Count soft-clip (non-IN) insertions per reference position from the tab file.

    Returns (insertion_counts, weighted_insertion_counts) as Counters.
    """
    insertion_counts = Counter()
    weighted_insertion_counts = Counter()

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            try:
                insertion_length = int(parts[3])
                ref_pos = int(parts[4])
            except ValueError:  # header line
                continue
            insert_type = parts[6]
            if insert_type == "IN":  # IN insertions are taken directly from the BAM file
                continue
            insertion_counts[ref_pos] += 1
            weighted_insertion_counts[ref_pos] += insertion_length

    return insertion_counts, weighted_insertion_counts


def _usable(read):
    return not (read.is_unmapped or read.is_secondary) and read.cigartuples


def extract_deletions_from_bam(bam_file_path, deletion_threshold):
    """Deletions (D) >= threshold: count at the deletion start; weighted count spread over
    every deleted base (weight = deletion length)."""
    deletion_counts = Counter()
    weighted_deletion_counts = Counter()

    with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if not _usable(read):
                continue
            ref_pos = read.reference_start
            for operation, length in read.cigartuples:
                if operation == CIGAR_D and length >= deletion_threshold:
                    deletion_counts[ref_pos] += 1
                    for offset in range(length):
                        weighted_deletion_counts[ref_pos + offset] += length
                if operation in REF_CONSUMING:
                    ref_pos += length

    return deletion_counts, weighted_deletion_counts


def extract_insertions_from_bam(bam_file_path, insertion_threshold):
    """CIGAR insertions (I) >= threshold, counted at their reference position."""
    insertion_counts = Counter()
    weighted_insertion_counts = Counter()

    with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if not _usable(read):
                continue
            ref_pos = read.reference_start
            for operation, length in read.cigartuples:
                if operation == CIGAR_I and length >= insertion_threshold:
                    insertion_counts[ref_pos] += 1
                    weighted_insertion_counts[ref_pos] += length
                if operation in REF_CONSUMING:
                    ref_pos += length

    return insertion_counts, weighted_insertion_counts


#Plotting
def plot_insertion_deletion_graphs(insertion_counts, deletion_counts, name, output):
    logger.info("%s_%s is processing", output, name)
    sorted_insertion_positions = sorted(insertion_counts.keys())
    sorted_insertion_counts = [insertion_counts[pos] for pos in sorted_insertion_positions]

    sorted_deletion_positions = sorted(deletion_counts.keys())
    sorted_deletion_counts = [deletion_counts[pos] for pos in sorted_deletion_positions]

    fig, ax1 = plt.subplots(figsize=(10, 6))
    try:
        ax1.set_xlabel('Reference Position')
        ax1.set_ylabel(f"Deletions [{name}]")
        ax1.plot(sorted_deletion_positions, sorted_deletion_counts, label=f"Deletions [{name}]",
                 color='m', marker='x')
        ax1.legend(loc='upper left')
        ax1.grid(True)
        ax1.minorticks_on()

        ax2 = ax1.twinx()
        ax2.set_ylabel(f"Insertions [{name}]")
        ax2.plot(sorted_insertion_positions, sorted_insertion_counts, label=f"Insertions [{name}]",
                 color='g', marker='s')
        ax2.legend(loc='upper right')
        ax2.minorticks_on()

        ax1.set_title('Insertions and Deletions')
        fig.savefig(f"{output}_{name}.pdf")
    finally:
        plt.close(fig)
    logger.info("%s_%s.pdf saved", output, name)

    df_ins = pd.DataFrame({"Insertion Position": sorted_insertion_positions,
                           "Insertion Count": sorted_insertion_counts})
    df_del = pd.DataFrame({"Deletion Position": sorted_deletion_positions,
                           "Deletion Count": sorted_deletion_counts})
    df_ins.to_csv(f'{output}_{name}_insertion_table.tsv', sep='\t', index=False)
    df_del.to_csv(f'{output}_{name}_deletion_table.tsv', sep='\t', index=False)


def main_in_del_plot(tab_file_path, bam_file_path, output_name, insertion_threshold=2, deletion_threshold=5):
    """Outputs are written next to the BAM file: <bam dir>/<output_name>_*."""
    dir_path = os.path.dirname(os.path.abspath(bam_file_path))
    output = os.path.join(dir_path, output_name)

    # SC insertions from the tab file + IN insertions from the BAM file (B43: add, not overwrite)
    insertion_counts, weighted_insertion_counts = parse_insertion_tab_file(tab_file_path)
    insertion_counts_in, weighted_insertion_counts_in = extract_insertions_from_bam(bam_file_path, insertion_threshold)
    insertion_counts = insertion_counts + insertion_counts_in
    weighted_insertion_counts = weighted_insertion_counts + weighted_insertion_counts_in

    deletion_counts, weighted_deletion_counts = extract_deletions_from_bam(bam_file_path, deletion_threshold)

    plot_insertion_deletion_graphs(insertion_counts, deletion_counts, 'non-weighted', output)
    plot_insertion_deletion_graphs(weighted_insertion_counts, weighted_deletion_counts, 'weighted', output)
    print(f"In/del plots and tables written to {output}_*")
