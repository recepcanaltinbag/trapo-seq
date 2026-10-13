#-------------------------------------
#This code used for kde visualization of how many basepairs in reads come from genome and plasmid
#It can be used to understand how much length of insertions
#Read coordinates of genome alignments are taken in original-read orientation, including hard clips (B44)

#INPUT: BAM FILES of Genome and plasmid maps
#OUTPUT: Graphical outputs

#Requirements: matplotlib, seaborn

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------
import logging
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pysam
import seaborn as sns

logger = logging.getLogger(__name__)

CIGAR_QUERY_ALIGNED = {0, 1, 7, 8}  # M, I, =, X
CIGAR_CLIP = {4, 5}                 # S, H


def _ensure_index(bam_file):
    if not (os.path.exists(bam_file + ".bai") or os.path.exists(bam_file + ".csi")
            or os.path.exists(os.path.splitext(bam_file)[0] + ".bai")):
        logger.info("Index file for %s not found, creating index...", bam_file)
        pysam.index(bam_file)


def read_query_interval(read):
    """Aligned interval of `read` in original-read coordinates (B44).

    Hard clips are counted (supplementary records usually hard-clip the rest of the
    read) and the interval is flipped for reverse-strand alignments, so every record of
    one read uses the same 5'->3' coordinate system. Returns (start, end, read_length).
    """
    cigar = read.cigartuples
    if not cigar:
        return None
    lead = 0
    for op, length in cigar:
        if op in CIGAR_CLIP:
            lead += length
        else:
            break
    aligned = sum(length for op, length in cigar if op in CIGAR_QUERY_ALIGNED)
    total = sum(length for op, length in cigar if op in CIGAR_QUERY_ALIGNED or op in CIGAR_CLIP)
    start, end = lead, lead + aligned
    if read.is_reverse:
        start, end = total - end, total - start
    return start, end, total


def read_alignment(bam_file, genome_name):
    _ensure_index(bam_file)

    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)
    for read in alignment.fetch():
        if not read.is_unmapped and not read.is_secondary:
            reads[read.query_name].append((genome_name, read.reference_start, read.reference_end))
    alignment.close()
    return reads


def calculate_overlap(start1, end1, start2, end2):
    """
    Calculate the overlap percentage between two genomic regions.
    
    Parameters:
    start1, end1: Start and end positions of the first region.
    start2, end2: Start and end positions of the second region.
    
    Returns:
    float: The overlap percentage (0-100%) between the two regions.
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    # No overlap
    if overlap_start >= overlap_end:
        return 0.0
    
    # Calculate the overlap length
    overlap_length = overlap_end - overlap_start
    region1_length = end1 - start1
    region2_length = end2 - start2
    
    # Return the minimum overlap as a percentage of the total region length
    return max(overlap_length / region1_length, overlap_length / region2_length) * 100

def read_alignment_plasmid(bam_file, genome_name, overlap_threshold=80.0):
    """
    Reads alignments from a BAM file, excluding redundant alignments with overlap 
    greater than the specified threshold.

    Parameters:
    bam_file (str): Path to the BAM file.
    genome_name (str): Name of the genome being aligned to.
    overlap_threshold (float): Maximum allowable overlap percentage (default is 80%).

    Returns:
    dict: A dictionary with read names as keys and lists of non-redundant alignments as values.
    """
    # Ensure BAM index exists
    _ensure_index(bam_file)

    # Open BAM file
    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)

    for read in alignment.fetch():
        if not read.is_unmapped and not read.is_secondary and read.reference_end is not None:
            query_name = read.query_name
            ref_start, ref_end = read.reference_start, read.reference_end
            is_redundant = False

            # Check for significant overlap with existing alignments for this read
            for existing_genome, existing_start, existing_end in reads[query_name]:
                overlap_percentage = calculate_overlap(ref_start, ref_end, existing_start, existing_end)
                
                # If overlap exceeds the threshold, skip this alignment
                if overlap_percentage > overlap_threshold:
                    is_redundant = True
                    break

            # Only add alignment if it's not redundant
            if not is_redundant:
                reads[query_name].append((genome_name, ref_start, ref_end))
                
    alignment.close()
    return reads




def calculate_overlap_qstart_end(q_start1, q_end1, q_start2, q_end2):
    """
    Calculate the overlap percentage between two query regions based on query start/end.
    
    Parameters:
    q_start1, q_end1: Start and end positions of the first query region.
    q_start2, q_end2: Start and end positions of the second query region.
    
    Returns:
    float: The overlap percentage (0-100%) between the two query regions.
    """
    overlap_start = max(q_start1, q_start2)
    overlap_end = min(q_end1, q_end2)
    
    # No overlap
    if overlap_start >= overlap_end:
        return 0.0
    
    # Calculate the overlap length
    overlap_length = overlap_end - overlap_start
    query1_length = q_end1 - q_start1
    query2_length = q_end2 - q_start2
    
    # Return the minimum overlap as a percentage of the total query length
    return max(overlap_length / query1_length, overlap_length / query2_length) * 100

def read_unique_alignment_genome(bam_file, genome_name, overlap_threshold=80.0):
    """
    Reads alignments from a BAM file, excluding redundant alignments with overlap 
    greater than the specified threshold based on query start and end.

    Parameters:
    bam_file (str): Path to the BAM file.
    genome_name (str): Name of the genome being aligned to.
    overlap_threshold (float): Maximum allowable overlap percentage (default is 80%).

    Returns:
    dict: A dictionary with read names as keys and lists of non-redundant alignments as values.
    """
    # Ensure BAM index exists
    _ensure_index(bam_file)

    # Open BAM file
    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)

    for read in alignment.fetch():
        if read.is_unmapped or read.is_secondary:
            continue
        interval = read_query_interval(read)
        if interval is None or interval[1] <= interval[0]:
            continue
        query_name = read.query_name
        q_start, q_end, _ = interval
        is_redundant = False
        # Check for significant overlap with existing alignments for this read
        for existing_genome, existing_q_start, existing_q_end in reads[query_name]:
            overlap_percentage = calculate_overlap_qstart_end(q_start, q_end, existing_q_start, existing_q_end)
            if overlap_percentage > overlap_threshold:
                is_redundant = True
                break

        # Only add alignment if it's not redundant
        if not is_redundant:
            reads[query_name].append((genome_name, q_start, q_end))
                
    alignment.close()
    return reads


def merge_intervals(intervals):
    """Verilen interval listesinde örtüşenleri birleştirir."""
    if not intervals:
        return []
    
    # Intervals'ı start'a göre sıralıyoruz
    intervals.sort(key=lambda x: x[1])
    
    # İlk bölgeyi başlat
    merged = [intervals[0]]
    
    for current in intervals[1:]:
        prev_genome, prev_start, prev_end = merged[-1]
        _, curr_start, curr_end = current
        
        # Eğer mevcut bölge öncekiyle örtüşüyorsa
        if curr_start <= prev_end:
            # Örtüşen bölgeleri birleştir
            merged[-1] = (prev_genome, prev_start, max(prev_end, curr_end))
        else:
            # Eğer örtüşmüyorsa, yeni bir bölge ekle
            merged.append(current)
    
    return merged

def analyze_basepair_distribution(plasmid_reads, genome_reads):
    read_distribution = defaultdict(lambda: {"plasmid": 0, "genome": 0})
    
    for read_name, alignments in plasmid_reads.items():
        merged_plasmid_alignments = merge_intervals(alignments)
        read_distribution[read_name]["plasmid"] += sum(end - start for _, start, end in merged_plasmid_alignments)

    for read_name, alignments in genome_reads.items():
        merged_genome_alignments = merge_intervals(alignments)
        read_distribution[read_name]["genome"] += sum(end - start for _, start, end in merged_genome_alignments)
    
    return read_distribution


def write_read_distribution_to_file(read_distribution, file_path):
    # Open the file in write mode
    with open(file_path + '.tab', 'w') as file:
        # Write the header to the file (optional)
        file.write("Read ID\tPlasmid Base Pairs\tGenome Base Pairs\n")
        
        # Iterate through the dictionary and write each read's information
        for read_id, base_pairs_info in read_distribution.items():
            plasmid_bps = base_pairs_info.get('plasmid', 0)  # Default to 0 if key does not exist
            genome_bps = base_pairs_info.get('genome', 0)    # Default to 0 if key does not exist
            
            # Write each entry in tab-separated format
            file.write(f"{read_id}\t{plasmid_bps}\t{genome_bps}\n")


def _filter_by_plasmid(read_distribution, plasmid_threshold):
    pairs = [(info["plasmid"], info["genome"]) for info in read_distribution.values()
             if info["plasmid"] > plasmid_threshold]
    if not pairs:
        return np.array([]), np.array([])
    p, g = zip(*pairs)
    return np.asarray(p), np.asarray(g)


def plot_read_distribution(read_distribution, output, plasmid_threshold):
    """KDE (<output>.pdf) and KDE + scatter (<output>_scatter.pdf) of plasmid vs genome
    base pairs per read, for reads with more than `plasmid_threshold` plasmid bp."""
    plasmid_bps, genome_bps = _filter_by_plasmid(read_distribution, plasmid_threshold)
    n = len(plasmid_bps)

    fig, ax = plt.subplots(figsize=(12, 6))
    try:
        if n == 0:  # B46
            logger.warning("No read has more than %d plasmid bp; writing empty KDE plots.", plasmid_threshold)
            ax.text(0.5, 0.5, f"No read with > {plasmid_threshold} plasmid bp", ha="center",
                    va="center", transform=ax.transAxes)
        elif n < 3 or np.ptp(plasmid_bps) == 0 or np.ptp(genome_bps) == 0:
            logger.warning("Too few / degenerate points (%d reads) for a 2D KDE; only the scatter is drawn.", n)
        else:
            try:
                sns.kdeplot(x=plasmid_bps, y=genome_bps, cmap="coolwarm", fill=True, thresh=0.01, ax=ax)
            except (ValueError, np.linalg.LinAlgError) as err:
                logger.warning("KDE failed (%s); only the scatter is drawn.", err)
        ax.minorticks_on()
        ax.set_xlabel("Base pairs from plasmid")
        ax.set_ylabel("Base pairs from genome")
        ax.set_title(f"Read Base Pair Distribution between Plasmid and Genome (KDE, n = {n})")
        fig.savefig(f"{output}.pdf")
        if n:
            ax.scatter(plasmid_bps, genome_bps, color='black', alpha=0.5, s=1, label="Reads")
        fig.savefig(f"{output}_scatter.pdf")
    finally:
        plt.close(fig)
    logger.info("Figures saved as %s.pdf and %s_scatter.pdf", output, output)
    return n


#Get reference lengths from bam file
def get_reference_lengths(bam_file_path):
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        return list(bam_file.lengths)


# Hotspot detection (identify_hotspots / plot_hotspots) was removed: it was never
# called and treated read (query) coordinates as genome coordinates. Genome coverage
# per contig is plotted by c_map_dist.py.


def find_file_with_prefix(directory, prefix, enfix):
    """First file (sorted by name) in `directory` that starts with `prefix` and ends with `enfix`."""
    for filename in sorted(os.listdir(directory)):
        if filename.startswith(prefix) and filename.endswith(enfix):
            return os.path.join(directory, filename)
    return None


def main_kde_mobile(plasmid_bam, output_kde, plasmid_threshold):
    """Writes <bam dir>/<output_kde>.tab, .pdf and _scatter.pdf."""
    dir_path = os.path.dirname(os.path.abspath(plasmid_bam))
    output = os.path.join(dir_path, output_kde)

    genome_bam = find_file_with_prefix(dir_path, "07_", ".bam")
    if genome_bam is None:
        raise FileNotFoundError(f"No genome BAM (07_*.bam) found next to {plasmid_bam}")
    logger.info("Corresponding genome bam: %s", genome_bam)

    plasmid_reads = read_alignment_plasmid(plasmid_bam, "plasmid")
    genome_reads = read_unique_alignment_genome(genome_bam, "other_genome")

    read_distribution = analyze_basepair_distribution(plasmid_reads, genome_reads)

    write_read_distribution_to_file(read_distribution, output)
    plot_read_distribution(read_distribution, output, plasmid_threshold)
    print(f"KDE outputs written to {output}.tab, {output}.pdf, {output}_scatter.pdf")
