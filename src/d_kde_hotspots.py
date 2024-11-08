#-------------------------------------
#This code used for kde visualization of how many basepairs in reads come from genome and plasmid
#It can be used to understand how much length of insertions
#Also shows the hotspots positions of genome

#INPUT: BAM FILES of Genome and plasmid maps
#OUTPUT: Graphical outputs

#Requirements: matplotlib, seaborn

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------
import pysam
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import os
import seaborn as sns


def read_alignment(bam_file, genome_name):
    if not os.path.exists(bam_file + ".bai"):
        print(f"Index file for {bam_file} not found, creating index...")
        pysam.index(bam_file)

    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)
    for read in alignment.fetch():
        if not read.is_unmapped:
            reads[read.query_name].append((genome_name, read.reference_start, read.reference_end))
    alignment.close()
    return reads


def analyze_basepair_distribution(plasmid_reads, genome_reads):
    read_distribution = defaultdict(lambda: {"plasmid": 0, "genome": 0})
    
    for read_name, alignments in plasmid_reads.items():
        for genome_name, start, end in alignments:
            read_distribution[read_name]["plasmid"] += (end - start)
    
    for read_name, alignments in genome_reads.items():
        for genome_name, start, end in alignments:
            read_distribution[read_name]["genome"] += (end - start)
    
    return read_distribution


def plot_read_distribution(read_distribution, output):
    """
    Plot the read base pair distribution between plasmid and genome as a KDE plot.
    
    Args:
        read_distribution (dict): A dictionary where each key is a read ID and 
                                  the value is a dictionary with 'plasmid' and 'genome' base pair information.
    """
    # Extract plasmid and genome base pair information
    plasmid_bps = [info["plasmid"] for info in read_distribution.values()]
    genome_bps = [info["genome"] for info in read_distribution.values()]
    
    # Set up the plot
    plt.figure(figsize=(12, 6))
    
    # Plot KDE for both plasmid and genome base pairs
    sns.kdeplot(x=plasmid_bps, y=genome_bps, cmap="coolwarm", fill=True, thresh=0.01)
    
    plt.minorticks_on()

    # Labeling the plot
    plt.xlabel("Base pairs from plasmid")
    plt.ylabel("Base pairs from genome")
    plt.title("Read Base Pair Distribution between Plasmid and Genome (KDE)")
    plt.savefig(f"{output}.pdf")
    plt.scatter(plasmid_bps, genome_bps, color='black', alpha=0.5, s=1, label="Reads")
    plt.savefig(f"{output}_scatter.pdf")
    print(f'Figure saved as {output}')
    # Show the plot
    plt.show()

#Get reference lengths from bam file
def get_reference_lengths(bam_file_path):
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    reference_lengths = []
    for reference in bam_file.references:
        length = bam_file.get_reference_length(reference)
        reference_lengths.append(length)

    bam_file.close()
    return reference_lengths


def identify_hotspots(read_distribution, genome_reads, genome_length, window_size=1000):
    coverage = [0] * genome_length

    for read_name, alignments in genome_reads.items():
        for genome_name, start, end in alignments:
            for i in range(start, end):
                if i < genome_length:
                    coverage[i] += 1

    hotspots = []
    for i in range(0, genome_length, window_size):
        window_coverage = sum(coverage[i:i+window_size])
        hotspots.append((i, i+window_size, window_coverage))
    
    hotspots = sorted(hotspots, key=lambda x: x[2], reverse=True)
    return hotspots

def plot_hotspots(hotspots, output, top_n=10):
    hotspots = hotspots[:top_n]
    regions = [f"{start}-{end}" for start, end, _ in hotspots]
    values = [coverage for _, _, coverage in hotspots]

    plt.figure(figsize=(12, 6))
    plt.bar(regions, values)
    plt.xlabel("Genomic Region")
    plt.ylabel("Read Coverage")
    plt.title("Top Hotspots in Genome")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output}.pdf")
    print(f'Figure saved as {output}')
    plt.show()


#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
plasmid_bam = "sorted_mapped_to_plasmid.bam"
genome_bam = "sorted_aligned_to_other_genome.bam"

# OUTPUTS
output_kde = "kde_plot_plasmid_and_genome"
output_hotspot = "genome_hotspots"
#----------------------------------------------------


plasmid_length = get_reference_lengths(plasmid_bam)[0]  
genome_length = get_reference_lengths(genome_bam)[0]

plasmid_reads = read_alignment(plasmid_bam, "plasmid")
genome_reads = read_alignment(genome_bam, "other_genome")

read_distribution = analyze_basepair_distribution(plasmid_reads, genome_reads)
plot_read_distribution(read_distribution, output_kde)

hotspots = identify_hotspots(read_distribution, genome_reads, genome_length)
plot_hotspots(hotspots, output_hotspot)
