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
    if not os.path.exists(bam_file + ".bai"):
        print(f"Index file for {bam_file} not found, creating index...")
        pysam.index(bam_file)

    # Open BAM file
    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)

    for read in alignment.fetch():
        if not read.is_unmapped:
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
    if not os.path.exists(bam_file + ".bai"):
        print(f"Index file for {bam_file} not found, creating index...")
        pysam.index(bam_file)

    # Open BAM file
    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)

    for read in alignment.fetch():
        if not read.is_unmapped:
            query_name = read.query_name
            q_start, q_end = read.query_alignment_start, read.query_alignment_end
            q_len = read.query_length
            is_redundant = False
            is_reverse = False
            if read.is_reverse:
                is_reverse = True
            if q_len == 0:
                continue
            # Check for significant overlap with existing alignments for this read
            for existing_genome, existing_q_start, existing_q_end in reads[query_name]:
                
                overlap_percentage = calculate_overlap_qstart_end(q_start, q_end, existing_q_start, existing_q_end)
                #if 'dd1582bb' in query_name:
                    #print(read)
                #    if read.is_reverse:
                #        print('Reverse')
                #    print(q_len, query_name, q_start, q_end, overlap_percentage)
                #    input()
                # If overlap exceeds the threshold, skip this alignment
                if overlap_percentage > overlap_threshold:
                    is_redundant = True
                    break

            # Only add alignment if it's not redundant
            if not is_redundant:
                reads[query_name].append((genome_name, q_start, q_end))
                
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


def plot_read_distribution(read_distribution, output, plasmid_threshold, genome_threshold):
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
    #plt.show()

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
    #plt.show()


def find_file_with_prefix(directory, prefix, enfix):
    """
    Searches for a file that starts with the specified prefix in the given directory.
    
    Parameters:
    directory (str): Path to the directory to search in.
    prefix (str): The prefix to search for (default is "07").
    
    Returns:
    str: Full path of the first file found with the specified prefix, or None if not found.
    """
    # Iterate through files in the directory
    for filename in os.listdir(directory):
        # Check if the file starts with the specified prefix
        if filename.startswith(prefix) and filename.endswith(enfix):
            # Return the full path of the matching file
            return os.path.join(directory, filename)
    # If no file is found, return None
    return None

#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
'''
plasmid_bam = "data/5bac/03_sorted_mapped_to_plasmid.bam"
genome_bam = "data/5bac/07_sorted_aligned_to_other_genome.bam"

# OUTPUTS
output_kde = "kde_plot_plasmid_and_genome"
#output_hotspot = "genome_hotspots"
'''
#--
#--------------------------------------------------

def main_kde_mobile(plasmid_bam, output_kde, plasmid_threshold, genome_threshold):
    dir_path = os.path.dirname(os.path.abspath(plasmid_bam))
    output = os.path.join(dir_path, output_kde)

    genome_bam = find_file_with_prefix(dir_path, "07_", ".bam")
    print('Corresponding genome bam: ', genome_bam)
    plasmid_length = get_reference_lengths(plasmid_bam)[0]  
    genome_length = get_reference_lengths(genome_bam)[0]

    #plasmid_reads = read_alignment(plasmid_bam, "plasmid")
    plasmid_reads = read_alignment_plasmid(plasmid_bam, "plasmid")

    #genome_reads = read_alignment(genome_bam, "other_genome")
    genome_reads = read_unique_alignment_genome(genome_bam, "other_genome")

    read_distribution = analyze_basepair_distribution(plasmid_reads, genome_reads)

    write_read_distribution_to_file(read_distribution, output)
    print("output file tabular: ", output + '.tab')

    plot_read_distribution(read_distribution, output, plasmid_threshold, genome_threshold)

    #Hotspots need to be develpoed better, homologous regions in the genome can be problem.
    #hotspots = identify_hotspots(read_distribution, genome_reads, genome_length)
    #plot_hotspots(hotspots, output_hotspot)
