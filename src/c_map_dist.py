#-------------------------------------
#This code used for stats with visualization of mapped reads to the plasmid and the genome
#It can be used to understand which part of the genome is mobile and how many reads cover the plasmid
#Restriction sites and deletions from the plasmid can be seen

#INPUT: BAM FILES of Genome and plasmid maps, fastq file
#OUTPUT: Graphical outputs and statistics txt file

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import pysam
import matplotlib.pyplot as plt
from collections import defaultdict
import os

def count_reads_in_fastq(fastq_file):
    """
    Count the number of reads in a FASTQ file.
    
    A FASTQ read consists of four lines: 
    1. The header line starting with '@'.
    2. The sequence line.
    3. The '+' separator line.
    4. The quality score line.
    
    Args:
        fastq_file (str): Path to the FASTQ file.
        
    Returns:
        int: The number of reads in the FASTQ file.
    """
    read_count = 0
    
    with open(fastq_file, 'r') as file:
        # Count the total number of lines and divide by 4
        for i, line in enumerate(file):
            if i % 4 == 0:  # Every 4th line is a new read
                read_count += 1
    
    return read_count



#Reading bam file
def read_alignment(bam_file):
    if not os.path.exists(bam_file + ".bai"):
        print(f"Index file for {bam_file} not found, creating index...")
        pysam.index(bam_file)

    alignment = pysam.AlignmentFile(bam_file, "rb")
    reads = defaultdict(list)
    for read in alignment.fetch():
        if not read.is_unmapped and read.query_length != 0:
            reads[read.query_name].append((read.reference_start, read.reference_end))
            
    alignment.close()
    return reads

#Plotting the read distribution on genome and plasmid, genome mapps are intersection with plasmids
def plot_read_distribution(plasmid_reads, genome_reads, plasmid_length, genome_length, output):
    plasmid_counts = [0] * plasmid_length
    genome_counts = [0] * genome_length

    for read_alignments in plasmid_reads.values():
        for start, end in read_alignments:
            for i in range(start, end):
                if i < plasmid_length:
                    plasmid_counts[i] += 1

    for read_alignments in genome_reads.values():
        for start, end in read_alignments:
            for i in range(start, end):
                if i < genome_length:
                    genome_counts[i] += 1

    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    ax[0].plot(range(plasmid_length), plasmid_counts, label="Plasmid")
    ax[0].set_title("Read Distribution on Plasmid")
    ax[0].set_xlabel("Plasmid Position")
    ax[0].set_ylabel("Read Count")
    ax[0].minorticks_on()
    
    ax[1].plot(range(genome_length), genome_counts, label="Other Genome", color="red")
    ax[1].set_title("Read Distribution on Genome")
    ax[1].set_xlabel("Genomic Position")
    ax[1].set_ylabel("Read Count")
    ax[1].minorticks_on()
    
    plt.tight_layout()
    plt.savefig(f"{output}.pdf")
    print(f'Figure saved as {output}')
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

#To write some statistics about mapping
def summarize_mappings(plasmid_reads, genome_reads, output_txt, fastq_file):
    total_reads = set(plasmid_reads.keys()).union(set(genome_reads.keys()))
    plasmid_only = set(plasmid_reads.keys()) - set(genome_reads.keys())
    genome_only = set(genome_reads.keys()) - set(plasmid_reads.keys())
    both = set(plasmid_reads.keys()).intersection(set(genome_reads.keys()))
    total_read_all = count_reads_in_fastq(fastq_file)
    print(f"Total reads: {total_read_all}")
    print(f"Total mapped reads: {len(total_reads)}")
    print(f"Plasmid only: {len(plasmid_only)}")
    #print(f"Other genome only: {len(genome_only)}")
    print(f"Both: {len(both)}")
    print(f"Transposition Rate:{(len(both)/(len(both)+len(plasmid_only))):.2f}\n----")

    with open(output_txt, 'w') as f:   
        f.write(f"Total reads: {total_read_all}\n")
        f.write(f"Total mapped reads: {len(total_reads)}\n")
        f.write(f"Plasmid only: {len(plasmid_only)}\n")
        f.write(f"Both: {len(both)}\n")
        f.write(f"Transposition Rate:{(len(both)/(len(both)+len(plasmid_only))):.2f}\n")
    print(f"Mapping statistics written to {output_txt}")



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
#plasmid_bam = "data/5bac/03_sorted_mapped_to_plasmid.bam"

# OUTPUTS
#output_p = "read_dist_plasmid_and_genome"
#----------------------------------------------------

def main_map_dist(plasmid_bam, output_p):

    dir_path = os.path.dirname(os.path.abspath(plasmid_bam))
    output = os.path.join(dir_path, output_p)

    genome_bam = find_file_with_prefix(dir_path, "07_", ".bam")
    print('Corresponding genome bam: ', genome_bam)

    plasmid_length = get_reference_lengths(plasmid_bam)[0]  
    genome_length = get_reference_lengths(genome_bam)[0]

    print('\n\nDistribution of reads\n')
    print(f'Len of plasmid: {plasmid_length}')
    print(f'Len of genome: {genome_length}\n----')

    plasmid_reads = read_alignment(plasmid_bam)
    genome_reads = read_alignment(genome_bam)

    plot_read_distribution(plasmid_reads, genome_reads, plasmid_length, genome_length, output)


