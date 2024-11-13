#-------------------------------------
#This code used for visualization of insertions and deletions based on reference positions
#Also, This code can deal with fragmented alignments, insertions can be taken from tabular file.

#INPUT: BAM FILE, and  TAB FILE of insertions in the format: "{read_id}\t{query_start}\t{query_end}\t{insertion_length}\t{ref_pos}\t{int(quality[read_id])}\t{insert_type}\n"
#OUTPUT: Graphical outputs

#Recep Can Altınbağ, 24 10 2024, v0.0
#-------------------------------------

import matplotlib.pyplot as plt
from collections import defaultdict
import pysam
import os
import pandas as pd

def parse_insertion_tab_file(file_path):
    idct = defaultdict(list)
    insertion_counts = defaultdict(int)
    weighted_insertion_counts = defaultdict(int)

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            read_id = parts[0]
            query_start = int(parts[1])
            query_end = int(parts[2])
            insertion_length = int(parts[3])
            ref_pos = int(parts[4])
            quality = int(parts[5])
            insert_type = parts[6]
            if insert_type == "IN":        #Insert type IN already can be taken from bam files, so do not need to
                continue
            insertion_counts[ref_pos] += 1
            weighted_insertion_counts[ref_pos] += 1 * insertion_length
    
    return insertion_counts, weighted_insertion_counts



# Deletion extraction from bam file
def extract_deletions_from_bam(bam_file_path, deletion_threshold):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    
    deletion_counts = defaultdict(int)
    weighted_deletion_counts = defaultdict(int)
    
    for read in bamfile.fetch():
        read_id = read.query_name
        cigar =  read.cigartuples
        ref_pos = read.reference_start
        
        for operation, length in cigar:
            if operation == 2:  # Deletion (CIGAR 'D')
                if length >= deletion_threshold:  # Threshold 
                    deletion_counts[ref_pos] += 1
                    for leng in range (0,length):
                        weighted_deletion_counts[ref_pos + leng] += 1 * length  # for weighted
            ref_pos += length if operation in {0, 2} else 0  # Deletion or match, update reference position
         
    return deletion_counts, weighted_deletion_counts


# Insertion extraction from bam file
def extract_insertions_from_bam(bam_file_path, insertion_threshold):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    
    insertion_counts = defaultdict(int)
    weighted_insertion_counts = defaultdict(int)

    for read in bamfile.fetch():
        read_id = read.query_name
        cigar =  read.cigartuples
        ref_pos = read.reference_start
        
        for operation, length in cigar:
            if operation == 1:  # Insertion (CIGAR 'I')
                if length >= insertion_threshold:  # Threshold
                    insertion_counts[ref_pos] += 1
                    weighted_insertion_counts[ref_pos] += 1 * length  # for weighted
            ref_pos += length if operation in {0, 2} else 0  # Deletion or match, update reference position
            
    return insertion_counts, weighted_insertion_counts


#Plotting
def plot_insertion_deletion_graphs(insertion_counts, deletion_counts, name, output):
    print(f"{output}_{name} is processing")
    sorted_insertion_positions = sorted(insertion_counts.keys())
    sorted_insertion_counts = [insertion_counts[pos] for pos in sorted_insertion_positions]

    sorted_deletion_positions = sorted(deletion_counts.keys())
    sorted_deletion_counts = [deletion_counts[pos] for pos in sorted_deletion_positions]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.set_xlabel('Reference Position')
    ax1.set_ylabel(f"Deletions [{name}]")

    ax1.plot(sorted_deletion_positions, sorted_deletion_counts, label=f"Deletions [{name}]", color='m', marker='x')
    ax1.legend(loc='upper left')
    ax1.grid(True)
    ax1.minorticks_on()
    
    ax2 = ax1.twinx()
    ax2.set_ylabel(f"Insertions [{name}]")
    ax2.plot(sorted_insertion_positions, sorted_insertion_counts, label=f"Insertions [{name}]", color='g', marker='s')
    ax2.legend(loc='upper right')
    ax2.minorticks_on()
    
    plt.title('Insertions and Deletions')
    plt.savefig(f"{output}_{name}.pdf")
    print(f"{output}_{name}.pdf saved")
    #plt.show()

    # Combine insertion and deletion data into a DataFrame
    data_insertions = {
        "Insertion Position": sorted_insertion_positions,
        "Insertion Count": sorted_insertion_counts,
    }
    data_deletions = {
        "Deletion Position": sorted_deletion_positions,
        "Deletion Count": sorted_deletion_counts
    }
    df_ins = pd.DataFrame(data_insertions)
    df_del = pd.DataFrame(data_deletions)

    # Save to a tab-delimited file
    df_ins.to_csv(f'{output}_{name}_insertion_table.tsv', sep='\t', index=False)
    df_del.to_csv(f'{output}_{name}_deletion_table.tsv', sep='\t', index=False)
    print(f"{output}_{name}.insertion_table.tsv saved")
    print(f"{output}_{name}.deletion_table.tsv saved")

#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
'''
tab_file_path = 'data/5bac/insertions_from_bam.tab'
bam_file_path = 'data/5bac/03_sorted_mapped_to_plasmid.bam'
output_name = "in_del_plot"
deletion_threshold = 5
insertion_threshold = 2
'''


def main_in_del_plot(tab_file_path, bam_file_path, output_name, insertion_threshold=2, deletion_threshold=5):

    dir_path = os.path.dirname(os.path.abspath(bam_file_path))
    output = os.path.join(dir_path, output_name)


    # Insertions
    insertion_counts, weighted_insertion_counts = parse_insertion_tab_file(tab_file_path) #SC types are taken from tabular file
    insertion_counts_in, weighted_insertion_counts_in = extract_insertions_from_bam(bam_file_path, insertion_threshold) #IN types are taken from bam file

    insertion_counts.update(insertion_counts_in)
    weighted_insertion_counts.update(weighted_insertion_counts_in) #Combine bam file inputs and tabular file inputs

    # Deletions
    deletion_counts, weighted_deletion_counts = extract_deletions_from_bam(bam_file_path, deletion_threshold)

    # Plotting
    plot_insertion_deletion_graphs(insertion_counts, deletion_counts, 'non-weighted', output)
    plot_insertion_deletion_graphs(weighted_insertion_counts, weighted_deletion_counts, 'weighted', output)
