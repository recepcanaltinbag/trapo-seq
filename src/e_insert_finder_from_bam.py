#-------------------------------------
#This code used for minimap map-ont based on Nanopore long-read data to find insertions inside the reads with using CIGAR strings from bam file (I=1)
#Also, This code can deal with fragmented alignments if there is a gap in start, end or middle of the total alingment of the specific reads to the reference
#For reverse aligned reads, their CIGAR tuples are reversed to calculate query positions 

#It is also tested for bwa-mem. But, for bwa-mem outputs I saw that insertions are not soft clipped but splitting the alignments to pieces. But, thats not a problem.

#INPUT: BAM FILE, and threshold for filtering insertions
#OUTPUT: TAB FILE in the format: "{read_id}\t{query_start}\t{query_end}\t{insertion_length}\t{ref_pos}\t{int(quality[read_id])}\t{insert_type}\t{is_reverse}\n"

#Recep Can Altınbağ, 22 10 2024, v0.2
#-------------------------------------

import pysam
from collections import defaultdict
import sys
import os

#Average qualities of reads based on read_id, return a dict {read_id:quality}
def calculate_average_quality(bamfile):
    average_quality = {}
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            
            read_id = read.query_name
            qualities = read.query_qualities
            avg_quality = sum(qualities) / len(qualities) if qualities else 0.0
            average_quality[read_id] = avg_quality
    return average_quality




# for the array as 1,1,0,0,0,0,1,1,1
# this function returns index of start of 0, and length, like (2,5,4)
# so the insertions can be found 
def find_zero_sequences(arr, insertion_threshold, ref):
    zero_sequences = []
    start = None  
    length = 0    

    for i, num in enumerate(arr):
        if num == 0:
            if start is None:
                start = i  
            length += 1    
        else:
            if start is not None:
                print('Start',start)
                if length > insertion_threshold:
                    zero_sequences.append((start, i - 1, length, ref))
                start = None
                length = 0

    if start is not None and length > insertion_threshold:
        zero_sequences.append((start, len(arr) - 1, length, ref))

    return zero_sequences

#Finding insertions in CIGAR, in some other alignment programs (bwa-mem) can assign big insertions as soft-clips!
#rev_type for if reads are reversed the ref and query must be moved differently!
def find_large_insertions(cigar_tuples, query_start, ref_start, insertion_threshold, rev_type):
    insertion_positions = []
    query_pos = query_start 
    ref_pos = ref_start   
 
    for cigar_type, length in cigar_tuples:
        if cigar_type == 0:  # Match (alignment)
            ref_pos = ref_pos + rev_type * length
            query_pos += length
        elif cigar_type == 1:  # Insertion 
            if length > insertion_threshold:  # If insertion is higher than threshold
                insertion_positions.append((query_pos, length, ref_pos))
        elif cigar_type == 2:  # Deletion 
            ref_pos = ref_pos + rev_type * length
        elif cigar_type == 4 or cigar_type == 5:
            query_pos += length

    return insertion_positions



# To deal with inserts in CIGAR strings
def get_middle_inserts(bamfile, insertion_threshold):
    reads = defaultdict(list)
    read_len_dict = {}
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                read_id = read.query_name
                ref_start = read.reference_start
                ref_end = read.reference_end

                if read.is_reverse:
                    cigar_tuples = list(reversed(read.cigartuples)) #for reverse reads, reference or query must be reversed!
                    query_start = 0
                    inserts = find_large_insertions(cigar_tuples, query_start, ref_end, insertion_threshold, -1)
                else:
                    cigar_tuples = read.cigartuples
                    query_start = 0
                    inserts = find_large_insertions(cigar_tuples, query_start, ref_start, insertion_threshold, 1)              

                if inserts != []:
                    insert_list = []
                    for insert in inserts:
                        insert_list.append((insert[0], insert[0] + insert[1], insert[1], insert[2], 'IN', read.is_reverse)) #insert 2 is ref position
                    reads[read_id] = insert_list
    return reads


# To deal with fragmented insertions
def get_read_info(bamfile, insertion_threshold):
    reads = defaultdict(list)
    read_len_dict = {}
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if not read.is_unmapped: 
                read_id = read.query_name
                #print(read_id)
                
                ref_start = read.reference_start
                ref_end = read.reference_end
                
                #IF the read is reversed, than to calculate actual query pos, we need to reverse the cigar tuples
                if read.is_reverse:
                    cigar_tuples = list(reversed(read.cigartuples))
                else:
                    cigar_tuples = read.cigartuples
                
                query_start = 0
                query_end = sum(length for op, length in cigar_tuples if op in {0, 1, 7, 8})  # Matches with insertions (1:I)
                
                
   
                #print(read.cigartuples[0],read.cigartuples[-1][1])


                if cigar_tuples[0][0] == 4 or cigar_tuples[0][0] == 5:  #If there is soft or hard clips in start or end.
                     query_start = cigar_tuples[0][1]
                
                query_end += query_start
                
                #if "037606da" in read_id:
                #    print(ref_start, ref_end)
                #    print(query_start, query_end)
                #    print(read.is_reverse)
                #    input()

                reads[read_id].append({
                    "reference_start": ref_start,
                    "reference_end": ref_end,
                    "query_start": query_start,
                    "query_end": query_end,
                    "aligned_length": sum(length for op, length in cigar_tuples if op in {0,2,7,8}),  # Matches with deletions (2:D)
                    "reverse": read.is_reverse
                })
                read_len_dict[read_id] = read.infer_read_length()


    return reads, read_len_dict


#function to write output file.
def filter_and_write_read_ids(read_data, output_file, threshold, quality):
    with open(output_file, 'w') as f:
        for read_id, sequences in read_data.items():
            for seq in sequences:
                start, end, length, ref_pos, insert_type, is_reverse = seq
                if length > threshold:
                    f.write(f"{read_id}\t{start}\t{end}\t{length}\t{ref_pos}\t{int(quality[read_id])}\t{insert_type}\t{is_reverse}\n")


#In splitted alignments, it is not easy to understand reference position, so this function gives the closest point!
def find_nearest_reference(align_table, start, end):
    nearest_reference = None
    min_diff = float('inf')
    
    for ref_start, ref_end, query_start, query_end, reverse in align_table:
       
        diff_start = abs(query_start - start)
        diff_end = abs(query_end - end)

        diff_e_start = abs(query_end - start)        
        diff_s_end = abs(query_start - end) 

        if diff_start < min_diff:
            min_diff = diff_start
            if reverse:
                nearest_reference = ref_end
            else:
                nearest_reference = ref_start

        if diff_s_end < min_diff:
            min_diff = diff_s_end
            if reverse:
                nearest_reference = ref_end
            else:
                nearest_reference = ref_start 
      

        if diff_end < min_diff:
            min_diff = diff_end
            if reverse:
                nearest_reference = ref_start
            else:
                nearest_reference = ref_end
        
        if diff_e_start < min_diff:
            min_diff = diff_e_start
            if reverse:
                nearest_reference = ref_start
            else:
                nearest_reference = ref_end 
    print('nearest: ',nearest_reference)
    return nearest_reference


def simple_loading_bar(current, total):
    progress = int((current / total) * 100)
    bar = '=' * (progress // 2) + ' ' * (50 - (progress // 2))
    print(f"\rProcessing: [{bar}] {progress}% ({current}/{total})", end='', flush=True)


def main_insert_finder_from_bam(bamfile, output_file, insertion_threshold = 500):

    #----------------------------------------------------
    # EXAMPLE USAGE -------------------------------------
    # INPUTS
    #bamfile = "sorted.bam"
    #insertion_threshold = 500
    # OUTPUTS
    #output_file = "filtered_read_insertions_minimap2_with_refs_v3.txt"
    #----------------------------------------------------

    original_stdout = sys.stdout #For debugging 
    

    read_info, read_len_dict = get_read_info(bamfile, insertion_threshold)
    reads_insertions = defaultdict(list)

    len_reads_items = read_info.items()

    for i, (read_id, alignments) in enumerate(read_info.items()):

        zero_list = [0] * read_len_dict[read_id]

        simple_loading_bar(i + 1, len(len_reads_items)) 
        sys.stdout = open(os.devnull, 'w')
        print(f"Read ID: {read_id}")
        print(len(zero_list))
        
        align_table = []

        for alignment in alignments:
            zero_list[alignment['query_start']:alignment['query_end']] = [1] * (alignment['query_end'] - alignment['query_start'])
            #print(zero_list)
            print(f"  Reference Start: {alignment['reference_start']}")
            print(f"  Reference End: {alignment['reference_end']}")
            print(f"  Query Start: {alignment['query_start']}")
            print(f"  Query End: {alignment['query_end']}")
            print(f"  Aligned Length: {alignment['aligned_length']}")
            print(f"  Aligned Reverse: {alignment['reverse']}")
            align_table.append((alignment['reference_start'], alignment['reference_end'], alignment['query_start'], alignment['query_end'], alignment['reverse']))
        
        temp_ins = find_zero_sequences(zero_list, insertion_threshold, -1)
        new_ins = []
        if temp_ins != []:
            print(f"Read ID: {read_id}")
            for each_temp_ins in temp_ins:
                #print(f"Read ID: {read_id}", alignment['reverse'])
                ref = find_nearest_reference(align_table, each_temp_ins[0], each_temp_ins[1])
                new_ins.append((each_temp_ins[0],each_temp_ins[1]+1,each_temp_ins[2],ref,'SC', alignment['reverse'])) #+1 for 0-based
            reads_insertions[read_id] = new_ins
            print(reads_insertions[read_id], alignment['reverse'])
            
            
            #if '037606da' in read_id:
            #    input()
        sys.stdout = original_stdout


    #Extend if there are same read IDs there will be no problem
    for key, value in get_middle_inserts(bamfile, insertion_threshold).items():
        reads_insertions[key].extend(value) 

    #reads_insertions.update(get_middle_inserts(bamfile, insertion_threshold))
    filter_and_write_read_ids(reads_insertions, output_file, insertion_threshold, calculate_average_quality(bamfile))

    print('\nInserts are found and writed to the folder ', output_file)