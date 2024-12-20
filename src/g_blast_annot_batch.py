#-------------------------------------
#This code used for annotate insertions, It uses manually curated ISes from genome, and the genome, to understand the origin. If there
#are sequences outside from IS DB, it will report, even there is no manually curated ISes, it can be extracted with using this 
#Also, it will report the contaminations that are not in the genome, can coem from different source or problems with library prep.

#NEEDS: executable Blast installation from PATH and BioPython

#INPUT: TAB FILE in the format: "{read_id}\t{query_start}\t{query_end}\t{insertion_length}\t{ref_pos}\t{int(quality[read_id])}\t{insert_type}\t{is_reverse}\n"

#tab_file = "filtered_read_insertions_minimap2_with_refs_v4.txt"
#mapped_fasta = "mapped_reads.fasta"
#genom_fasta = "genome.fasta"
#is_fasta = "cleaned_IS.fasta"
#temp_dir = "temp"

#OUTPUT: TAB FILE in the format: "Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\tSubject Start\tSubject End\tNote\tExplained\tref_pos\tis_Reverse\n"


#Recep Can Altınbağ, 07 11 2024, v0.2
#-------------------------------------

import time
import sys
import os
import shutil
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from concurrent.futures import ProcessPoolExecutor, as_completed


def read_tab_file(file_path):
    tab_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            query_start = int(fields[1])
            query_end = int(fields[2])
            insertion_length = int(fields[3])
            ref_pos = fields[4]
            quality = int(fields[5])
            insert_type = fields[6]
            if fields[7] == "True":
                is_reverse = True
            else:
                is_reverse = False
            tab_data.append((read_id, query_start, query_end, insertion_length, ref_pos, quality, insert_type, is_reverse))
    return tab_data

def get_sequence_length(fasta_file, target_id):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == target_id:  
            return len(record.seq)  
    return None


def extract_sequence(mapped_fasta_path, read_id, query_start, query_end):
    for record in SeqIO.parse(mapped_fasta_path, "fasta"):
        if record.id == read_id:
            return record.seq[query_start:query_end]
    return None

def save_to_temp_fasta(sequence, temp_dir, read_id, query_start, query_end):
    file_name = f"{read_id}_{query_start}_{query_end}.fasta"
    file_path = os.path.join(temp_dir, file_name)
    SeqIO.write([SeqIO.SeqRecord(sequence, id=f"{read_id}_{query_start}_{query_end}", description="")], file_path, "fasta")
    return file_path

def overlap(hsp1, hsp2):
    return hsp1.query_start <= hsp2.query_end and hsp2.query_start <= hsp1.query_end


def process_blast_results(xml_output, type_i, is_fasta, partial_threshold, partial_len):
    best_alignments = []
    alignments = []
    with open(xml_output) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            alignments_to_keep = []
            query_id = blast_record.query
            query_length = blast_record.query_length  
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    
                    if hsp.expect != 0:
                        #print('Not a good alignment')
                        continue  # Skip if E-value is zero
                    to_continue_from_partial_IS = False
                    if type_i == "is":    
                        if abs(hsp.sbjct_end - hsp.sbjct_start + 1) / alignment.length * 100 < partial_threshold and alignment.length > partial_len:
                            to_continue_from_partial_IS = True
                    
                    if to_continue_from_partial_IS:
                        continue # Skip if there is no enough coverage for IS DB elements
                          
                    #print(hsp)
                    #print(hsp.align_length, hsp.sbjct_end,  hsp.sbjct_start)
                    
                    #coverage_percentage = abs(hsp.sbjct_end - hsp.sbjct_start + 1) / hsp.align_length * 100
                    #if coverage_percentage < 70: #Eliminate lower than 80 coverage
                    #    continue

                    keep = True
                    for kept_alignment in alignments_to_keep:
                        if overlap(hsp, kept_alignment['HSP']):
                            # If there's overlap, keep the better one
                            if hsp.score > kept_alignment['HSP'].score:
                                kept_alignment.update({
                                    'Subject ID': alignment.title.split(" ")[-1],
                                    'Identity (%)': hsp.identities / hsp.align_length * 100,
                                    'Score': hsp.score,
                                    'E-value': hsp.expect,
                                    'Query Start': hsp.query_start,
                                    'Query End': hsp.query_end,
                                    'Subject Start': hsp.sbjct_start,
                                    'Subject End': hsp.sbjct_end,
                                    'Query Length':query_length,
                                    'HSP': hsp
                                })
                            keep = False
                            break
                    if keep:
                        alignments_to_keep.append({
                            'Query ID': query_id,
                            'Subject ID': alignment.title.split(" ")[-1],
                            'Identity (%)': hsp.identities / hsp.align_length * 100,
                            'Score': hsp.score,
                            'E-value': hsp.expect,
                            'Query Start': hsp.query_start,
                            'Query End': hsp.query_end,
                            'Subject Start': hsp.sbjct_start,
                            'Subject End': hsp.sbjct_end,
                            'Query Length':query_length,
                            'HSP': hsp  # Store the HSP to compare later for overlaps
                        })
    
            for alignment in alignments_to_keep:
                alignment.pop('HSP')
                best_alignments.append(alignment)
    #If there is no alignment, add 'no blast hit'
    if not best_alignments:
        best_alignments.append({
            'Query ID': query_id,
            'Subject ID': 'no blast hit',
            'Identity (%)': 'N/A',
            'Score': 'N/A',
            'E-value': 'N/A',
            'Query Start': 0,
            'Query End': 0,
            'Subject Start': 0,
            'Subject End': 0,
            'Query Length':query_length
        })
    return best_alignments


#For writing to the file
def align_to_csv(out_str_list, alignment, note, coverage, ref_pos, is_reverse):
    
    out_str_list.append(f"{alignment['Query ID']}\t{alignment['Subject ID']}\t{alignment['Identity (%)']}\t"
                   f"{alignment['Score']}\t{alignment['E-value']}\t{alignment['Query Start']}\t"
                   f"{alignment['Query End']}\t{alignment['Subject Start']}\t{alignment['Subject End']}\t{note}\t{coverage}\t{ref_pos}\t{is_reverse}\n")

    return out_str_list




def perform_blast(query_fasta, db_fasta, db_name, xml_output, tab_output):
    os.system(f"blastn -query {query_fasta} -db {db_name} -out {xml_output} -outfmt 5")
    #os.system(f"blastn -query {query_fasta} -db {db_name} -out {tab_output} -outfmt 6") uncommented, just for debugging :)

# Function to handle each read_id’s processing, including perform_blast
def process_read(read_id, query_start, query_end, q_len, ref_pos, qual, i_type, is_reverse, tab_data, genom_fasta, temp_dir, mapped_fasta, is_fasta, is_db, genom_db):
    sequence = extract_sequence(mapped_fasta, read_id, query_start, query_end)
    if sequence:
        #print('.', end='', flush=True)
        fasta_path = save_to_temp_fasta(sequence, temp_dir, read_id, query_start, query_end)
        genom_blast_xml = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_genom_blast.xml")
        genom_blast_tab = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_genom_blast.tab")
        perform_blast(fasta_path, genom_fasta, genom_db, genom_blast_xml, genom_blast_tab)
        is_blast_xml = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_is_blast.xml")
        is_blast_tab = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_is_blast.tab")
        perform_blast(fasta_path, is_fasta, is_db, is_blast_xml, is_blast_tab)


def simple_loading_bar(current, total):
    progress = int((current / total) * 100)
    bar = '=' * (progress // 2) + ' ' * (50 - (progress // 2))
    print(f"\rProcessing: [{bar}] {progress}% ({current}/{total})", end='', flush=True)


def main_annot_single(tab_file, mapped_fasta, genom_fasta, is_fasta, temp_dir, output_csv, threshold=70, max_workers = 2, debug=False, temp_value=False, partial_threshold = 80, partial_len = 8000):
    print(f'Settings: Threshold: {threshold}, Threads: {max_workers}, Debugging: {debug}, Partial Threshold: {partial_threshold}, Partial Len: {partial_len}')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    tab_data = read_tab_file(tab_file)
    out_str_list = []
    len_of_tab = len(tab_data)
    print('Blast Process of ', len_of_tab)
    genom_dir = os.path.dirname(genom_fasta)
    genom_db = os.path.join(genom_dir, "genom_db")
    is_dir = os.path.dirname(is_fasta)
    is_db = os.path.join(is_dir, "is_db")

    # Veritabanının olup olmadığını kontrol et, yoksa oluştur
    db_files = [f"{genom_db}.{ext}" for ext in ['nin', 'nhr', 'nsq']]
    if not all(os.path.isfile(file) for file in db_files):
        os.system(f"makeblastdb -in {genom_fasta} -dbtype nucl -out {genom_db}")
    db_files = [f"{is_db}.{ext}" for ext in ['nin', 'nhr', 'nsq']]
    if not all(os.path.isfile(file) for file in db_files):
        os.system(f"makeblastdb -in {is_fasta} -dbtype nucl -out {is_db}")
    
    #print_progress_bar(index, len_of_tab, prefix='Processing', length=40)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_read, read_id, query_start, query_end, q_len, ref_pos, qual, i_type, is_reverse, tab_data, genom_fasta, temp_dir, mapped_fasta, is_fasta, is_db, genom_db)
            for read_id, query_start, query_end, q_len, ref_pos, qual, i_type, is_reverse, *_ in tab_data
        ]

        # Wait for all tasks to complete
        for i, future in enumerate(futures):
            future.result()  # This will raise any exception that occurred
            simple_loading_bar(i + 1, len(futures)) 


    print('\nBlast Processes were Ended')
    print('Time to process results')

    for i, (read_id, query_start, query_end, q_len, ref_pos, qual, i_type, is_reverse, *_) in enumerate(tab_data):

        if not debug:
            simple_loading_bar(i + 1, len_of_tab)  # Yükleme çubuğunu güncelle

        #fasta_path = save_to_temp_fasta(sequence, temp_dir, read_id, query_start, query_end)
            
        # Blasting and temp files etc.
        genom_blast_xml = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_genom_blast.xml")
        genom_blast_tab = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_genom_blast.tab")
        #perform_blast(fasta_path, genom_fasta, "genom_db", genom_blast_xml, genom_blast_tab)
                
        is_blast_xml = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_is_blast.xml")
        is_blast_tab = os.path.join(temp_dir, f"{read_id}_{query_start}_{query_end}_is_blast.tab")
        #perform_blast(fasta_path, is_fasta, "is_db", is_blast_xml, is_blast_tab)
                
        # Proces blast outputs
        genom_alignments = process_blast_results(genom_blast_xml, 'genome', None, partial_threshold, partial_len)
        is_alignments = process_blast_results(is_blast_xml, 'is', is_db, partial_threshold, partial_len)
                
        genom_dict = {alignment['Query ID']: alignment for alignment in genom_alignments} # dict for fast access

        # Main Logic: Look to the IS DB, if not enough look genome, if not enough classify as contamination!
        for alignment in is_alignments:

            if (abs(alignment['Query Start'] - alignment['Query End'])/ alignment['Query Length'] )*100 < threshold or alignment['Score'] == 'N/A':
                genom_start = genom_dict[alignment['Query ID']]['Query Start']
                genom_end = genom_dict[alignment['Query ID']]['Query End']
                exp_coverage = abs(genom_start - genom_end)/alignment['Query Length'] * 100
                if genom_dict[alignment['Query ID']]['Score'] == 'N/A' or genom_dict[alignment['Query ID']] == 'no blast hit':
                    if debug:
                        print(read_id, 'No genome alignment, probably a contamination or wrong barcode or sequencing artifact!, Read Quality: ', qual)
                    out_str_list = align_to_csv(out_str_list, genom_dict[alignment['Query ID']], 'Contamination', 0, ref_pos, is_reverse)
                else:
                    if debug:
                        print(read_id, 'Explained with Genome: ', exp_coverage )
                    out_str_list = align_to_csv(out_str_list, genom_dict[alignment['Query ID']], 'Genome', exp_coverage, ref_pos, is_reverse)
                    #print('In genome: ', genom_dict[alignment['Query ID']]['Subject Start'], genom_dict[alignment['Query ID']]['Subject End'])
                    #print('No align on IS', genom_dict[alignment['Query ID']])
            else:
                out_str_list = align_to_csv(out_str_list, alignment, 'IS_DB', (abs(alignment['Query Start'] - alignment['Query End'])/ alignment['Query Length'] )*100, ref_pos, is_reverse)
                if debug:
                    print(read_id, 'On manually curated IS DB')
            #if "3a95b9e1" in alignment['Query ID']:
            #    print(alignment['Query Start'], alignment['Query End'], alignment['Query Length'])
            #    (abs(alignment['Query Start'] - alignment['Query End'])/ alignment['Query Length'] )*100
            #    input()
   
    with open(output_csv, 'w') as out_file:
        out_file.write("Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\tSubject Start\tSubject End\tNote\tExplained\tref_pos\tis_Reverse\n")    
        out_file.write(''.join(out_str_list))
    
    if not temp_value:
        # Delet temp
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    print(f"\nOutputs {output_csv} were recorded.")

#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
'''
tab_file = "filtered_read_insertions_minimap2_with_refs_v4.txt"
mapped_fasta = "100ctab.fasta"
genom_fasta = "BIOMIG1BAC.fasta"
is_fasta = "BIOMIG_cleaned_IS.fasta"
temp_dir = "temp"
threads = 16
threshold = 70
debug = False
temp_value = True
partial_threshold = 80
partial_len = 8000
# OUTPUTS
output_csv = "best_alignments.tab"
'''
#----------------------------------------------------

#main_annot(tab_file, mapped_fasta, genom_fasta, is_fasta, temp_dir, output_csv, threshold, threads, debug, temp_value, partial_threshold, partial_len)

def find_file_with_prefix(directory, prefix):
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
        if filename == prefix:
            # Return the full path of the matching file
            return os.path.join(directory, filename)
    # If no file is found, return None
    return None

def main_annot_batch(data_dir, genom_fasta, is_fasta, temp_dir, threshold, threads, debug, temp_value, partial_threshold, partial_len):

    tab_file_name = "insertion_from_bam.tab"
    output_csv = "best_alignment.tab"

    original_stdout = sys.stdout #For debugging 
    current_dir = data_dir
    folders = [f for f in os.listdir(current_dir) if os.path.isdir(os.path.join(current_dir, f))]

    print('\n')
    print(folders)
    for folder in folders:
        if folder == temp_dir:
            continue
        print('\nProcessing: ', folder)
        output_file = os.path.join(current_dir, folder, output_csv)
        mapped_fasta = os.path.join(current_dir, folder, f"{folder}.fasta")
        tab_file = find_file_with_prefix(os.path.join(current_dir, folder), tab_file_name)
        if tab_file == None:
            print("No tab insertions file in ", folder, ", skipping!")
            continue
        print('Corresponding tab file: ', tab_file)
        if not os.path.exists(mapped_fasta):
            print("No mapped fasta file in ", folder, ", skipping!")
            continue
        
        main_annot_single(tab_file, mapped_fasta, genom_fasta, is_fasta, temp_dir, output_file, threshold, threads, debug, temp_value, partial_threshold, partial_len)
