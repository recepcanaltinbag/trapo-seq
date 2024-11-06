#-------------------------------------
#This code used for annotate insertions, It uses manually curated ISes from genome, and the genome, to understand the origin. If there
#are sequences outside from IS DB, it will report, even there is no manually curated ISes, it can be extracted with using this 
#Also, it will report the contaminations that are not in the genome, can coem from different source or problems with library prep.

#NEEDS: executable Blast installation from PATH and BioPython

#INPUT: TAB FILE in the format: "{read_id}\t{query_start}\t{query_end}\t{insertion_length}\t{ref_pos}\t{int(quality[read_id])}\t{insert_type}\t{is_reverse}\n"

#Recep Can Altınbağ, 06 11 2024, v0.0
#-------------------------------------


import os
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

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
            if fields[7] == True:
                is_reverse = True
            else:
                is_reverse = False
            tab_data.append((read_id, query_start, query_end, insertion_length, ref_pos, quality, insert_type, is_reverse))
    return tab_data

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

def perform_blast(query_fasta, db_fasta, db_name, xml_output, tab_output):
    # Veritabanının olup olmadığını kontrol et, yoksa oluştur
    db_files = [f"{db_name}.{ext}" for ext in ['nin', 'nhr', 'nsq']]
    if not all(os.path.isfile(file) for file in db_files):
        os.system(f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_name}")

    # Hem XML hem de TAB formatında çıktı alarak BLAST çalıştır
    os.system(f"blastn -query {query_fasta} -db {db_name} -out {xml_output} -outfmt 5")
    os.system(f"blastn -query {query_fasta} -db {db_name} -out {tab_output} -outfmt 6")


def process_blast_results(xml_output):
    best_alignments = []
    alignments = []
    with open(xml_output) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            alignments_to_keep = []
            query_id = blast_record.query
            query_length = blast_record.query_length  # Query'nin total len
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect != 0:
                        print('Not a good alignment')
                        continue  # Skip if E-value is zero
                    
                    #print(hsp)
                    #print(hsp.align_length, hsp.sbjct_end,  hsp.sbjct_start)
                    
                    coverage_percentage = abs(hsp.sbjct_end - hsp.sbjct_start + 1) / hsp.align_length * 100
                    if coverage_percentage < 50: #Eliminate lower than 50 coverage
                        continue

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
def align_to_csv(out_file, alignment, note, coverage, ref_pos, is_reverse, sequence):
    
    out_file.write(f"{alignment['Query ID']}\t{alignment['Subject ID']}\t{alignment['Identity (%)']}\t"
                   f"{alignment['Score']}\t{alignment['E-value']}\t{alignment['Query Start']}\t"
                   f"{alignment['Query End']}\t{alignment['Subject Start']}\t{alignment['Subject End']}\t{note}\t{coverage}\t{ref_pos}\t{is_reverse}\t{sequence}\n")

    return


def main(tab_file, mapped_fasta, genom_fasta, is_fasta, temp_dir, output_csv, threshold=80):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    tab_data = read_tab_file(tab_file)
    
    with open(output_csv, 'w') as out_file:
        out_file.write("Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\tSubject Start\tSubject End\tNote\tExplained\tref_pos\tis_Reverse\tSequence\n")
        for read_id, query_start, query_end, q_len, ref_pos, qual, i_type, is_reverse, *_ in tab_data:
            sequence = extract_sequence(mapped_fasta, read_id, query_start, query_end)
            if sequence:
                fasta_path = save_to_temp_fasta(sequence, temp_dir, read_id, query_start, query_end)
                
                # Blasting and temp files etc.
                genom_blast_xml = os.path.join(temp_dir, f"{read_id}_genom_blast.xml")
                genom_blast_tab = os.path.join(temp_dir, f"{read_id}_genom_blast.tab")
                perform_blast(fasta_path, genom_fasta, "genom_db", genom_blast_xml, genom_blast_tab)
                
                is_blast_xml = os.path.join(temp_dir, f"{read_id}_is_blast.xml")
                is_blast_tab = os.path.join(temp_dir, f"{read_id}_is_blast.tab")
                perform_blast(fasta_path, is_fasta, "is_db", is_blast_xml, is_blast_tab)
                
                # Proces blast outputs
                genom_alignments = process_blast_results(genom_blast_xml)
                is_alignments = process_blast_results(is_blast_xml)
                

                genom_dict = {alignment['Query ID']: alignment for alignment in genom_alignments} # dict for fast access

                # Main Logic: Look IS DB, if not enough look genome, if not enough classify as contamination!
                for alignment in is_alignments:
                    
                    if (abs(alignment['Query Start'] - alignment['Query End'])/ alignment['Query Length'] )*100 < threshold or alignment['Score'] == 'N/A':
                        if genom_dict[alignment['Query ID']]['Score'] == 'N/A':
                            print('No genome alignment, probably a contamination or wrong barcode!')
                            align_to_csv(out_file, genom_dict[alignment['Query ID']], 'Contamination', 0, ref_pos, is_reverse, sequence)
                        else:
                            genom_start = genom_dict[alignment['Query ID']]['Query Start']
                            genom_end = genom_dict[alignment['Query ID']]['Query End']
                            exp_coverage = abs(genom_start - genom_end)/alignment['Query Length'] * 100
                            align_to_csv(out_file, genom_dict[alignment['Query ID']], 'Genome', exp_coverage, ref_pos, is_reverse, sequence)
                            print('Explained: ', abs(genom_start - genom_end)/alignment['Query Length'] * 100 )
                            print('In genome: ', genom_dict[alignment['Query ID']]['Subject Start'], genom_dict[alignment['Query ID']]['Subject End'])
                            print('No align on IS', genom_dict[alignment['Query ID']])
                    else:
                        align_to_csv(out_file, alignment, 'IS_DB', (abs(alignment['Query Start'] - alignment['Query End'])/ alignment['Query Length'] )*100, ref_pos, is_reverse, sequence)
                        print('On manually curated IS DB')

    print(f"Outputs {output_csv} were recorded.")

#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
tab_file = "filtered_read_insertions_minimap2_with_refs_v3.txt"
mapped_fasta = "mapped.fasta"
genom_fasta = "02_VD2Genome.fasta"
is_fasta = "vd2_is_curated.fasta"
temp_dir = "temp"
# OUTPUTS
output_csv = "best_alignments.csv"
#----------------------------------------------------

main(tab_file, mapped_fasta, genom_fasta, is_fasta, temp_dir, output_csv)
