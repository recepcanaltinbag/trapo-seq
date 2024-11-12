#-------------------------------------
#This code used for sequence logo analysis of Direct Repeats caused by transposition

#INPUT: Insertion (.csv file)
#OUTPUT: Graphs: Sequence Logos and Insertion Coordinates
#Recep Can Altınbağ, 04 11 2024, v0.0

#Libs to Load: pandas, matplotlib, logomaker, BioPython, statistics
#MAFFT needed for MSA
#-------------------------------------
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from logomaker import Logo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import tempfile
from Bio import AlignIO, SeqIO
import subprocess
import logomaker
import matplotlib.pyplot as plt
import statistics
from collections import Counter


def get_sequence_length(fasta_file):
    """
    Reads a single-sequence FASTA file and returns the sequence length.
    
    Parameters:
    fasta_file (str): Path to the FASTA file.
    
    Returns:
    int: The length of the sequence.
    """
    # Parse the FASTA file and get the first (and only) sequence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return len(record.seq)


# Sometimes the naming can be problematic for file cration, so it is good to clean these type of characters
def clean_filename(name):
    return name.replace('::', '_').replace('/', '_').replace('\\', '_')

#MAFFT subprocess
def run_mafft_silently(input_file, output_file):
    print('MAFFT')
    # MAFFT komutunu sessiz bir şekilde çalıştır
    try:
        # subprocess.run ile komutu çalıştır
        result = subprocess.run(
            ["mafft", "--auto", input_file],
            stdout=open(output_file, 'w'),  # Çıktıyı dosyaya yönlendir
            stderr=subprocess.PIPE,  # Hataları yakala
            check=True,  # Hata durumunda istisna fırlat
            text=True  # Çıktıyı metin olarak al
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"Bir hata oluştu: {e.stderr}")  # Hata mesajını yazdır
        return None

# Creating temp files and alignment
def align_sequences_temp(sequences, gap_threshold):
    if len(sequences) < 2:
        print('one sequence, warning')
        return []
    seq_records = []
    for index, sequence in enumerate(sequences):
        seq_id = f"seq_{index + 1}"  # ID'leri 'seq_1', 'seq_2' şeklinde oluştur
        seq_record = SeqRecord(Seq(sequence), id=seq_id)  # Seq objesi oluştur
        seq_records.append(seq_record)  # Listeye ekle

    with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta") as temp_input, \
         tempfile.NamedTemporaryFile(mode='r', suffix=".aln") as temp_output:
        
        SeqIO.write(seq_records, temp_input.name, "fasta")
        
        run_mafft_silently(temp_input.name, temp_output.name)

        print('End of the MAFFT process')
        
        alignment = AlignIO.read(temp_output.name, "fasta")
    
    alignment = remove_high_gap_columns(alignment, gap_threshold)
    print(alignment)
    return alignment

def remove_high_gap_columns(alignment, gap_threshold):
    """
    :param alignment: seqs in list
    :param gap_threshold: Maximum acceptible gap threshold default 60
    :return: Cleaned alignment.
    """
    # Sütun sayısını belirle
    num_columns = len(alignment[0].seq)
    filtered_alignment = []

    for i in range(num_columns):
       
        gap_count = ( sum(1 for record in alignment if record.seq[i] == '-') / len(alignment) ) * 100
        if gap_count < gap_threshold:
            for j in range(len(alignment)):
                if len(filtered_alignment) <= j:
                    filtered_alignment.append('')
                filtered_alignment[j] += alignment[j].seq[i]
    
    return filtered_alignment


# Finding the length
def get_best_alignment_length(alignment):
    trimmed_sequences = []
    for record in alignment:
        trimmed_seq = record.strip('-')  # Baştaki ve sondaki boşlukları kaldır
        trimmed_sequences.append(trimmed_seq)
    
    trimmed_lengths = [len(seq) for seq in trimmed_sequences]
    most_common_length = Counter(trimmed_lengths).most_common(1)[0][0]
    return most_common_length

# 4. Logo creating
def create_sequence_logo(alignment, best_length, output_file, safe_name):

    sequences = [str(record)[:best_length].upper() for record in alignment]
    
    counts_df = logomaker.alignment_to_matrix(sequences=sequences, to_type='counts')
    
    color_scheme = {
        'A': 'green',
        'C': 'blue',
        'G': 'orange',
        'T': 'red'
    }
    
    plt.figure(figsize=(10, 3)) 
    logo = logomaker.Logo(counts_df, color_scheme=color_scheme)
    logo.style_spines(visible=False)  
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d') 
    logo.ax.set_title(f'Overlap Sequence Logo for {safe_name}')
    plt.title("Sequence Logo") 
    plt.savefig(output_file, dpi=300) 
    plt.close()




#----------------------------------------------------
# EXAMPLE USAGE -------------------------------------
# INPUTS
'''
plasmid_fasta = 'data/01_pMAT1_plasmid.fasta'
input_folder = 'data/insertions'
# OUTPUTS
output_folder = 'data/dr_logos'
gap_threshold = 60
'''
#----------------------------------------------------


def main_dr_logo(plasmid_fasta, input_folder, output_folder, gap_threshold=60):
    plasmid_len = get_sequence_length(plasmid_fasta)
    os.makedirs(output_folder, exist_ok=True)

    # All _insertions.csv files will be taken into account
    all_files = [f for f in os.listdir(input_folder) if f.endswith('_insertions.csv')]

    df_list = []
    for f in all_files:
        if f.endswith('.csv') and not f.startswith('.'):  # CSV dosyası ve gizli dosya olmamalı
            file_path = os.path.join(input_folder, f)
            try:
                df_list.append(pd.read_csv(file_path))
            except Exception as e:
                print(f"Error reading {f}: {e}")

    # All data frames will be combined
    all_data = pd.concat(df_list, ignore_index=True)
    grouped = all_data.groupby('Subject ID')

    for name, group in grouped:
        safe_name = clean_filename(name)
        print("Processing: ", safe_name)
        insertion_points = group['Insertion Point'].tolist()
        overlap_sequences = group.loc[group['Repeat Length'] > 0, 'Overlap Sequence'].tolist()


        # Histogram
        plt.figure(figsize=(10, 6))
        plt.hist(insertion_points, bins=60, histtype='stepfilled', edgecolor='orange', linewidth=1.5, color='white')

        plt.xlim(0, plasmid_len)  # X ekseni maksimum değeri
        plt.title(f'Insertion Points Histogram for {safe_name}')
        plt.xlabel('Insertion Point')
        plt.ylabel('Frequency')

        # Saving histograms
        plt.savefig(os.path.join(output_folder, f'{safe_name}_insertion_points_histogram.pdf'))
        plt.close()

        alignment = align_sequences_temp(overlap_sequences, gap_threshold)  # 2. MAFFT ile hizalama
        if len(alignment) < 2:
            print(f'len of alignment {alignment} is not enough! skipping {safe_name}')
        else:
            best_length = get_best_alignment_length(alignment)  # 3. En iyi uzunluk
            create_sequence_logo(alignment, best_length, os.path.join(output_folder, f'{safe_name}_overlap_sequence_logo.pdf'),safe_name)  # 4. Sekans logosu ve kaydetme

