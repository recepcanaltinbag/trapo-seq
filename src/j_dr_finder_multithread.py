import os
import pandas as pd
from Bio import SeqIO
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor


def save_to_csv(df, file_path):
    # Dosya mevcutsa, başlık olmadan ekle
    if os.path.isfile(file_path):
        # Eğer dosya zaten varsa, başlık eklemeden veriyi ekle
        df.to_csv(file_path, mode='a', header=False, index=False)
    else:
        # Eğer dosya mevcut değilse, başlık ile veriyi yaz
        df.to_csv(file_path, mode='w', header=True, index=False)



def split_sequences(input_fasta, batch_size, output_dir):
    """
    Split the input fasta file into smaller batches.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    batch_count = len(sequences) // batch_size + (1 if len(sequences) % batch_size else 0)
    
    batch_files = []
    for i in range(batch_count):
        batch_sequences = sequences[i * batch_size: (i + 1) * batch_size]
        batch_file = os.path.join(output_dir, f"batch_{i+1}.fasta")
        
        with open(batch_file, "w") as batch_fasta:
            SeqIO.write(batch_sequences, batch_fasta, "fasta")
        
        batch_files.append(batch_file)
    
    return batch_files


def get_total_sequences(input_fasta):
    """
    Get the total number of sequences in the input FASTA file.
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    return len(sequences)

def make_blast_db(batch_fasta, db_name):
    """
    Create a BLAST database for the given batch fasta file.
    """
    subprocess.run(
        ["makeblastdb", "-in", batch_fasta, "-dbtype", "nucl", "-out", db_name],
        check=True, 
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL
    )

def run_blast(query_fasta, db_name, output_file):
    """
    Run BLAST on the given query against the given database.
    """
    subprocess.run(
        ["blastn", "-query", query_fasta, "-db", db_name, "-outfmt", "6", "-out", output_file],
        check=True, 
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.DEVNULL
    )

def process_blast_m(input_fasta, query_fasta, batch_size, blast_db_dir, output_dir):
    """
    Split the database, create BLAST DBs for each batch, run BLAST for the query, and merge the results.
    """
    # Step 1: Split the input fasta into smaller batches
    print("Splitting database into smaller batches...")
    batch_files = split_sequences(input_fasta, batch_size, blast_db_dir)
    
    # Step 2: Run BLAST for each batch
    output_files = []
    for idx, batch_file in enumerate(batch_files):
        simple_loading_bar(idx + 1, len(batch_files))
        db_name = os.path.join(blast_db_dir, f"batch_{idx+1}_db")
        output_file = os.path.join(output_dir, f"batch_{idx+1}_blast.out")
        
        # Step 2.1: Create BLAST database for this batch
        #print(f"Creating BLAST DB for {batch_file}...")
        make_blast_db(batch_file, db_name)
        
        # Step 2.2: Run BLAST for the query against the current batch database
        #print(f"Running BLAST for query against {batch_file}...")
        run_blast(query_fasta, db_name, output_file)
        
        output_files.append(output_file)
    
    # Step 3: Merge the output files into one
    merged_output = os.path.join(output_dir, "merged_blast_results.out")
    with open(merged_output, "w") as outfile:
        for file in output_files:
            with open(file, "r") as infile:
                outfile.write(infile.read() + "\n")  # Append content from each file
    
    print(f"All BLAST results merged into {merged_output}")
    print('All DBs are ready')
    return merged_output





# BLAST sonuçlarını analiz etmek için fonksiyon
def analyze_blast_results(query_id, blast_df_raw, query_start_csv, query_end_csv, fasta_path, query_file, threshold=200, repeat_th=30):



    blast_df_raw.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    

    blast_df = blast_df_raw[blast_df_raw['sseqid'] == query_id].copy()

    # Eğer filtreleme sonucunda hiçbir satır kalmazsa None veya uyarı döndür
    if blast_df.empty:
        print(f"No rows found for sseqid: {query_id}")
        return []  # veya uygun bir değer döndürebilirsiniz
    
    for record in SeqIO.parse(query_file, "fasta"):
        sequence = str(record.seq)
        break  # Tek bir kayıt varsa, döngüden çık

    overlaps = []
    #print(blast_df)
    #input()
    #print(len(blast_df))

    if len(blast_df) == 1:
        # Tek bir satır varsa, CSV'deki start ve end ile subject start ve end arasındaki yakınlığa bak
        single_row = blast_df.iloc[0]
        print('SINGLE:',single_row)
        # Query Start ve End ile Subject Start ve End arasındaki farkları hesapla
        query_start_distance = abs(single_row['sstart'] - query_start_csv)
        query_end_distance = abs(single_row['send'] - query_end_csv)
        
        if query_start_distance < query_end_distance:
            insertion_point = single_row['qstart']
        else:
            insertion_point = single_row['qend']
        
        repeat_length = 0  # Tek bir satırda overlap hesaplanamaz

        # FASTA dosyasından nükleotidleri al
        if repeat_length > 0:
            overlap_sequence = sequence[int(single_row['qstart']) - 1 : int(single_row['qend']) -1]
        else:
            overlap_sequence = "N/A"
        print(insertion_point)
        overlaps.append({
            "Query Start CSV": query_start_csv,
            "Query End CSV": query_end_csv,
            "Start Subject": single_row['sstart'],
            "End Subject": single_row['send'],
            "Insertion Point": insertion_point,
            "Repeat Length": repeat_length,
            "Overlap Sequence": overlap_sequence
        })
    else:
        # En yakın subject start ve subject end pozisyonlarını bulma
        #closest_start_row = blast_df.iloc[(blast_df['sstart'] - query_start_csv).abs().argsort()[:1]]
        #closest_end_row = blast_df.iloc[(blast_df['send'] - query_end_csv).abs().argsort()[:1]]
        print('Query start:', query_start_csv, 'Query end:', query_end_csv)
        print(blast_df['sstart'])

  
        blast_df['sstart_qs_diff'] = (blast_df['sstart'] - query_start_csv).abs()
        blast_df['ssend_qs_diff'] = (blast_df['send'] - query_start_csv).abs()

        blast_df['sstart_qe_diff'] = (blast_df['sstart'] - query_end_csv).abs()
        blast_df['ssend_qe_diff'] = (blast_df['send'] - query_end_csv).abs()
        
        closest_row_idx = blast_df[['sstart_qs_diff', 'ssend_qs_diff']].min(axis=1).idxmin()
        closest_rowE_idx = blast_df[['sstart_qe_diff', 'ssend_qe_diff']].min(axis=1).idxmin()

        print(blast_df)

        print('Closest_row_idx:',closest_row_idx)
        print('closest_rowE_idx:',closest_rowE_idx)


        # Farkın hangi sütunda daha düşük olduğunu kontrol et
        if blast_df.loc[closest_row_idx, 'sstart_qs_diff'] > threshold and blast_df.loc[closest_row_idx, 'ssend_qs_diff'] > threshold:
            print('ERROR: Value is bigger than threshold!')
            
            
        else:
            if blast_df.loc[closest_row_idx, 'sstart_qs_diff'] < blast_df.loc[closest_row_idx, 'ssend_qs_diff']:
                pMAT1 = blast_df.loc[closest_row_idx, 'qstart']  # sstart_qs_diff küçükse qstart'ı al
            else:
                pMAT2 = blast_df.loc[closest_row_idx, 'qend']  # ssend_qs_diff küçükse qend'i al
                print('forward?')
        
        if blast_df.loc[closest_rowE_idx, 'sstart_qe_diff'] > threshold and blast_df.loc[closest_rowE_idx, 'ssend_qe_diff'] > threshold:
            print('ERROR: Value is bigger than threshold!')
            
            
        else:
            if blast_df.loc[closest_rowE_idx, 'sstart_qe_diff'] < blast_df.loc[closest_rowE_idx, 'ssend_qe_diff']:
                pMAT1 = blast_df.loc[closest_rowE_idx, 'qstart']  # sstart_qs_diff küçükse qstart'ı al
            else:
                pMAT2 = blast_df.loc[closest_rowE_idx, 'qend']  # ssend_qs_diff küçükse qend'i al
                print('reverse read?')
        # Sonucu göster
        try:
            print(pMAT1)
            print('Error pMAT1')
        except UnboundLocalError:
            pMAT1 = -1
        try:
            print(pMAT2)
        except UnboundLocalError:
            print('Error pMAT2')
            pMAT2 = -1    

        print(f"The lowest value: {pMAT1} : {pMAT2}")
        #print(blast_df)

        # Careful, Blast results are 1-based. But sequence and programming array indexes are 0-based!
        # So, for length it has to +1
        # For [] calculations -> [start -1 : end] will be true!
        repeat_length = max(0, pMAT2 - pMAT1 + 1)
        if repeat_length > 0 and pMAT1 >= 0 and pMAT2 >= 0 and repeat_length < repeat_th:
            overlap_sequence = sequence[pMAT1 - 1 : pMAT2]
        else:
            overlap_sequence = "N/A"
            repeat_length = 0
        print(pMAT1, repeat_length, overlap_sequence)
        
        if pMAT1 > 0:
            insertion_point = pMAT1
        elif pMAT2 > 0:
            insertion_point = pMAT2
        else:
            insertion_point = -1

        overlaps.append({
                "Query Start CSV": query_start_csv,
                "Query End CSV": query_end_csv,
                "Start Subject": pMAT1,
                "End Subject": pMAT2,
                "Insertion Point": insertion_point,
                "Repeat Length": repeat_length,
                "Overlap Sequence": overlap_sequence})
    return overlaps

def copy_fasta(input_fasta, output_fasta):
    """
    Verilen FASTA dosyasını alır ve başka bir dosyaya yazar.
    
    :param input_fasta: Okunacak FASTA dosyasının yolu
    :param output_fasta: FASTA verilerini yazacağımız dosyanın yolu
    """
    try:
        # Girdi dosyasını açma ve çıkış dosyasına yazma
        with open(input_fasta, 'r') as input_file:
            with open(output_fasta, 'w') as output_file:
                # FASTA verilerini kopyala
                for record in SeqIO.parse(input_file, "fasta"):
                    SeqIO.write(record, output_file, "fasta")
        print(f"FASTA data has been successfully written to {output_fasta}")
    
    except FileNotFoundError:
        print(f"Error: {input_fasta} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def simple_loading_bar(current, total):
    progress = int((current / total) * 100)
    bar = '=' * (progress // 2) + ' ' * (50 - (progress // 2))
    print(f"\rProcessing: [{bar}] {progress}% ({current}/{total})", end='', flush=True)



def perform_blast_with_os(query_fasta, subject_fasta, fasta_db, xml_output):
    os.system(f"makeblastdb -in {subject_fasta} -dbtype nucl -out {fasta_db} > /dev/null 2>&1")
    os.system(f"blastn -query {query_fasta} -db {fasta_db} -out {xml_output} -outfmt 6 -max_target_seqs 0") #For unlimited hits



def process_row(index, row, blast_out, fasta_path, plasmid_fasta, out_ins_file, len_reads_items, original_stdout):
    simple_loading_bar(index + 1, len_reads_items)
    sys.stdout = open(os.devnull, 'w')
    
    query_id = row['Query ID'].split('_')[0]
    q_st = int(row['Query ID'].split('_')[1])
    q_end = int(row['Query ID'].split('_')[2])
    subject_id = row['Subject ID']
    
    overlaps = analyze_blast_results(query_id, blast_out, q_st, q_end, fasta_path, plasmid_fasta)
    
    # Writing results to CSV
    for overlap in overlaps:
        result_df = pd.DataFrame({
            "Read ID": [query_id],
            "Subject ID": [subject_id],
            "Query Start CSV": [q_st],
            "Query End CSV": [q_end],
            "Start Subject": [overlap["Start Subject"]],
            "End Subject": [overlap["End Subject"]],
            "Insertion Point": [overlap["Insertion Point"]],
            "Repeat Length": [overlap["Repeat Length"]],
            "Overlap Sequence": [overlap["Overlap Sequence"]]
        })
        save_to_csv(result_df, out_ins_file)
    
    sys.stdout = original_stdout

# Function to handle multithreading
def process_with_threads(df, blast_out, fasta_path, plasmid_fasta, out_ins_file, len_reads_items):
    # Using ThreadPoolExecutor to process rows in parallel
    with ThreadPoolExecutor() as executor:
        futures = []
        for index, row in df.iterrows():
            futures.append(executor.submit(process_row, index, row, blast_out, fasta_path, plasmid_fasta, out_ins_file, len_reads_items, sys.stdout))
        
        # Wait for all threads to complete
        for future in futures:
            future.result()


def main_dr_finder(current_dir, plasmid_fasta, open_gap, force_w=True):
    
    original_stdout = sys.stdout #For debugging 

    #current_dir = "."
    # Temp klasörünü oluşturma (yoksa)
    temp_dir = os.path.join(current_dir, "temp-dr-finder")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    ins_dir = os.path.join(current_dir, "insertions")
    if not os.path.exists(ins_dir):
        os.makedirs(ins_dir)


    #open_gap = 500
    #query_file = "01_pMAT1_plasmid.fasta"




    # Klasörleri gezip her klasördeki işlemleri gerçekleştiriyoruz
    for folder in os.listdir(current_dir):
        folder_path = os.path.join(current_dir, folder)
        

        if os.path.isdir(folder_path):
            print(folder_path)
            out_ins_file = os.path.join(ins_dir, f"{folder}_dr_insertionsv2.csv")


            if os.path.exists(out_ins_file) and force_w:
                print(f"{out_ins_file} already exist, deleting...")
                os.remove(out_ins_file)
            elif os.path.exists(out_ins_file) and not force_w:
                print(f"{out_ins_file} already exist, skipping...")
                continue
            else:
                print(f"{out_ins_file} not exist, will write...")
            
            # Best alignment CSV dosyasını ve FASTA dosyasını al
            best_alignment_csv_files = [f for f in os.listdir(folder_path) if f.endswith("best_alignment.tab")]
            
            if not best_alignment_csv_files:
                print(f"{folder_path} can not be found 'best_alignments.tab' , Skipping")
                continue  # Bu klasör için işlemi atla
            
            best_alignment_csv = best_alignment_csv_files[0]  # Dosya varsa ilkini al
            fasta_file = f"{folder}.fasta"
            fasta_path = os.path.join(folder_path, fasta_file)
            temp_fasta = os.path.join(temp_dir, fasta_file)
            #blast_out = os.path.join(temp_dir, 'blast_out')

            copy_fasta(fasta_path, temp_fasta)

            temp_fasta_db = os.path.splitext(temp_fasta)[0] + '_db'

            batch_size = 100
            
            blast_out = process_blast_m(temp_fasta, plasmid_fasta, batch_size, temp_dir, temp_dir)
            #blast_out = process_large_blast(plasmid_fasta, batch_size, temp_fasta, temp_fasta_db, temp_dir)
            #perform_blast_with_os(plasmid_fasta, temp_fasta, temp_fasta_db, blast_out)

            best_alignment_path = os.path.join(folder_path, best_alignment_csv)
            
            # CSV dosyasını oku
            df = pd.read_csv(best_alignment_path, sep='\t')
            
            len_reads_items = len(df)
            # Her bir Query ID için işlem yap
            for index, row in df.iterrows():
                simple_loading_bar(index + 1, len_reads_items)
                sys.stdout = open(os.devnull, 'w')
                
                query_id = row['Query ID'].split('_')[0]
                q_st = int(row['Query ID'].split('_')[1])
                q_end = int(row['Query ID'].split('_')[2])
                subject_id = row['Subject ID']

                try:
                    df = pd.read_csv(blast_out)
                except pd.errors.EmptyDataError:
                    print(f"{blast_out} boş ya da okunabilir veri yok.")
                    df = None
                    continue    
    
                blast_df_raw = pd.read_csv(blast_out, sep='\t', header=None)
                overlaps = analyze_blast_results(query_id, blast_df_raw, q_st, q_end, fasta_path, plasmid_fasta)
                
                # Örnek sonuçları CSV'ye yazma
                for overlap in overlaps:
                    result_df = pd.DataFrame({
                        "Read ID": [query_id],
                        "Subject ID": [subject_id],
                        "Query Start CSV": [q_st],
                        "Query End CSV": [q_end],
                        "Start Subject": [overlap["Start Subject"]],
                        "End Subject": [overlap["End Subject"]],
                        "Insertion Point": [overlap["Insertion Point"]],
                        "Repeat Length": [overlap["Repeat Length"]],
                        "Overlap Sequence": [overlap["Overlap Sequence"]]
                    })
                    save_to_csv(result_df, out_ins_file)
                sys.stdout = original_stdout
            print('Control the file')
            



main_dr_finder("data", "data/01_pMAT1_plasmid.fasta", 500)