import os
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio import SeqIO


def save_to_csv(df, file_path):
    # Dosya mevcutsa, başlık olmadan ekle
    if os.path.isfile(file_path):
        # Eğer dosya zaten varsa, başlık eklemeden veriyi ekle
        df.to_csv(file_path, mode='a', header=False, index=False)
    else:
        # Eğer dosya mevcut değilse, başlık ile veriyi yaz
        df.to_csv(file_path, mode='w', header=True, index=False)


# BLAST sonuçlarını analiz etmek için fonksiyon
def analyze_blast_results(blast_output, query_start_csv, query_end_csv, fasta_path, threshold=100, max_diff=2000):
    
    try:
        df = pd.read_csv(blast_output)
    except pd.errors.EmptyDataError:
        print(f"{blast_output} boş ya da okunabilir veri yok.")
        df = None
        return []
    
    
    blast_df = pd.read_csv(blast_output, sep='\t', header=None)
    blast_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    
    
    
    for record in SeqIO.parse(query_file, "fasta"):
        sequence = str(record.seq)
        break  # Tek bir kayıt varsa, döngüden çık

    overlaps = []
    #print(blast_df)
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
        print(query_start_csv, query_end_csv)
        print(blast_df['sstart'])

  
        blast_df['sstart_qs_diff'] = (blast_df['sstart'] - query_start_csv).abs()
        blast_df['ssend_qs_diff'] = (blast_df['send'] - query_start_csv).abs()

        blast_df['sstart_qe_diff'] = (blast_df['sstart'] - query_end_csv).abs()
        blast_df['ssend_qe_diff'] = (blast_df['send'] - query_end_csv).abs()
        
        closest_row_idx = blast_df[['sstart_qs_diff', 'ssend_qs_diff']].min(axis=1).idxmin()
        closest_rowE_idx = blast_df[['sstart_qe_diff', 'ssend_qe_diff']].min(axis=1).idxmin()


        # Farkın hangi sütunda daha düşük olduğunu kontrol et
        if blast_df.loc[closest_row_idx, 'sstart_qs_diff'] > threshold and blast_df.loc[closest_row_idx, 'ssend_qs_diff'] > threshold:
            print('ERROR: Value is bigger than threshold!')
            print(blast_df)
            
        else:
            if blast_df.loc[closest_row_idx, 'sstart_qs_diff'] < blast_df.loc[closest_row_idx, 'ssend_qs_diff']:
                pMAT1 = blast_df.loc[closest_row_idx, 'qstart']  # sstart_qs_diff küçükse qstart'ı al
            else:
                pMAT2 = blast_df.loc[closest_row_idx, 'qend']  # ssend_qs_diff küçükse qend'i al
                print('forward?')
        
        if blast_df.loc[closest_rowE_idx, 'sstart_qe_diff'] > threshold and blast_df.loc[closest_rowE_idx, 'ssend_qe_diff'] > threshold:
            print('ERROR: Value is bigger than threshold!')
            print(blast_df)
            
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

        print(f"En düşük farkın olduğu değere karşılık gelen değer: {pMAT1} : {pMAT2}")
        #print(blast_df)

        repeat_length = max(0, pMAT2 - pMAT1)
        
        # FASTA dosyasından nükleotidleri al
        if repeat_length > 0 and pMAT1 >= 0 and pMAT2 >= 0:
            overlap_sequence = sequence[pMAT1 - 1 : pMAT2 -1]
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

def extract_sequence(mapped_fasta_path, read_id, query_start, query_end, open_gap=200):
    # Ensure query_start and query_end are integers
    query_start = int(query_start)
    query_end = int(query_end)
    
    for record in SeqIO.parse(mapped_fasta_path, "fasta"):
        if record.id == read_id:
            # Adjust the start and end positions by 200 base pairs
            #adjusted_start = max(0, query_start - open_gap)  # Ensure start doesn't go below 0
            #adjusted_end = min(len(record.seq), query_end + open_gap)  # Ensure end doesn't exceed the sequence length
            
            # Adjust the sequence and return the modified record
            #record.seq = record.seq[adjusted_start:adjusted_end]
            return record
    
    return None


current_dir = "."
# Temp klasörünü oluşturma (yoksa)
temp_dir = os.path.join(current_dir, "temp")
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

open_gap = 500
query_file = "01_pMAT1_plasmid.fasta"

# Klasörleri gezip her klasördeki işlemleri gerçekleştiriyoruz
for folder in os.listdir(current_dir):
    folder_path = os.path.join(current_dir, folder)
    if os.path.isdir(folder_path):
        print(folder_path)


        if os.path.exists(f"{folder}_insertions.csv"):
            print(f"{folder} zaten var, atlanıyor...")
            continue
        

        # Best alignment CSV dosyasını ve FASTA dosyasını al
        best_alignment_csv_files = [f for f in os.listdir(folder_path) if f.endswith("_best_alignments.tab")]
        
        if not best_alignment_csv_files:
            print(f"{folder_path} klasöründe '_best_alignments.tab' dosyası bulunamadı, işlem atlanıyor.")
            continue  # Bu klasör için işlemi atla
        
        best_alignment_csv = best_alignment_csv_files[0]  # Dosya varsa ilkini al
        fasta_file = f"{folder}.fasta"
        
        best_alignment_path = os.path.join(folder_path, best_alignment_csv)
        fasta_path = os.path.join(folder_path, fasta_file)
        
        # CSV dosyasını oku
        df = pd.read_csv(best_alignment_path, sep='\t')
        
        # Her bir Query ID için işlem yap
        for index, row in df.iterrows():
            query_id = row['Query ID'].split('_')[0]
            q_st = int(row['Query ID'].split('_')[1])
            q_end = int(row['Query ID'].split('_')[2])
            subject_id = row['Subject ID']


            # Query ID'ye karşılık gelen sekansı FASTA dosyasından çek
            temp_fasta = os.path.join(temp_dir, f"{query_id}.fasta")
            
            with open(temp_fasta, "w") as temp_fasta_file:
                SeqIO.write(extract_sequence(fasta_path, query_id, q_st, q_end, open_gap), temp_fasta_file, "fasta")
            
            # Query için BLAST taraması yap
            blast_output = os.path.join(temp_dir, f"{query_id}_blast.out")
            
            #Ekim116-2024de comment yaptım gerekirse açarsın
            blastn_cline = NcbiblastnCommandline(query=query_file, subject=temp_fasta, outfmt=6, out=blast_output)
            stdout, stderr = blastn_cline()
            
            # BLAST sonuçlarını analiz et ve overlap hesapla (threshold ile)
            overlaps = analyze_blast_results(blast_output, q_st, q_end, fasta_path)
            
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
                save_to_csv(result_df, f"{folder}_insertions.csv")
                
            print('control the file')
        
