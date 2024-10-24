

# 1 - FOLDER PREPARATION 

#After changing names of folders: barcode numbers -> sample or condition names
#This will combine unzip files and concatanate in to single *_merged.fastq

for dir in */; do
    cd "$dir" || continue
    for file in *.gz; do
        gunzip "$file"
    done
    cat * > "${dir%/}_merged.fastq"
    cd ..
done

# 1 - ALIGNMENT 

# For each folder which has merged fastq files, this script will map them to the plasmid, and genome, and creates
# just mapped fastq files and fasta file 
# Reads mapped to -> plasmid -> mapped to genome -> 

# Statistics of reads:
#  Total Reads: A
#  Plasmid Only: B
#  Genome and Plasmid: C 
#   ->Transposition Ratio: C/(C+B)


for dir in */; do
    cd "$dir" || continue

    QUERY_FASTA="../plasmid.fasta"
    GENOME="../genome.fasta"
    
    # Belirtilen klasördeki ilk .fastq dosyasını al
    FILTERED_FASTQ=$(ls *.fastq | head -n 1)
    
    ALIGNED_PL="01_aligned_to_plasmid.sam"
    MAPPED_PL="02_mapped_to_plasmid.bam"
    SORTED_MAPPED_PL="03_sorted_mapped_to_plasmid.bam"
    MAPPED_PL_FQ="04_mapped_to_plasmid.fastq"
    ALIGNED_GN="05_aligned_to_other_genome.sam"
    ALIGNED_GN_BM="06_aligned_to_other_genome.bam"
    ALIGNED_GN_ST="07_sorted_aligned_to_other_genome.bam"

    # Run the pipeline using variables
    [ ! -f "$ALIGNED_PL" ] && minimap2 -ax map-ont "$QUERY_FASTA" "$FILTERED_FASTQ" > "$ALIGNED_PL"
    [ ! -f "$MAPPED_PL" ] && samtools view -bS -F 4 "$ALIGNED_PL" > "$MAPPED_PL"
    [ ! -f "$SORTED_MAPPED_PL" ] && samtools sort "$MAPPED_PL" -o "$SORTED_MAPPED_PL"
    [ ! -f "$MAPPED_PL_FQ" ] && samtools fastq "$MAPPED_PL" > "$MAPPED_PL_FQ"
    [ ! -f "$ALIGNED_GN" ] && minimap2 -ax map-ont "$GENOME" "$MAPPED_PL_FQ" > "$ALIGNED_GN"
    [ ! -f "$ALIGNED_GN_BM" ] && samtools view -bS "$ALIGNED_GN" > "$ALIGNED_GN_BM"
    [ ! -f "$ALIGNED_GN_ST" ] && samtools sort "$ALIGNED_GN_BM" -o "$ALIGNED_GN_ST"
    [ ! -f "${dir%/}.fasta" ] && sed -n '1~4s/^@/>/p;2~4p' "$MAPPED_PL_FQ" > "${dir%/}.fasta"

    cd ..
done

