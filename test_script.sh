#!/bin/bash


# Function to print with green tick, version, and alignment
print_step() {
    local step_message="$1"
    local success="$2"
    local version_info="$3"
    
    # Prints the message, version info, and the corresponding green tick or red cross aligned on the right side
    if [ "$success" -eq 1 ]; then
        printf "%-30s (%s) \033[0;32m✓\033[0m\n" "$step_message" "$version_info"  # Green checkmark
    else
        printf "%-30s \033[0;31m✗ (missing)\033[0m\n" "$step_message"  # Red cross
    fi
}


programs=True

echo "Controlling externals.."
echo "----------------------"

# Start of step 1: Check if BLAST is installed and display version
if command -v blastn &> /dev/null; then
    blast_version=$(blastn -version  2>&1  | head -n 1)  # Get the first line of the version info
    makeblastdb_version=$(makeblastdb -version   2>&1  | head -n 1)
    print_step "(1/6) BLAST was installed" 1 "$blast_version, $makeblastdb_version"
else
    print_step "(1/6) Warning: BLAST is not installed." 0 ""
    programs=False
fi

# Step 3: Check if minimap2 is installed and display version
if command -v minimap2 &> /dev/null; then
    minimap2_version=$(minimap2 --version  2>&1  | head -n 1)  # Get the first line of version info
    print_step "(3/6) minimap2 was installed" 1 "$minimap2_version"
else
    print_step "(3/6) Warning: minimap2 is not installed." 0 ""
    programs=False
fi

# Step 4: Check if samtools is installed and display version
if command -v samtools &> /dev/null; then
    samtools_version=$(samtools --version  2>&1  | head -n 1)  # Get the version info
    print_step "(4/6) samtools was installed" 1 "$samtools_version"
else
    print_step "(4/6) Warning: samtools is not installed." 0 ""
    programs=False
fi

# Step 5: Check if mafft is installed and display version
if command -v mafft &> /dev/null; then
    mafft_version=$(mafft --version 2>&1  | head -n 1 )  # Get the version info
    print_step "(5/6) mafft was installed" 1 "$mafft_version"
else
    print_step "(5/6) Warning: mafft is not installed." 0 ""
    programs=False
fi

# Step 6: Check if Python is installed and display version
if command -v python &> /dev/null; then
    python_version=$(python --version  2>&1 )  # Get the version info
    print_step "(6/6) Python was installed" 1 "$python_version"
else
    print_step "(6/6) Warning: Python is not installed." 0 ""
    programs=False
fi

echo ""
if [ "$programs" = False ]; then
    echo "Warning: Some required programs are missing. Please install the missing programs and rerun the script."
    # Kullanıcıya devam etmek isteyip istemediğini soruyoruz
    read -p "Do you want to continue despite missing programs? (yes/no): " user_choice
    [[ "$user_choice" == "yes" || "$user_choice" == "y" ]] && echo "Continuing..." || { echo "Exiting..."; exit 0; }
else
    echo "All required programs are installed."
fi


TESTLOG="test-data/log.txt"

echo -e "\nTEST LOG\n" > "$TESTLOG"
echo "$(date) - Process started" >> "$TESTLOG"
echo -e "\ncreated by Recep Can Altınbağ \n" >> "$TESTLOG"

#YOU CAN USE THIS IS A CHEAT SHEET FOR YOUR ANAYLSIS:

python trapo-seq.py --version >> "$TESTLOG" 2>&1

echo -e "\n Input: test-data/*/*_raw.fastq \n"

print_step "(1/9) Histogram" 1 "test-data/barcode01/read_len_hist.pdf"
python trapo-seq.py read_histogram -f test-data/barcode01/barcode01_raw.fastq -o test-data/barcode01/read_len_hist >> "$TESTLOG" 2>&1

print_step "(2/9) Filter" 1 "test-data/barcode01/barcode01.fastq"
python trapo-seq.py filter -f test-data/barcode01/barcode01_raw.fastq -o test-data/barcode01/barcode01.fastq -l 2000 >> "$TESTLOG" 2>&1

print_step "(3/9) Map" 1 "test-data/*/03.*.bam"
python trapo-seq.py map -d test-data -p test-data/plasmid.fasta -g test-data/genome.fasta --force >> "$TESTLOG" 2>&1

print_step "(4/9) Insert Finder" 1 "test-data/*/*insertion_from_bam.tab"
python trapo-seq.py insert_finder_batch -d test-data >> "$TESTLOG" 2>&1

print_step "(5/9) Annotation of Insertions" 1 "test-data/*/*best_alignment.tab"
python trapo-seq.py blast_annot_batch -d test-data -g  test-data/genome.fasta --is_fasta test-data/IS_curated.fasta --threads 4 >> "$TESTLOG" 2>&1

print_step "(6/9) Stats" 1 "test-data/is_stats.rcp"
python trapo-seq.py is_stat -d test-data -o test-data/is_stats.rcp >> "$TESTLOG" 2>&1

print_step "(7/9) Heatmap" 1 "test-data/heatmap.pdf"
python trapo-seq.py heatmap -r test-data/is_stats.rcp -o test-data/heatmap.pdf -t test-data/heatmap.tsv >> "$TESTLOG" 2>&1

print_step "(8/9) DR Finder" 1 "test-data/insertions/*.csv"
python trapo-seq.py dr_finder -d test-data -p test-data/plasmid.fasta >> "$TESTLOG" 2>&1

print_step "(9/9) DR Logo" 1 "test-data/dr_logos/*.pdf"
python trapo-seq.py dr_logo -d data/insertions -p test-data/plasmid.fasta -o test-data/dr_logos >> "$TESTLOG" 2>&1

echo -e "\n Optional Steps: \n"

print_step "(1/3) Kde Plots" 1 "test-data/barcode01/kde_plot.pdf"
python trapo-seq.py kde_mobile -b test-data/barcode01/03_sorted_mapped_to_plasmid.bam -o kde_plot >> "$TESTLOG" 2>&1

print_step "(2/3) In-del Plot" 1 "test-data/barcode01/in_del_plot.pdf"
python trapo-seq.py kde_mobile -i test-data/barcode01/insertion_from_bam.tab -b test-data/barcode01/03_sorted_mapped_to_plasmid.bam -o kde_plot >> "$TESTLOG" 2>&1

print_step "(3/3) Map Dist Plot" 1 "test-data/barcode01/map_dist.pdf"
python trapo-seq.py map_dist -b test-data/barcode01/03_sorted_mapped_to_plasmid.bam -o map_dist >> "$TESTLOG" 2>&1


# Additional checks or commands
echo ""
echo "Tests completed."