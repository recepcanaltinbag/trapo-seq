#!/bin/bash


for dir in */; do
    cd "$dir" || continue

    QUERY_FASTA="$1"
    GENOME="$2"
    

    FORCE=false  # Default is not force

    # Parse arguments (look for --force flag)
    for arg in "$@"; do
        if [[ "$arg" == "--force" ]]; then
            FORCE=true
            break
        fi
    done


    # Belirtilen klasördeki ilk .fastq dosyasını al
    FILTERED_FASTQ=$(ls ${dir%/}.fastq | head -n 1)

    if [[ -z "$FILTERED_FASTQ" ]]; then
        echo -e "\033[31;1mWARNING:\033[0m"
        echo -e "No $dir.fastq file found in \033[31;1m$dir\033[0m. Skipping..."
        echo -e "\033[31;1mIt is expected exact name of fastq file with the folder name such as <folder>/<folder>.fastq) \033[0m. Skipping..."
        cd ..  # Return to the parent directory
        continue
    fi

    ALIGNED_PL="01_aligned_to_plasmid.sam"
    MAPPED_PL="02_mapped_to_plasmid.bam"
    SORTED_MAPPED_PL="03_sorted_mapped_to_plasmid.bam"
    MAPPED_PL_FQ="04_mapped_to_plasmid.fastq"
    ALIGNED_GN="05_aligned_to_other_genome.sam"
    ALIGNED_GN_BM="06_aligned_to_other_genome.bam"
    ALIGNED_GN_ST="07_sorted_aligned_to_other_genome.bam"

    echo ""
    echo -e "\033[33;1m---------INPUTS----------\033[0m"
    echo -e "\033[33;1mProcessing file:\033[0m $FILTERED_FASTQ"
    echo -e "\033[33;1mQUERY_FASTA: $QUERY_FASTA\033[0m"
    echo -e "\033[33;1mGENOME: $GENOME\033[0m"
    echo -e "\033[33;1m-------------------------\033[0m"

    # Run the pipeline using variables
    echo "minimap2 with version: "
    minimap2 --version
    echo -e "\033[33;1m------------STEP 1---------------\033[0m"
    echo "Mapping the plasmid with minimap2"

    # Check if file exists and handle the force flag
    if [ -f "$ALIGNED_PL" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $ALIGNED_PL already exists. Overwriting due to --force flag."
            minimap2 -ax map-ont "$QUERY_FASTA" "$FILTERED_FASTQ" > "$ALIGNED_PL"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $ALIGNED_PL already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$ALIGNED_PL does not exist. Proceeding with the operation."
        minimap2 -ax map-ont "$QUERY_FASTA" "$FILTERED_FASTQ" > "$ALIGNED_PL"
        # Perform the actual operation (e.g., creating or processing the file)
    fi

    echo -e "\033[33;1m------------STEP 2---------------\033[0m"
    echo "Samtools..."

    if [ -f "$MAPPED_PL" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $MAPPED_PL already exists. Overwriting due to --force flag."
            samtools view -bS -F 4 "$ALIGNED_PL" > "$MAPPED_PL"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $MAPPED_PL already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$MAPPED_PL does not exist. Proceeding with the operation."
        samtools view -bS -F 4 "$ALIGNED_PL" > "$MAPPED_PL"
        # Perform the actual operation (e.g., creating or processing the file)
    fi


    if [ -f "$SORTED_MAPPED_PL" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $SORTED_MAPPED_PL already exists. Overwriting due to --force flag."
            samtools sort "$MAPPED_PL" -o "$SORTED_MAPPED_PL"
            echo "Samtools indexing the " "$SORTED_MAPPED_PL"
            samtools index "$SORTED_MAPPED_PL"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $SORTED_MAPPED_PL already exists. Use --force to overwrite."
            echo "Samtools indexing the " "$SORTED_MAPPED_PL"
            samtools index "$SORTED_MAPPED_PL"
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$SORTED_MAPPED_PL does not exist. Proceeding with the operation."
        samtools sort "$MAPPED_PL" -o "$SORTED_MAPPED_PL"
        echo "Samtools indexing the " "$SORTED_MAPPED_PL"
        samtools index "$SORTED_MAPPED_PL"
        # Perform the actual operation (e.g., creating or processing the file)
    fi


    if [ -f "$MAPPED_PL_FQ" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $MAPPED_PL_FQ already exists. Overwriting due to --force flag."
            samtools fastq "$MAPPED_PL" > "$MAPPED_PL_FQ"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $MAPPED_PL_FQ already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$MAPPED_PL_FQ does not exist. Proceeding with the operation."
        samtools fastq "$MAPPED_PL" > "$MAPPED_PL_FQ"
        # Perform the actual operation (e.g., creating or processing the file)
    fi
    
    
    echo -e "\033[33;1m------------STEP 3---------------\033[0m"

    echo "Mapping to genome with minimap2"

    if [ -f "$ALIGNED_GN" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $ALIGNED_GN already exists. Overwriting due to --force flag."
            minimap2 -ax map-ont "$GENOME" "$MAPPED_PL_FQ" > "$ALIGNED_GN"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $ALIGNED_GN already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$ALIGNED_GN does not exist. Proceeding with the operation."
        minimap2 -ax map-ont "$GENOME" "$MAPPED_PL_FQ" > "$ALIGNED_GN"
        # Perform the actual operation (e.g., creating or processing the file)
    fi

    if [ -f "$ALIGNED_GN_BM" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $ALIGNED_GN_BM already exists. Overwriting due to --force flag."
            samtools view -bS "$ALIGNED_GN" > "$ALIGNED_GN_BM"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $ALIGNED_GN_BM already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$ALIGNED_GN_BM does not exist. Proceeding with the operation."
        samtools view -bS "$ALIGNED_GN" > "$ALIGNED_GN_BM"
        # Perform the actual operation (e.g., creating or processing the file)
    fi

    if [ -f "$ALIGNED_GN_ST" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: $ALIGNED_GN_ST already exists. Overwriting due to --force flag."
            samtools sort "$ALIGNED_GN_BM" -o "$ALIGNED_GN_ST"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: $ALIGNED_GN_ST already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "$ALIGNED_GN_ST does not exist. Proceeding with the operation."
        samtools sort "$ALIGNED_GN_BM" -o "$ALIGNED_GN_ST"
        # Perform the actual operation (e.g., creating or processing the file)
    fi

    echo -e "\033[33;1m------------STEP FINAL---------------\033[0m"


    if [ -f "${dir%/}.fasta" ]; then
        if [ "$FORCE" = true ]; then
            echo "Warning: fasta already exists. Overwriting due to --force flag."
            sed -n '1~4s/^@/>/p;2~4p' "$MAPPED_PL_FQ" > "${dir%/}.fasta"
            # Proceed with the operation (e.g., overwrite or create new file)
        else
            echo "Warning: fasta already exists. Use --force to overwrite."
        fi
    else
        # File doesn't exist, proceed with the operation
        echo "fasta does not exist. Proceeding with the operation."
        sed -n '1~4s/^@/>/p;2~4p' "$MAPPED_PL_FQ" > "${dir%/}.fasta"
        # Perform the actual operation (e.g., creating or processing the file)
    fi

    cd ..
done