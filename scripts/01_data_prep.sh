#!/usr/bin/env bash
# Stand-alone data preparation (without the Python CLI).
#
# Run from the data directory, which must contain plasmid.fasta, genome.fasta
# and one folder per sample (barcode folders renamed to sample/condition names).
#
# 1 - FOLDER PREPARATION
#   In every folder, *.fastq.gz / *.fq.gz files are decompressed and all
#   FASTQ files are concatenated into <folder>_merged.fastq.
#
# 2 - ALIGNMENT
#   Reads -> plasmid -> plasmid-mapped reads -> genome, via
#   01_data_prep_from_py.sh (same outputs as `trapo-seq.py map`).
#
# Statistics of reads:
#   Total Reads: A
#   Plasmid Only: B
#   Genome and Plasmid: C
#     -> Transposition Ratio: C/(C+B)
#
# Usage: 01_data_prep.sh [--force]

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR=$(cd "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
QUERY_FASTA="$PWD/plasmid.fasta"
GENOME="$PWD/genome.fasta"
[[ -f "$QUERY_FASTA" ]] || { echo "ERROR: $QUERY_FASTA not found" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: $GENOME not found" >&2; exit 1; }

for dir in */; do
    name=${dir%/}
    merged="${name}_merged.fastq"
    (
        cd -- "$name"
        for gz in *.fastq.gz *.fq.gz; do
            gunzip -- "$gz"
        done
        inputs=()
        for f in *.fastq *.fq; do
            [[ "$f" == "$merged" || "$f" == .tmp_* ]] && continue
            inputs+=("$f")
        done
        if [[ ${#inputs[@]} -eq 0 ]]; then
            echo "WARNING: no FASTQ files in $name/" >&2
            exit 0
        fi
        if [[ ! -f "$merged" ]]; then
            tmp=".tmp_$$_${merged}"
            trap 'rm -f -- "$tmp"' EXIT
            cat -- "${inputs[@]}" > "$tmp"
            mv -f -- "$tmp" "$merged"
        fi
    )
done

exec bash "$SCRIPT_DIR/01_data_prep_from_py.sh" "$QUERY_FASTA" "$GENOME" \
    --fastq-suffix _merged.fastq "$@"
