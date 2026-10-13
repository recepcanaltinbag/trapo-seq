#!/usr/bin/env bash
# Map the reads of every sample folder to the plasmid, then the plasmid-mapped
# reads to the genome, and write mapping statistics.
#
# Usage (run from the data directory that holds one folder per sample):
#   01_data_prep_from_py.sh PLASMID.fasta GENOME.fasta [--force] [--fastq-suffix SUFFIX]
#
# Each folder <name>/ must contain <name><SUFFIX> (default SUFFIX: ".fastq").
# Outputs written into each folder:
#   01_aligned_to_plasmid.sam, 02_mapped_to_plasmid.bam,
#   03_sorted_mapped_to_plasmid.bam(.bai), 04_mapped_to_plasmid.fastq,
#   05_aligned_to_other_genome.sam, 06_aligned_to_other_genome.bam,
#   07_sorted_aligned_to_other_genome.bam(.bai), 08_mapping_stats.txt,
#   <name>.fasta
#
# Every output is written to a temporary file and renamed only when the step
# succeeded, so a failed run never leaves a partial file that a later run
# would reuse. Any failing command aborts the script with a non-zero status.
# Once a step is (re)built, all later steps of that folder are rebuilt too.

set -Eeuo pipefail
shopt -s nullglob

die() { echo "ERROR: $*" >&2; exit 1; }

abspath() {
    # Absolute path of an existing file (portable; no `realpath -f` on macOS).
    local d
    d=$(cd "$(dirname -- "$1")" && pwd -P) || return 1
    printf '%s/%s\n' "$d" "$(basename -- "$1")"
}

# ---------------------------------------------------------------- arguments
FORCE=false
FASTQ_SUFFIX=".fastq"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force) FORCE=true; shift ;;
        --fastq-suffix) [[ $# -ge 2 ]] || die "--fastq-suffix needs a value"
                        FASTQ_SUFFIX="$2"; shift 2 ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done
[[ ${#POSITIONAL[@]} -eq 2 ]] || die "usage: $0 PLASMID.fasta GENOME.fasta [--force] [--fastq-suffix SUFFIX]"
[[ -f "${POSITIONAL[0]}" ]] || die "plasmid FASTA not found: ${POSITIONAL[0]}"
[[ -f "${POSITIONAL[1]}" ]] || die "genome FASTA not found: ${POSITIONAL[1]}"
QUERY_FASTA=$(abspath "${POSITIONAL[0]}")
GENOME=$(abspath "${POSITIONAL[1]}")

for tool in minimap2 samtools; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool not found in PATH"
done

# ---------------------------------------------------------- temp-file helpers
CURRENT_TMP=""
cleanup() { if [[ -n "$CURRENT_TMP" ]]; then rm -f -- "$CURRENT_TMP"; fi; }
trap cleanup EXIT
trap 'echo "ERROR: command failed (line $LINENO) in $PWD" >&2' ERR

# tmp_for OUT -> prints a temporary name in the same folder (same extension,
# so samtools picks the right output format).
tmp_for() { printf '.tmp_%s_%s' "$$" "$1"; }

# commit TMP OUT -> atomically move a finished temp file into place.
commit() { mv -f -- "$1" "$2"; CURRENT_TMP=""; }

REBUILT=false
# needs_build OUT -> true if OUT must be (re)built. Once anything was rebuilt in
# this folder, every later output is rebuilt as well (it may depend on it).
needs_build() {
    local out=$1
    if [[ "$FORCE" == true ]]; then
        [[ -f "$out" ]] && echo "  $out exists; overwriting (--force)."
        REBUILT=true; return 0
    fi
    if [[ "$REBUILT" == true ]]; then
        [[ -f "$out" ]] && echo "  $out exists but an upstream file was rebuilt; rebuilding."
        return 0
    fi
    if [[ -f "$out" ]]; then
        echo "  $out exists; reusing (use --force to overwrite)."
        return 1
    fi
    REBUILT=true
    return 0
}

# ------------------------------------------------------------------ steps
count_reads_in_fastq() {
    local lines
    lines=$(wc -l < "$1")
    echo $(( lines / 4 ))
}

count_mapped_read_names() {
    # Distinct read names with at least one primary/supplementary alignment.
    samtools view -F 0x104 "$1" | cut -f1 | sort -u | wc -l | tr -d ' '
}

summarize_mappings() {
    local plasmid_bam=$1 genome_bam=$2 output_txt=$3 fastq_file=$4
    local total plasmid both only_plasmid rate tmp

    total=$(count_reads_in_fastq "$fastq_file")
    plasmid=$(count_mapped_read_names "$plasmid_bam")
    both=$(count_mapped_read_names "$genome_bam")
    only_plasmid=$(( plasmid - both ))
    # Reads mapped to the genome are a subset of the plasmid-mapped reads, so
    # the rate is both / plasmid. No plasmid reads -> NA (no division by zero).
    rate=$(awk -v b="$both" -v p="$plasmid" 'BEGIN { if (p > 0) printf "%.2f", b / p; else printf "NA" }')

    tmp=$(tmp_for "$output_txt"); CURRENT_TMP=$tmp
    {
        echo "Total reads: $total"
        echo "Plasmid: $plasmid"
        echo "Only Plasmid: $only_plasmid"
        echo "Both plasmid and genome: $both"
        echo "Transposition Rate: $rate"
    } > "$tmp"
    commit "$tmp" "$output_txt"
    cat -- "$output_txt"
    echo "Mapping statistics written to $output_txt"
}

process_folder() {
    local name=$1
    local fastq="${name}${FASTQ_SUFFIX}"
    local tmp

    local ALIGNED_PL="01_aligned_to_plasmid.sam"
    local MAPPED_PL="02_mapped_to_plasmid.bam"
    local SORTED_MAPPED_PL="03_sorted_mapped_to_plasmid.bam"
    local MAPPED_PL_FQ="04_mapped_to_plasmid.fastq"
    local ALIGNED_GN="05_aligned_to_other_genome.sam"
    local ALIGNED_GN_BM="06_aligned_to_other_genome.bam"
    local ALIGNED_GN_ST="07_sorted_aligned_to_other_genome.bam"
    local OUTPUT_STAT="08_mapping_stats.txt"
    local FASTA_OUT="${name}.fasta"

    if [[ ! -f "$fastq" ]]; then
        echo "WARNING: no $fastq in $name/ (expected <folder>/<folder>${FASTQ_SUFFIX}). Skipping." >&2
        return 0
    fi

    echo ""
    echo "--------- $name ---------"
    echo "Reads:   $fastq"
    echo "Plasmid: $QUERY_FASTA"
    echo "Genome:  $GENOME"
    REBUILT=false

    echo "STEP 1: mapping to the plasmid (minimap2 $(minimap2 --version))"
    if needs_build "$ALIGNED_PL"; then
        tmp=$(tmp_for "$ALIGNED_PL"); CURRENT_TMP=$tmp
        minimap2 -ax map-ont "$QUERY_FASTA" "$fastq" > "$tmp"
        commit "$tmp" "$ALIGNED_PL"
    fi

    echo "STEP 2: keeping mapped reads, sorting, indexing"
    if needs_build "$MAPPED_PL"; then
        tmp=$(tmp_for "$MAPPED_PL"); CURRENT_TMP=$tmp
        samtools view -b -F 4 -o "$tmp" "$ALIGNED_PL"
        commit "$tmp" "$MAPPED_PL"
    fi
    if needs_build "$SORTED_MAPPED_PL"; then
        tmp=$(tmp_for "$SORTED_MAPPED_PL"); CURRENT_TMP=$tmp
        samtools sort -o "$tmp" "$MAPPED_PL"
        commit "$tmp" "$SORTED_MAPPED_PL"
    fi
    if [[ "$REBUILT" == true || ! -f "$SORTED_MAPPED_PL.bai" ]]; then
        tmp=$(tmp_for "$SORTED_MAPPED_PL.bai"); CURRENT_TMP=$tmp
        samtools index "$SORTED_MAPPED_PL" "$tmp"
        commit "$tmp" "$SORTED_MAPPED_PL.bai"
    fi
    if needs_build "$MAPPED_PL_FQ"; then
        # Primary records only (one per read), in the original read orientation.
        tmp=$(tmp_for "$MAPPED_PL_FQ"); CURRENT_TMP=$tmp
        samtools fastq -F 0x904 "$MAPPED_PL" > "$tmp"
        commit "$tmp" "$MAPPED_PL_FQ"
    fi

    echo "STEP 3: mapping plasmid reads to the genome"
    if needs_build "$ALIGNED_GN"; then
        tmp=$(tmp_for "$ALIGNED_GN"); CURRENT_TMP=$tmp
        minimap2 -ax map-ont "$GENOME" "$MAPPED_PL_FQ" > "$tmp"
        commit "$tmp" "$ALIGNED_GN"
    fi
    if needs_build "$ALIGNED_GN_BM"; then
        tmp=$(tmp_for "$ALIGNED_GN_BM"); CURRENT_TMP=$tmp
        samtools view -b -o "$tmp" "$ALIGNED_GN"
        commit "$tmp" "$ALIGNED_GN_BM"
    fi
    if needs_build "$ALIGNED_GN_ST"; then
        tmp=$(tmp_for "$ALIGNED_GN_ST"); CURRENT_TMP=$tmp
        samtools sort -o "$tmp" "$ALIGNED_GN_BM"
        commit "$tmp" "$ALIGNED_GN_ST"
    fi
    if [[ "$REBUILT" == true || ! -f "$ALIGNED_GN_ST.bai" ]]; then
        tmp=$(tmp_for "$ALIGNED_GN_ST.bai"); CURRENT_TMP=$tmp
        samtools index "$ALIGNED_GN_ST" "$tmp"
        commit "$tmp" "$ALIGNED_GN_ST.bai"
    fi

    echo "STEP 4: FASTA of the plasmid-mapped reads"
    if needs_build "$FASTA_OUT"; then
        # Same records as $MAPPED_PL_FQ (primary only, original orientation).
        tmp=$(tmp_for "$FASTA_OUT"); CURRENT_TMP=$tmp
        samtools fasta -F 0x904 "$MAPPED_PL" > "$tmp"
        commit "$tmp" "$FASTA_OUT"
    fi

    echo "Summarizing the mapping stats to $OUTPUT_STAT"
    summarize_mappings "$SORTED_MAPPED_PL" "$ALIGNED_GN_ST" "$OUTPUT_STAT" "$fastq"
}

# ------------------------------------------------------------------- main
found=0
for dir in */; do
    name=${dir%/}
    found=$((found + 1))
    pushd -- "$name" >/dev/null
    process_folder "$name"
    popd >/dev/null
done
[[ $found -gt 0 ]] || die "no sample folders found in $PWD"
echo "Mapping finished for $found folder(s)."
