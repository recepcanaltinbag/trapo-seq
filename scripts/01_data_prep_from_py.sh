#!/usr/bin/env bash
# Map the reads of every sample folder to the plasmid, then the plasmid-mapped
# reads to the genome, and write mapping statistics.
#
# Usage (run from the data directory that holds one folder per sample):
#   01_data_prep_from_py.sh PLASMID.fasta GENOME.fasta [--force] [--fastq-suffix SUFFIX] [--threads N] [--linear]
#       [--homology-min-len N] [--homology-min-identity F] [--insertion-threshold N]
#
# Each folder <name>/ must contain <name><SUFFIX> (default SUFFIX: ".fastq").
# Outputs written into each folder:
#   02_mapped_to_plasmid.bam, 03_sorted_mapped_to_plasmid.bam(.bai),
#   04_mapped_to_plasmid.fastq, 07_sorted_aligned_to_other_genome.bam(.bai),
#   08_mapping_stats.txt, <name>.fasta
# minimap2 output is streamed straight into samtools, so no SAM files
# (formerly 01_*.sam, 05_*.sam, 06_*.bam) are written.
#
# Circular plasmid (AUDIT B09): the reads are mapped to the plasmid as given
# (minimap2 defaults, as before) and to the plasmid rotated by half its length
# (with --no-long-join, so a short plasmid stretch between an insert and a read
# end is aligned as its own record instead of being lost in a long-joined
# chain). Per read, src/circular_mapping.py keeps the linear alignments unless
# the rotated ones align clearly more of the read, and converts rotated
# alignments back to the ORIGINAL plasmid coordinates (an alignment crossing
# the origin is split there into two records). So 02_/03_*.bam are ordinary
# BAMs against the plasmid as given, but plasmid sequence next to the origin
# is aligned instead of being left in a soft clip.
# --linear restores the plain linear mapping. The Python interpreter used for
# the merge is $TRAPOSEQ_PYTHON (default: python3); it needs pysam.
#
# Plasmid-genome shared segments (GitHub issue #5): before the samples are
# mapped, src/plasmid_genome_homology.py compares the plasmid with the genome
# (minimap2 -x asm20) and writes ./plasmid_genome_homology.tsv (segments of at
# least --homology-min-len bp, default 100, identity >= --homology-min-identity,
# default 0.8), with a WARNING for segments >= --insertion-threshold (default
# 500), which can create false insertion calls. 08_mapping_stats.txt does not
# count a plasmid read as genome-mapped when its genome alignments only cover
# plasmid-aligned sequence (a shared segment), and excludes genomic reads that
# reach the plasmid only through a shared segment; the uncorrected numbers are
# reported too. Both need $PYTHON with pysam; without it the homology table is
# not written and the uncorrected statistics are reported.
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
THREADS=4
CIRCULAR=true
HOMOLOGY_MIN_LEN=100
HOMOLOGY_MIN_IDENTITY=0.8
INSERTION_THRESHOLD=500
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force) FORCE=true; shift ;;
        --fastq-suffix) [[ $# -ge 2 ]] || die "--fastq-suffix needs a value"
                        FASTQ_SUFFIX="$2"; shift 2 ;;
        --threads) [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || die "--threads needs a positive integer"
                   THREADS="$2"; shift 2 ;;
        --linear) CIRCULAR=false; shift ;;
        --homology-min-len) [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || die "--homology-min-len needs a positive integer"
                            HOMOLOGY_MIN_LEN="$2"; shift 2 ;;
        --homology-min-identity) [[ $# -ge 2 && "$2" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]] || die "--homology-min-identity needs a value in [0, 1]"
                                 HOMOLOGY_MIN_IDENTITY="$2"; shift 2 ;;
        --insertion-threshold) [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || die "--insertion-threshold needs a non-negative integer"
                               INSERTION_THRESHOLD="$2"; shift 2 ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done
[[ ${#POSITIONAL[@]} -eq 2 ]] || die "usage: $0 PLASMID.fasta GENOME.fasta [--force] [--fastq-suffix SUFFIX] [--threads N] [--linear] [--homology-min-len N] [--homology-min-identity F] [--insertion-threshold N]"
[[ -f "${POSITIONAL[0]}" ]] || die "plasmid FASTA not found: ${POSITIONAL[0]}"
[[ -f "${POSITIONAL[1]}" ]] || die "genome FASTA not found: ${POSITIONAL[1]}"
QUERY_FASTA=$(abspath "${POSITIONAL[0]}")
GENOME=$(abspath "${POSITIONAL[1]}")

for tool in minimap2 samtools; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool not found in PATH"
done
SCRIPT_DIR=$(cd "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
CIRCULAR_PY="$SCRIPT_DIR/../src/circular_mapping.py"
HOMOLOGY_PY="$SCRIPT_DIR/../src/plasmid_genome_homology.py"
PYTHON=${TRAPOSEQ_PYTHON:-python3}
if [[ "$CIRCULAR" == true ]]; then
    [[ -f "$CIRCULAR_PY" ]] || die "helper not found: $CIRCULAR_PY"
    command -v "$PYTHON" >/dev/null 2>&1 || die "$PYTHON not found (set TRAPOSEQ_PYTHON or use --linear)"
fi
# The homology helpers need Python with pysam; without them the uncorrected
# statistics are written (as before).
HOMOLOGY=false
if [[ -f "$HOMOLOGY_PY" ]] && command -v "$PYTHON" >/dev/null 2>&1 \
        && "$PYTHON" -c "import pysam" >/dev/null 2>&1; then
    HOMOLOGY=true
else
    echo "WARNING: $PYTHON with pysam not available: plasmid-genome homology is not checked and" \
         "08_mapping_stats.txt is not corrected for shared segments." >&2
fi

# ---------------------------------------------------------- temp-file helpers
CURRENT_TMP=""
WORK_TMPDIR=""
cleanup() {
    if [[ -n "$CURRENT_TMP" ]]; then rm -f -- "$CURRENT_TMP"; fi
    if [[ -n "$WORK_TMPDIR" ]]; then rm -rf -- "$WORK_TMPDIR"; fi
}
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
    if [[ "$HOMOLOGY" == true ]]; then
        # Corrected for plasmid-genome shared segments (see the header).
        tmp=$(tmp_for "$output_txt"); CURRENT_TMP=$tmp
        "$PYTHON" "$HOMOLOGY_PY" mapstats "$plasmid_bam" "$genome_bam" "$total" "$tmp" \
            --genome "$GENOME" --threads "$THREADS" >/dev/null
        commit "$tmp" "$output_txt"
        cat -- "$output_txt"
        echo "Mapping statistics written to $output_txt"
        return 0
    fi
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

    local MAPPED_PL="02_mapped_to_plasmid.bam"
    local SORTED_MAPPED_PL="03_sorted_mapped_to_plasmid.bam"
    local MAPPED_PL_FQ="04_mapped_to_plasmid.fastq"
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

    echo "STEP 1: mapping to the plasmid (minimap2 $(minimap2 --version)), keeping mapped reads"
    if needs_build "$MAPPED_PL"; then
        tmp=$(tmp_for "$MAPPED_PL"); CURRENT_TMP=$tmp
        if [[ "$CIRCULAR" == true ]]; then
            # Both mappings keep unmapped records so the two streams stay in step.
            minimap2 -t "$THREADS" -ax map-ont "$QUERY_FASTA" "$fastq" \
                | samtools view -@ "$THREADS" -u -o "$WORK_TMPDIR/linear.bam" -
            minimap2 -t "$THREADS" -ax map-ont --no-long-join "$ROTATED_FASTA" "$fastq" \
                | samtools view -@ "$THREADS" -u -o "$WORK_TMPDIR/rotated.bam" -
            "$PYTHON" "$CIRCULAR_PY" merge "$QUERY_FASTA" \
                "$WORK_TMPDIR/linear.bam" "$WORK_TMPDIR/rotated.bam" "$tmp"
            rm -f -- "$WORK_TMPDIR/linear.bam" "$WORK_TMPDIR/rotated.bam"
        else
            minimap2 -t "$THREADS" -ax map-ont "$QUERY_FASTA" "$fastq" \
                | samtools view -@ "$THREADS" -b -F 4 -o "$tmp" -
        fi
        commit "$tmp" "$MAPPED_PL"
    fi

    echo "STEP 2: sorting, indexing"
    if needs_build "$SORTED_MAPPED_PL"; then
        tmp=$(tmp_for "$SORTED_MAPPED_PL"); CURRENT_TMP=$tmp
        samtools sort -@ "$THREADS" -o "$tmp" "$MAPPED_PL"
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
    if needs_build "$ALIGNED_GN_ST"; then
        tmp=$(tmp_for "$ALIGNED_GN_ST"); CURRENT_TMP=$tmp
        minimap2 -t "$THREADS" -ax map-ont "$GENOME" "$MAPPED_PL_FQ" \
            | samtools sort -@ "$THREADS" -o "$tmp" -
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
if [[ "$CIRCULAR" == true ]]; then
    WORK_TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/trapo-seq-map.XXXXXX")
    ROTATED_FASTA="$WORK_TMPDIR/plasmid_rotated.fasta"
    "$PYTHON" "$CIRCULAR_PY" rotate "$QUERY_FASTA" "$ROTATED_FASTA"
fi
if [[ "$HOMOLOGY" == true ]]; then
    echo "Plasmid-genome shared segments (minimap2 -x asm20, >= $HOMOLOGY_MIN_LEN bp, identity >= $HOMOLOGY_MIN_IDENTITY)"
    "$PYTHON" "$HOMOLOGY_PY" detect "$QUERY_FASTA" "$GENOME" "$PWD/plasmid_genome_homology.tsv" \
        --min-len "$HOMOLOGY_MIN_LEN" --min-identity "$HOMOLOGY_MIN_IDENTITY" \
        --warn-len "$INSERTION_THRESHOLD" --threads "$THREADS"
else
    rm -f -- plasmid_genome_homology.tsv     # never leave a table from an earlier run
fi
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
