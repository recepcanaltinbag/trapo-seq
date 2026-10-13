#!/usr/bin/env bash
# End-to-end test of trapo-seq on the bundled test data (test-data/).
# Works from any directory. Set PYTHON to choose the interpreter
# (default: python3, or python if python3 is missing).
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"
PYTHON="${PYTHON:-$(command -v python3 || command -v python)}"
CLI=("$PYTHON" "$REPO_DIR/trapo-seq.py")
TD="test-data"

# Prints a step result: message, success (1/0), extra info
print_step() {
    if [ "$2" -eq 1 ]; then
        printf "%-34s (%s) \033[0;32m✓\033[0m\n" "$1" "$3"
    else
        printf "%-34s \033[0;31m✗ %s\033[0m\n" "$1" "$3"
    fi
}

echo "Checking external programs.."
echo "----------------------------"
missing=0
check_tool() {
    local label="$1" tool="$2"; shift 2
    if command -v "$tool" &> /dev/null; then
        print_step "$label" 1 "$("$tool" "$@" 2>&1 | head -n 1)"
    else
        print_step "$label" 0 "(missing)"
        missing=1
    fi
}
check_tool "(1/6) blastn" blastn -version
check_tool "(2/6) makeblastdb" makeblastdb -version
check_tool "(3/6) minimap2" minimap2 --version
check_tool "(4/6) samtools" samtools --version
check_tool "(5/6) mafft" mafft --version
check_tool "(6/6) python" "$PYTHON" --version

echo ""
if [ "$missing" -eq 1 ]; then
    echo "Some required programs are missing. Install them (e.g. 'conda env create -f environment.yml') and rerun."
    exit 1
fi
echo "All required programs are installed."

TESTLOG="$TD/log.txt"
{
    echo "TEST LOG"
    echo "$(date) - Process started"
    "${CLI[@]}" --version
} > "$TESTLOG" 2>&1

# run_step <label> <expected output (glob)> <command...>
run_step() {
    local label="$1" expected="$2"; shift 2
    echo "$(date) - $label: $*" >> "$TESTLOG"
    if ! "$@" >> "$TESTLOG" 2>&1; then
        print_step "$label" 0 "failed"
        echo -e "\033[0;31mSee $TESTLOG for details\033[0m"
        exit 1
    fi
    # shellcheck disable=SC2086
    if ! compgen -G "$expected" > /dev/null; then
        print_step "$label" 0 "no output: $expected"
        exit 1
    fi
    print_step "$label" 1 "$expected"
}

echo -e "\nInput: simulated ONT reads in $TD/barcode*/ (truth in $TD/truth/)\n"

for bc in barcode01 barcode02; do
    run_step "(1/10) Histogram $bc" "$TD/$bc/read_len_hist.pdf" \
        "${CLI[@]}" read_histogram -f "$TD/$bc/${bc}_raw.fastq" -o "$TD/$bc/read_len_hist"
    run_step "(2/10) Filter $bc" "$TD/$bc/$bc.fastq" \
        "${CLI[@]}" filter -f "$TD/$bc/${bc}_raw.fastq" -o "$TD/$bc/$bc.fastq" -l 2000
done
run_step "(3/10) Map" "$TD/barcode01/03_sorted_mapped_to_plasmid.bam" \
    "${CLI[@]}" map -d "$TD" -p "$TD/plasmid.fasta" -g "$TD/genome.fasta" --force
run_step "(4/10) Insert Finder" "$TD/barcode01/insertion_from_bam.tab" \
    "${CLI[@]}" insert_finder_batch -d "$TD"
run_step "(5/10) Annotation of Insertions" "$TD/barcode01/best_alignment.tab" \
    "${CLI[@]}" blast_annot_batch -d "$TD" -g "$TD/genome.fasta" --is_fasta "$TD/IS_curated.fasta" --threads 4
run_step "(6/10) Stats" "$TD/is_stats.rcp" \
    "${CLI[@]}" is_stat -d "$TD" -o "$TD/is_stats.rcp"
run_step "(7/10) Heatmap" "$TD/heatmap.pdf" \
    "${CLI[@]}" heatmap -r "$TD/is_stats.rcp" -o "$TD/heatmap.pdf" -t "$TD/heatmap.tsv"
run_step "(8/10) DR Finder" "$TD/insertions/*.csv" \
    "${CLI[@]}" dr_finder -d "$TD" -p "$TD/plasmid.fasta" --is_fasta "$TD/IS_curated.fasta"
run_step "(9/10) DR Logo" "$TD/dr_logos/*" \
    "${CLI[@]}" dr_logo -d "$TD/insertions" -p "$TD/plasmid.fasta" -o "$TD/dr_logos"
run_step "(10/10) Novel elements" "$TD/novel_elements/novel_elements.tsv" \
    "${CLI[@]}" novel_elements -d "$TD" -g "$TD/genome.fasta" --is_fasta "$TD/IS_curated.fasta" --threads 4

echo -e "\nOptional steps:\n"

BAM="$TD/barcode01/03_sorted_mapped_to_plasmid.bam"
run_step "(1/3) KDE plots" "$TD/barcode01/kde_plot*.pdf" \
    "${CLI[@]}" kde_mobile -b "$BAM" -o kde_plot
run_step "(2/3) In-del plot" "$TD/barcode01/in_del_plot_*weighted.pdf" \
    "${CLI[@]}" in_del_plot -i "$TD/barcode01/insertion_from_bam.tab" -b "$BAM" -o in_del_plot
run_step "(3/3) Map distribution" "$TD/barcode01/map_dist*.pdf" \
    "${CLI[@]}" map_dist -b "$BAM" -o map_dist -f "$TD/barcode01/barcode01_raw.fastq"

echo -e "\nAccuracy against the simulation truth:\n"
run_step "Evaluate" "$TD/evaluation/summary.tsv" \
    "$PYTHON" "$REPO_DIR/simulate/evaluate.py" -s "$TD/truth" -r "$TD" -o "$TD/evaluation" \
    -p "$TD/plasmid.fasta" -i "$TD/IS_curated.fasta" -g "$TD/genome.fasta" -q
# Minimum accuracy expected on this dataset; a regression below these fails the test.
if ! "$PYTHON" - "$TD/evaluation/summary.tsv" <<'PYEOF'
import csv, sys
m = {r["metric"]: r["overall"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
checks = [
    ("detection_finder.recall_unfragmented_eligible", ">=", 0.95),
    ("detection_finder.false_call_rate_empty_plasmid", "<=", 0.01),
    ("position.error.pct_within_20bp", ">=", 95.0),
    ("classification.is_insertion.accuracy", ">=", 0.90),
    ("sites.event_recall", ">=", 0.80),
]
ok = True
for name, op, limit in checks:
    try:
        v = float(m[name])
    except (KeyError, ValueError):
        print(f"  {name}: missing"); ok = False; continue
    good = v >= limit if op == ">=" else v <= limit
    ok &= good
    print(f"  {'PASS' if good else 'FAIL'}  {name} = {v:g} (expected {op} {limit:g})")
sys.exit(0 if ok else 1)
PYEOF
then
    echo -e "\033[0;31mAccuracy below the expected minimum\033[0m"
    exit 1
fi

echo -e "\nAll tests completed.\n"
