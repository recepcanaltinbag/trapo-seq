# Code audit — trapo-seq v0.1.0 (main @ 980d5e6)

Audit date: 2026-10-13. Every finding has an ID that is referenced from
`ROADMAP.md` and `CHANGELOG.md` when it is fixed. Batch modules
(`*_batch.py`) are copies of their single-mode counterparts, so a finding in
one applies to both unless stated otherwise.

Severity: **S** = scientifically wrong results, **C** = crash / pipeline stops,
**U** = usability / portability, **Q** = code quality.

## Mapping and insertion finding

| ID | Sev | Location | Finding |
|----|-----|----------|---------|
| B01 | C | `scripts/01_data_prep*.sh` | `sed '1~4…'` is GNU-only; on macOS an empty `<sample>.fasta` is written and `blast_annot_batch` later crashes. |
| B02 | S | `e_insert_finder_from_bam.py:109` | `reads[read_id] = insert_list` overwrites CIGAR insertions from earlier alignment records of the same read. |
| B03 | S | `e_insert_finder_from_bam.py:20-31,89-119` | Secondary alignments are not skipped: quality becomes 0, coverage is inflated, insertions duplicated. |
| B04 | S | `e_insert_finder_from_bam.py:274-276` | Strand of soft-clip (SC) insertions comes from a leftover loop variable, not the flanking alignment. |
| B05 | S | `e_insert_finder_from_bam.py:69-79` | `=`/`X` CIGAR ops are ignored, so coordinates break with `--eqx` alignments. |
| B06 | C | `a_data_mapping.py`, `01_data_prep_from_py.sh` | No `set -e`; mapping failures exit 0 and leave partial SAM/BAM that later runs silently reuse. |
| B07 | C | `01_data_prep_from_py.sh:33` | Division by zero when no plasmid reads map. |
| B08 | S | `e_insert_finder_from_bam.py` | *(found during fixing)* For gaps next to the plasmid origin, the flank ending at the linear reference end was used, giving ref_pos 0 / plasmid length instead of the true site. |
| B09 | S | mapping (design) | *(found during fixing, **fixed**)* The plasmid is mapped as a linear reference: plasmid sequence just past the origin is not aligned, so the extracted insert includes ~hundreds of bp of plasmid. Needs a circular-aware reference (rotated/doubled plasmid). |

## Annotation (BLAST)

| ID | Sev | Location | Finding |
|----|-----|----------|---------|
| B10 | S | `g_blast_annot.py:89-91` | Only HSPs with **E-value == 0** are kept (≈ < 1e-180). Real, shorter or lower-identity inserts become "no blast hit" → labelled *Contamination*. |
| B11 | S | `g_blast_annot.py:255` | Genome row keeps the *last* (weakest) HSP instead of the best one. |
| B12 | S | `g_blast_annot.py:258-275` | One insertion can emit several duplicate/contradictory rows (IS_DB and Genome); downstream counts are inflated. |
| B13 | S | `g_blast_annot.py:113,129` | Subject ID taken as the *last word* of the FASTA header (`…complete genome` → `genome`). |
| B14 | S | `g_blast_annot.py:94,260-275` | Off-by-one in coverage / "Explained" (`abs(e-s)` without `+1`; minus-strand 2 bp short). |
| B15 | C | `trapo-seq.py:68,181,187-190` | `--is_fasta` fallback to the dummy file is computed but never passed; single mode crashes with `None`. |
| B16 | U | `trapo-seq.py:74,86` | `--no_temp` (store_false) does the opposite of what the README says. |
| B17 | S | `g_blast_annot.py:206-217` | BLAST DBs have fixed names and are never rebuilt when the FASTA changes → stale results. |
| B18 | C | `g_blast_annot.py:174,214,217` | `os.system` with unquoted paths, return codes ignored (breaks on paths with spaces). |
| B19 | Q | `g_blast_annot.py:60-64` | Whole FASTA re-parsed per insertion (O(N_ins × N_reads)); one blastn call per insertion. |
| B20 | U | `trapo-seq.py:70` | Single-mode output `best_alignments.tab` (plural) is ignored by `is_stat`/`dr_finder`, which expect `best_alignment.tab`. |

## Direct repeats (TSD) and logos

| ID | Sev | Location | Finding |
|----|-----|----------|---------|
| B30 | S | `j_dr_finder_multithread.py:137-161` | Single-flank branch pairs the wrong flank ends, mixes read and plasmid coordinates, no distance threshold. |
| B31 | S | `j_dr_finder_multithread.py:225-238` | Circular plasmid not handled: DRs near the origin are lost or mis-sized. |
| B32 | C | `j_dr_finder_multithread.py:389` | `pd.concat([])` crash when a sample yields no DR rows. |
| B33 | U | `trapo-seq.py:102,208` | `dr_finder --gap` is accepted but never used. |
| B34 | S | `j_dr_finder_multithread.py:367-369` | Read IDs containing `_` are mis-parsed (`split('_')` instead of `rsplit('_', 2)`). |
| B35 | S | `j_dr_finder_multithread.py:226` | Strict `<` makes a DR equal to the threshold count as "none"; zero/negative overlaps indistinguishable. |
| B36 | S | `k_dr_logo.py:114-127` | Alignment truncated from column 0 by modal length; leading gap columns shift the logo. |
| B37 | C | `k_dr_logo.py:48-63,181` | MAFFT errors unhandled; `pd.concat([])` crash when no input; leaked file handles; stray figures. |
| B38 | Q | `j_dr_finder.py` | Dead duplicate of the "multithread" module (which has no threads). |

## Statistics and plots

| ID | Sev | Location | Finding |
|----|-----|----------|---------|
| B40 | C | `h_is_stats.py:30` | `str.split('_', expand=True)` crashes on names with ≠ 1 underscore or no IS hits. |
| B41 | S | `h_is_stats.py:39-50` | Counts BLAST rows, not insertions; denominator excludes Genome/Contamination without saying so. |
| B42 | C | `i_heatmap.py:117-126` | Crash when every IS falls on one side of the 15 % split. |
| B43 | S | `f_in_del_plot.py:156-157` | `dict.update` overwrites SC counts with IN counts instead of adding them. |
| B44 | S | `d_kde_mobile.py:155-183` | Hard-clipped supplementary alignments use wrong read coordinates → genome segments dropped. |
| B45 | S | `c_map_dist.py:53-75,175` | Multi-contig genomes: coverage of all contigs summed into contig 1's array. |
| B46 | C | `d_kde_mobile.py:258` | Crash when no read passes `--plasmid_threshold`. |
| B47 | U | `trapo-seq.py:144-148` | `read_histogram -v` never draws the histogram. |
| B48 | C | `a_read_histograms.py:51,70` | Zero bin width / empty FASTQ crash; `-o x.pdf` → `x.pdf.pdf`. |

## Portability and repository hygiene

| ID | Sev | Location | Finding |
|----|-----|----------|---------|
| B50 | U | `trapo-seq.py`, `a_data_mapping.py` | All helper paths are relative to the cwd; the tool only works from the repo root. |
| B51 | Q | repo | `__pycache__/*.pyc` committed (case-colliding `k_DR_logo`/`k_dr_logo`). |
| B52 | Q | `g_blast_annot*.py` | `Bio.Blast.Applications` imported (deprecated/removed in recent Biopython). |
| B53 | U | `test_script.sh` | Calls `python`, unquoted command expansion, writes some outputs to cwd. |
| B54 | U | README | Python ≥ 3.5 claimed (code needs ≥ 3.7); several CLI examples don't match the code. |
