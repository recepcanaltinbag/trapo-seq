# Changelog

All notable changes to trapo-seq. Finding IDs (B01…) refer to
[`docs/AUDIT.md`](docs/AUDIT.md).

## [Unreleased] — `dev` branch

### Added
- `simulate/simulate_traposome.py`: ground-truth simulator for traposome
  long reads (IS and novel genome-segment insertions with target-site
  duplications, skewed clone abundances, restriction-site linearisation,
  fragmentation, strand, empty plasmids, deletions, genomic contamination,
  substitution/indel errors). Writes FASTQ per barcode plus `truth.tsv`,
  `clones.tsv` and `simulation.json`. Standard library only.
- `docs/AUDIT.md`: full code audit of v0.1.0.
- `ROADMAP.md`, `CHANGELOG.md`, `environment.yml`.
- `blast_annot*`: `--evalue` (default 1e-10), `--blast_task`
  (`megablast` | `dc-megablast` | `blastn`), `--keep_temp`.
- Test suite (`pytest`): simulator, CLI smoke tests, annotation unit tests.

### Changed
- **Annotation output now has exactly one row per insertion** (B12). Each
  insert is classified once: best curated-IS subject if it covers ≥
  `--threshold` % of the insert (`IS_DB`), otherwise best genome subject
  (`Genome`), otherwise `Contamination`. Inserts made of two different IS
  that each cover less than the threshold are no longer reported as `IS_DB`.
- One `blastn` run per sample and database (tabular output) instead of one
  run per insertion (B19); much faster on large samples.
- BLAST databases are built in the temp directory, named after the FASTA
  and rebuilt when the FASTA changes (B17).
- `--no_temp` replaced by `--keep_temp`; temp files are deleted by default (B16).
- Batch annotation temp directory defaults to `<input_dir>/temp`.
- Single-mode annotation writes `best_alignment.tab` (was
  `best_alignments.tab`, which downstream steps ignored) (B20).
- Version reported by `--version` comes from `src/__init__.py`.

### Fixed
- **BLAST hits were kept only when E-value == 0** (≈ < 1e-180), so genuine
  inserts from lower-quality reads were labelled *Contamination* (B10).
- Genome classification used the weakest HSP instead of the best (B11).
- Subject IDs were the last word of the FASTA header (B13).
- Coverage / "Explained" off-by-one on both strands (B14).
- `--is_fasta` fallback to the placeholder database never worked (B15).
- External commands broke on paths with spaces; failures were ignored (B18).
- Deprecated `Bio.Blast.Applications` removed (B52).
- `read_histogram -v` did not draw the histogram (B47).
- The CLI only worked when started from the repository root (B50).
- `__pycache__` files were committed (B51).
- `is_stat` crashed on IS names with zero or several underscores and on
  samples without IS hits (B40); it now counts insertions (unique Query ID),
  not BLAST rows, and records the IS_DB / Genome / Contamination breakdown
  and the percentage denominator as `#` lines in the `.rcp` file (B41).
- `heatmap` crashed when all IS fell on one side of the 15 % split (B42).
- `in_del_plot` overwrote soft-clip insertion counts with CIGAR insertion
  counts instead of adding them (B43); `=`, `X`, `N` CIGAR ops handled.
- `kde_mobile` used hard-clip-blind, strand-blind read coordinates and
  dropped genome segments of split reads (B44); no longer crashes when no read
  passes the threshold, and the scatter shows only filtered reads (B46).
- `map_dist` summed coverage of all contigs into the first contig (B45);
  mapping statistics are now written (`--fastq` optional).
- `read_histogram` crashed on uniform/empty input and wrote `x.pdf.pdf` (B48).
- **Mapping on macOS produced an empty read FASTA** (GNU-only `sed`), which
  later crashed annotation (B01); FASTA now written by `samtools fasta`.
- Mapping failures exited with status 0 and left partial files that later
  runs silently reused (B06); steps now write to temp files, fail loudly, and
  rebuild all downstream outputs when an upstream step is rebuilt.
- Division by zero in the transposition rate when no reads map (B07).
- **Insertion finder**: CIGAR insertions from earlier alignment records of a
  read were overwritten (B02); secondary alignments corrupted quality,
  coverage and positions (B03); soft-clip insertion strand came from the
  wrong alignment (B04); `=`/`X` CIGAR ops ignored (B05); insertions near the
  plasmid origin were placed at position 0 / plasmid end (B08). On simulated
  data, insertion reads placed within 20 bp rose from 66 % to 92 %, with no
  calls on empty-plasmid reads.

## [0.1.0]
- Initial release.
