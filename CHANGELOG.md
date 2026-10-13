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

## [0.1.0]
- Initial release.
