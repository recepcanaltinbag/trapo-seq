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
- `ROADMAP.md`, `CHANGELOG.md`.

## [0.1.0]
- Initial release.
