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
- **`novel_elements` subcommand** (`src/l_novel_elements.py`): clusters
  inserts labelled *Genome* across samples by their genomic interval, defines
  element boundaries from read alignment ends refined by **terminal inverted
  repeats (TIRs)** and the **genomic direct repeat**, reports the longest ORF,
  plasmid insertion sites, genomic copy number, similarity to curated IS, and a
  verdict (`mobile_element_candidate` / `partial_evidence` / `genome_segment`).
  On 30 simulated datasets (60 synthetic elements): 59/60 recovered as one
  cluster each, boundaries exact for 53/59 and within 2 bp for all, TIR found
  for 58/59, genomic DR exact for 53/59; random genome intervals give 1.1 %
  `mobile_element_candidate`.
- **Artefact detection** (`src/artefacts.py`, see `docs/ARTEFACTS.md`): per-call
  evidence (`insertion_from_bam_artefacts.tsv`: flank presence, gap, strands,
  MAPQ, distance to the cut site), automatic cut-site / enzyme inference from
  read ends (`insertion_from_bam_cut_site.tsv`), and flags for non-adjacent
  flanks, strand flips, non-plasmid flanks, junctions at the cut site,
  overhang-sized DRs, one-sided calls, internal adapters, fold-backs,
  plasmid-self inserts, single-molecule support and barcode bleed-through
  (`sites_across_barcodes.tsv`). `best_alignment.tab` gains an `Artefact`
  column; `--artefact_action exclude` relabels flagged calls *Contamination*.
  Simulated chimeras / ligation chimeras / dimers / fold-backs: false calls
  310 → 0 (uniform) and 323 → 0 (ONT r9.4.1); 0 of 2162 true insertions
  flagged. New CLI options: `insert_finder(_batch) --plasmid --enzyme
  --cut_site --overhang --adapter_fasta`, `blast_annot(_batch)
  --artefact_action`, `--bleed_pct`. Simulator: `--chimera_rate`,
  `--ligation_chimera_rate`, `--overhang`, `--dimer_rate`, `--foldback_rate`.
- `scripts/build_element_db.py`: builds a non-redundant `--is_fasta` database
  from locally downloaded ISfinder, TnCentral/INTEGRALL and PICOTA files
  (repairs malformed records, merges identical/reverse-complement sequences,
  maps families to the `Family_Name` header format, selects element types).
- **Plasmid–genome shared segments** (issue #5): `map` detects them
  (`plasmid_genome_homology.tsv`, warning for segments ≥ the insertion
  threshold); the transposition rate in `08_mapping_stats.txt` no longer counts
  reads that touch the genome only through a shared segment (the bundled
  pMAT1 / E. coli K-12 pair shares 230 bp, which inflated the rate from ~0.27
  to 0.99; uncorrected values are kept as extra lines); insertion calls whose
  flanks lie in a shared segment and whose read continues collinearly into the
  genome are flagged (`insertion_from_bam_homology_flags.tsv`) and relabelled
  *Contamination* (`--homology_action keep` only flags). `best_alignment.tab`
  gains a trailing `Homology` column. Simulated 1.5–3 kb shared segments:
  all 144–173 homology-driven false calls removed, no true insertion flagged.
  Simulator `--homology_len` plants shared segments; `evaluate.py` reports a
  homology section.
- **`Unannotated` label**: inserts with plasmid on both sides but no IS or
  genome hit are no longer called *Contamination*. *Contamination* is kept for
  read-end clips that do not resemble any both-flanked insert of the sample.
  On the bundled test data 34/41 inserts move from Contamination to
  Unannotated. `.rcp` files gain a `#Unannotated` line; DR logos are drawn
  separately for Unannotated inserts.
- `novel_elements` also characterises **Unannotated inserts from the reads
  alone**: all-vs-all clustering, consensus (MAFFT + pileup polishing), TIR,
  longest ORF, plasmid DR from `dr_finder`. On simulations where the genome
  lacks the elements, every element is one cluster with a 99.9–100 % consensus
  and its TIR found (uniform, ONT r10.4.1, ONT r9.4.1).
- `evaluate.py --novel_absent` scores Unannotated as correct for novel inserts.
- Simulator `--error_model ont` with `--ont_profile r9.4.1|r10.4.1`:
  per-read accuracy from a Beta distribution, indel-dominant composition,
  homopolymer length errors, error-correlated qualities; all parameters
  overridable (presets are approximations). Runs in time proportional to the
  number of errors: 20 k reads in 21.5 s (r9.4.1) / 9.1 s (r10.4.1) vs ~180 s
  for the default uniform model, whose output is unchanged.
- Simulator `--is_divergence F` and `--is_truncate MIN MAX`
  (`--is_truncate_frac`): curated-IS clones carry a divergent and/or truncated
  copy of the element; defaults keep output byte-identical.
- Benchmark: `run_benchmark.py --reuse_work / --from_step / --snapshot_after`
  with a stamp that refuses mismatched reuse, and per-step CPU time;
  `run_grid.py` shares mapping + insert finding across annotation variants
  (identical metrics, ~35 % less work) and adds IS divergence/truncation sweeps
  (full preset: 549 cells).
- Simulator `--novel_mode element`: synthetic mobile elements with TIRs, an
  ORF and a genomic DR, inserted into a copy of the genome and then into the
  plasmid; truth in `elements.tsv`. Default output is unchanged (byte-identical).

### Changed
- **Bundled test data replaced by a simulated dataset with known truth**
  (E. coli K-12 host, pMAT1, 2 barcodes × 200 ONT r10.4.1-like reads,
  insertion rates 0.3 / 0.6; truth in `test-data/truth/`). The previous
  46 reads came from a different host than the bundled genome, so most
  inserts could not be annotated. `trapo-seq test` now also scores the run
  against the truth and fails if accuracy drops below set minimums.
- Mapping streams minimap2 output straight into samtools: the intermediate
  `01_aligned_to_plasmid.sam`, `05_aligned_to_other_genome.sam` and
  `06_aligned_to_other_genome.bam` are no longer written (none was used
  downstream); disk use per sample 2.1 GB → 1.3 GB on 20 k reads.
- `filter` copies 4-line FASTQ records byte for byte with a plain parser
  (17.7 s → 1.2 s on 20 k reads), accepts `.gz`, writes atomically, and falls
  back to Biopython for multi-line FASTQ.
- `map --threads` (default 4) for minimap2 and samtools: on 20 k simulated
  reads, mapping 138 s → 58 s and the whole pipeline 192 s → 106 s with 8
  threads, identical results.
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
- **DR (target-site duplication) finder rewritten** (B30–B38): correct flank
  pairing, strand-aware, plasmid coordinates throughout; circular plasmid
  handled (DRs across the origin); signed repeat length distinguishes DR,
  blunt insertion and target-site deletion; `--gap` implemented; read IDs
  with `_` parsed correctly; one BLAST run per sample. Optional `--is_fasta`
  uses IS lengths to sharpen boundaries. On simulated data: exact TSD length
  100 % (error-free) and 78.8 % (default errors; 98.1 % within ±1 bp).
- DR logos: columns no longer shifted by modal-length trimming; MAFFT only
  for mixed lengths; robust subprocess and figure handling (B36, B37).
- **Circular plasmid handled in mapping** (B09): reads are also mapped to the
  plasmid rotated by half its length (`--no-long-join`), and per read the
  mapping that aligns clearly more of the read is kept, converted back to
  original coordinates (split at the origin). Inserts near the origin no
  longer include hundreds of bp of plasmid. Simulated data (default / error-free
  / ONT r9.4.1): IS-name accuracy 97.4 → 100 / 100 → 100 / 93.8 → 98.8 %;
  spurious sites 4 → 0 / 1 → 1 / 21 → 10; TSD exact length 62.5 → 75.6 /
  87.0 → 96.0 / 42.9 → 52.3 %. `map --linear` restores the old mapping.
  `03_sorted_mapped_to_plasmid.bam` may contain more supplementary records;
  `ref_pos` is always in [0, plasmid length).
- `src/j_dr_finder_multithread.py` (which had no threads) replaced
  `src/j_dr_finder.py` (dead duplicate) (B38).

## [0.1.0]
- Initial release.
