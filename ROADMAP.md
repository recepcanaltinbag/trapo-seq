# trapo-seq roadmap

Goal: turn trapo-seq into a validated, installable tool and publish it.
Work happens on the `dev` branch; `main` stays at the last release until a
milestone is merged. Finding IDs (B01…) refer to [`docs/AUDIT.md`](docs/AUDIT.md).

Status legend: ✅ done · 🔄 in progress · ⬜ planned

GitHub issues: #9 (macOS FASTA) → B01 ✅ · #8 (temp dir) → B17/batch temp ✅ ·
#7 (multi-contig map_dist) → B45 ✅ · #3 (better DR finder) → B30–B38 ✅ ·
#4 (TIR logo) → M4 · #6 (cluster insertions) → M3/M4 · #5 (plasmid–genome
homology) → M1 · #1 (insertion-level quality) → M4. Fixed issues are closed by
the `dev` → `main` pull request.

---

## M0 — Foundations

| | Task | Acceptance criterion |
|---|------|----------------------|
| ✅ | Code audit of all modules | `docs/AUDIT.md`, 40+ findings with IDs |
| ✅ | Ground-truth read simulator (`simulate/`) | Error-free simulation reproduces IS + TSD exactly for 100 % of reads |
| ✅ | `environment.yml` (conda) with all external tools | `conda env create -f environment.yml` then `trapo-seq test` passes |
| ✅ | Remove committed `__pycache__`, update `.gitignore` (B51) | clean `git status` after a run |
| ✅ | Paths independent of the working directory (B50) | pipeline runs from any directory and on paths with spaces |

## M1 — Correctness (blocks everything below)

Priority order: results that are scientifically wrong first, then crashes.

| | Group | IDs |
|---|-------|-----|
| ✅ | Annotation: real E-value cutoff, best hit, one row per insertion, header IDs, off-by-one | B10–B14, B17, B18 |
| ✅ | Insertion finder: secondary alignments, multi-record reads, strand, `=`/`X` | B02–B05 |
| ✅ | Mapping: portable FASTA conversion, fail loudly, no stale outputs | B01, B06, B07 |
| ✅ | DR finder: flank logic, circular origin, read-ID parsing, thresholds | B30–B35 |
| ✅ | Stats/plots: counting per insertion, crashes, multi-contig | B36–B48 |
| ✅ | CLI: `--is_fasta` fallback, `--keep_temp`, `--gap`, output names | B15, B16, B20, B33, B47 |
| ✅ | Circular plasmid in mapping (inserts near the origin carry plasmid sequence) | B09 |
| ✅ | Inserts with plasmid flanks but no IS/genome hit are labelled *Contamination*; relabel as `Unannotated` and reserve *Contamination* for reads without plasmid flanks (seen on bundled test data: 34/41 inserts with a shared 4 bp DR) | — |
| ⬜ | **Chimera / artefact flags per insertion**: true insertions have adjacent plasmid flanks (overlap = TSD), same strand, both flanks present; ONT chimeras and ligation artefacts show non-adjacent flanks, strand flips, one-sided junctions or internal adapters. Add `flanks`/`flank_gap`/`strand_consistent` evidence to the insertion table, an `artefact_flag`, and report support from independent reads (both read strands, distinct read ends). Seen on test data: E. coli segment at 3243 is one-sided in all 3 reads (artefact), IS10 at 4839 has a 9 bp DR on reads from both strands (genuine). | — |
| ⬜ | Simulator: ONT chimeras (plasmid–plasmid, plasmid–genome, plasmid–foreign), internal adapters, barcode bleed-through, so the flags can be benchmarked | — |
| ⬜ | **Element boundary vs TSD ambiguity** when element termini resemble the target site (VD2 ISPs1: termini `TAGG…CCTA`, plasmid DR `CCTA`; curated 913 bp vs read consensus 904 bp; DR finder −2 with curated bounds). Resolve with per-site junction consensus across reads and multi-copy genomic comparison; report boundary confidence | — |
| ✅ | Plasmid–genome homologous regions cause false insertions: detect shared segments (plasmid vs genome BLAST), mask or flag them; add homologous segments to the simulator to test | #5 |

Acceptance: the original `test-data` and a simulated dataset run end to end
with no errors; each fix has a regression test.

## M2 — Validation (what makes it publishable)

| | Task | Acceptance criterion |
|---|------|----------------------|
| ✅ | `simulate/evaluate.py`: score pipeline output against `truth.tsv` | per-read and per-site precision, recall, F1; position error (bp); IS-name accuracy; TSD-length accuracy; contamination classification |
| ⬜ | Benchmark grid | read error 1–10 %; insertion rate 1–50 %; depth 100–10 000 reads; novel-insert fraction; fragmentation |
| ⬜ | E-value / identity threshold sweep | justified defaults (replaces E = 0, B10) |
| ⬜ | Comparison with a general long-read SV caller (e.g. Sniffles2) on the same simulated data | table + figure |
| ⬜ | Real-data case study | trapo-seq on the authors' trap-plasmid runs; findings vs. Sanger/PCR-confirmed insertions if available |

## M3 — Engineering quality

| | Task |
|---|------|
| ⬜ | `pyproject.toml`, `trapo-seq` console entry point, version from one place |
| ⬜ | `pytest` suite (unit tests + small end-to-end test on simulated data) |
| ⬜ | GitHub Actions CI (Linux + macOS) |
| ✅ | Merge single/batch duplicates; one BLAST call per sample instead of per insertion (B19, B38) |
| ⬜ | Speed at real-data scale (50–100k reads/barcode): profile every step; cluster identical inserts before BLAST (#6) |
| ⬜ | `logging` instead of prints and stdout redirection |
| ⬜ | Bioconda recipe; Zenodo DOI for releases |

## M4 — Features

| | Task |
|---|------|
| ⬜ | Insertion-site clustering across reads (unique events vs. read counts) |
| ⬜ | Per-sample summary report (HTML/TSV): reads, insertion frequency, IS spectrum, hotspots, TSDs |
| ⬜ | Condition comparison (e.g. stress vs. control) with statistics |
| ⬜ | Automatic extraction and clustering of novel (non-database) inserts |
| ✅ | **Novel element characterisation**: consensus sequence from all reads of a cluster; exact element boundaries; **TIR detection** (terminal inverted repeats at the element ends); DR check both in the plasmid and at the element's genomic copy; family hints from TIR/DR length and transposase presence |
| ⬜ | Event types beyond simple insertions: cointegrates (replicative transposition, two element copies), one-ended / IS26-type translocatable units, elements without TSD (IS91, IS110, IS200/605) |
| 🔄 | Simulator: add TIRs to novel elements (done; cointegrates pending) and cointegrate events so the new modules can be benchmarked |
| ⬜ | TIR sequence logos, like DR logos (#4) |
| ⬜ | Insertion-level base quality instead of whole-read average, to flag sequencing artefacts (#1) |
| ✅ | Realistic ONT error model in the simulator (indel-dominant, homopolymer-aware, per-read accuracy; R9.4.1 / R10.4.1 profiles) |

## M5 — Paper

| | Task |
|---|------|
| ⬜ | Choose venue (see below) and format |
| ⬜ | Figures: workflow, benchmark accuracy, SV-caller comparison, case study |
| ⬜ | Manuscript draft, data availability (SRA + Zenodo), software availability |

Candidate venues, considering Boğaziçi's Read & Publish agreements
(free or discounted APC — confirm with the library before submission):
OUP (*Bioinformatics*, *NAR Genomics and Bioinformatics*), Wiley
(*Molecular Ecology Resources*, *MicrobiologyOpen*), Royal Society
(*Royal Society Open Science*). Free fallbacks: *JOSS*, or a hybrid
journal through the non-OA route.
