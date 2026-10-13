# What can look like an insertion, a DR or an element boundary

A catalogue of every known source of false or distorted calls in long-read
traposome data, the signature that separates it from a genuine transposition,
and how trapo-seq handles it. Status: ✅ handled · 🟡 partly · ⬜ planned.

A **genuine simple insertion** has all of the following:

1. plasmid sequence on **both** sides of the insert, on the **same strand**;
2. the two plasmid flanks are **adjacent** on the plasmid, overlapping by the
   target-site duplication (TSD/DR) of the element family (0–14 bp);
3. the insert ends at the element's termini (often terminal inverted repeats, TIRs);
4. the same site, DR and insert are seen in **independent molecules**
   (reads from both strands, different read start/end positions).

Every artefact below violates at least one of these.

## A. Sample and library

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| A1 | Host genomic DNA contamination | No plasmid flank, or only a short one | ✅ not called (needs plasmid alignment); labelled *Contamination* if one-sided |
| A2 | Plasmid–genome shared segments (issue #5) | Insert flank lies inside a shared segment; insert is genome DNA continuing the shared segment | 🟡 in progress (#5) |
| A3 | Barcode bleed-through / demultiplexing errors | Same site + insert as in another barcode, at a much lower read count; no barcode-specific events | ⬜ cross-barcode site table with count ratios |
| A4 | Plasmid dimers and concatemers | Insert maps to the plasmid itself; read longer than twice the plasmid; junctions at the restriction site | ⬜ label inserts that map to the plasmid (`Plasmid_self`) |
| A5 | Cointegrates (replicative transposition, real biology) | Donor replicon between two direct copies of the element, each with a TSD | ⬜ event type planned (ROADMAP M4) |
| A6 | Elements from native plasmids/phages of the host | Genuine insert, absent from the genome FASTA | ✅ *Unannotated* + read consensus (`novel_elements`) |
| A7 | Passage through a cloning host (e.g. E. coli) | Elements of the cloning host (IS1, IS10, IS2, IS5, …) in a non-enteric host's data | 🟡 reported as IS_DB/Unannotated; interpretation is the user's |

## B. Restriction digest / linearisation

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| B1 | Insert close to the cut site | One flank very short or missing; insert reaches a read end | 🟡 one-sided calls are flagged by the DR finder (`left/right_flank_only`) |
| B2 | **Sticky-end ligation chimeras** | Junction exactly at the cut site; apparent "DR" equals the enzyme overhang (e.g. 4 bp `AATT` for EcoRI); foreign or plasmid DNA between | ⬜ flag calls with `ref_pos` within ±overhang of the cut site whose DR equals the overhang (needs the enzyme/cut site as input or inferred from read ends) |
| B3 | Restriction site inside the insert | Read ends inside the insert; plasmid on one side only; two reads cover the two halves | ⬜ pair one-sided calls at the same site with complementary halves |
| B4 | Additional or star-activity sites | Read ends piling up at secondary positions; shorter reads | 🟡 visible in read-length histogram / map_dist |
| B5 | Incomplete digestion, nicked circles | Reads start anywhere; reads longer than the plasmid that wrap around | ✅ circular-aware mapping (B09) |

The cut site can be inferred from the data: most read ends of unfragmented
plasmid reads pile up at it. Planned: infer it automatically and report it.

## C. Sequencing (Oxford Nanopore)

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| C1 | Chimeric reads (two molecules read as one, ligation of adapters/molecules) | Plasmid flanks **not adjacent** (gap or large overlap), or strand flip, or unrelated sequences joined | ⬜ flank adjacency/strand check per call (`artefact_flag`) |
| C2 | Internal adapter / barcode sequence | Adapter motif inside the read at the junction | ⬜ scan inserts and junctions for adapter/barcode motifs |
| C3 | Fold-back (inverted) reads | Read = sequence + its reverse complement; "insert" is plasmid in inverted orientation | ⬜ detect self-reverse-complement reads |
| C4 | Untrimmed low-quality or adapter ends | Soft clips at read ends | 🟡 read-end clips need resemblance to a both-flanked insert to be kept (*Unannotated* rule), else *Contamination* |
| C5 | Homopolymer/indel errors at junctions | DR length ±1, boundary jitter | 🟡 per-read DR; per-site consensus planned |
| C6 | Basecaller read splitting/joining | Read IDs containing parent IDs; truncated/merged reads | 🟡 read IDs with `_` are handled; joined reads behave like C1 |

## D. Analysis and biology that blur boundaries

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| D1 | Repeats inside the plasmid | Multi-mapping flanks; ambiguous `ref_pos` | ⬜ report MAPQ of flanks; flag ambiguous sites |
| D2 | **Element termini resemble the target site** (e.g. VD2 ISPs1 `TAGG…CCTA`, DR `CCTA`) | Boundary and DR length trade off; curated vs consensus lengths differ by a few bp | ⬜ per-site junction consensus + multi-copy genomic comparison, report boundary confidence |
| D3 | Chance flank extension | DR 1–2 bp longer than the family's typical TSD | 🟡 IS bounds (`--is_fasta`) clip flanks; novel elements +1 bp typical |
| D4 | Linear reference for a circular plasmid | Inserts near the origin carry plasmid sequence | ✅ fixed (B09) |
| D5 | Nested / composite elements | Inner element hit covers part of the insert | 🟡 one row per insert; inner elements not reported separately |
| D6 | Deletions next to insertions (adjacent deletions, excision footprints, IS26/Tn3 resolution) | Negative overlap (`target_site_deletion`) | ✅ reported as signed repeat length |
| D7 | Elements that leave no TSD (IS91, IS110, IS200/IS605 families, integron cassette recombination) | Blunt junction, length 0 | ✅ reported as `no_duplication`, not as an error |
| D8 | Insertion-site hotspots | Many independent events at one site | ✅ distinguish from clonal amplification by read start/end diversity (planned in C1 work) |

## Benchmarking

Each ⬜/🟡 item should get a simulator option so its false-positive rate and
the flag's sensitivity can be measured (ROADMAP: chimera simulation).
