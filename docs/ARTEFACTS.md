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
| A2 | Plasmid–genome shared segments (issue #5) | Insert flank lies inside a shared segment; insert is genome DNA continuing the shared segment | ✅ detected in `map`; transposition rate corrected; calls flagged in `Homology` column (#5) |
| A3 | Barcode bleed-through / demultiplexing errors | Same site + insert as in another barcode, at a much lower read count; no barcode-specific events | 🟡 `sites_across_barcodes.tsv` (reads per barcode per site); calls at a site with < 5 % of the reads it has in another barcode get the caution flag `possible_bleed_through`. Not simulated yet |
| A4 | Plasmid dimers and concatemers | Insert maps to the plasmid itself; read longer than twice the plasmid; junctions at the restriction site | ✅ `plasmid_self` (insert interior ≥ 50 % covered by plasmid 13-mers); junctions at the cut site → `cut_site_junction`. Simulated head-to-tail/inverted dimers align as two plasmid copies and give no call |
| A5 | Cointegrates (replicative transposition, real biology) | Donor replicon between two direct copies of the element, each with a TSD | ⬜ event type planned (ROADMAP M4) |
| A6 | Elements from native plasmids/phages of the host | Genuine insert, absent from the genome FASTA | ✅ *Unannotated* + read consensus (`novel_elements`) |
| A7 | Passage through a cloning host (e.g. E. coli) | Elements of the cloning host (IS1, IS10, IS2, IS5, …) in a non-enteric host's data | 🟡 reported as IS_DB/Unannotated; interpretation is the user's |

## B. Restriction digest / linearisation

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| B1 | Insert close to the cut site | One flank very short or missing; insert reaches a read end | ✅ `one_sided_cut_site`; kept as `one_sided_cut_site_is` (caution) when the clip is ≥ 90 % explained by a curated IS (element + short stub) |
| B2 | **Sticky-end ligation chimeras** | Junction exactly at the cut site; apparent "DR" equals the enzyme overhang (e.g. 4 bp `AATT` for EcoRI); foreign or plasmid DNA between | ✅ cut site inferred from read ends (or `--enzyme` / `--cut_site --overhang`); two-sided calls within \|end offset\| + 20 bp of it: `cut_site_junction`, plus `overhang_dr` when the flank overlap equals the end signature; one-sided: `one_sided_cut_site`. Simulated ligation chimeras: 100 % flagged |
| B3 | Restriction site inside the insert | Read ends inside the insert; plasmid on one side only; two reads cover the two halves | 🟡 caution flag `cut_in_insert` (recognition site of the given/guessed enzyme at the read end of a one-sided call; seen on VD2: ISPen2 carries an XbaI site). Pairing of complementary halves planned |
| B4 | Additional or star-activity sites | Read ends piling up at secondary positions; shorter reads | 🟡 visible in read-length histogram / map_dist; `*_cut_site.tsv` gives the support of the dominant site; an enzyme with several sites is warned about |
| B5 | Incomplete digestion, nicked circles | Reads start anywhere; reads longer than the plasmid that wrap around | ✅ circular-aware mapping (B09) |
| B6 | Untrimmed adapters/barcodes at read ends | Alignment extends a few bp into the adapter where its bases happen to continue the plasmid; molecule ends shifted | ✅ tolerated (enzyme guess ±3 bp; cut-site window +20 bp). Seen on VD2: XbaI start 3239 / filled-in end 3243 observed as 3237 / 3246 |

The cut site is inferred from the data (src/artefacts.py): the molecule start
and end with most unfragmented plasmid reads (aligned ≥ 80 % of the plasmid,
≤ 200 bp unaligned at each read end) within ±3 bp, reported when ≥ 20 % of
those reads (and ≥ 5) support it. `end − start` is the end offset: +k for a
5′ overhang filled in by end repair, −k for a 3′ overhang removed by it, 0 for
a blunt cut. With the plasmid sequence, a table enzyme whose site explains both
ends (±3 bp) is reported. Simulation (12 runs, uniform and ONT r9.4.1, blunt
and +4 overhang, 10/60 % fragmented molecules): start and offset exact in all
8 runs with intact molecules; no site reported when every molecule was
fragmented. Written to `<sample>/insertion_from_bam_cut_site.tsv`.

## C. Sequencing (Oxford Nanopore)

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| C1 | Chimeric reads (two molecules read as one, ligation of adapters/molecules) | Plasmid flanks **not adjacent** (gap or large overlap), or strand flip, or unrelated sequences joined | ✅ `nonadjacent_flanks` (flank gap outside [−100, +40] bp), `strand_flip`, `unaligned_flank` (flank identity < 0.65, or a < 50 bp "flank" between two unaligned stretches: aligner bridged the junction), `one_sided_unsupported` (one-sided Genome clip with no two-sided Genome call at that genome locus) |
| C2 | Internal adapter / barcode sequence | Adapter motif inside the read at the junction | ✅ `internal_adapter`: insert ± 100 bp scanned (not the outer 150 bp of the read) for the ONT ligation Y-adapter top strand (≤ 18 % edits); other kits' adapters/barcodes via `--adapter_fasta` |
| C3 | Fold-back (inverted) reads | Read = sequence + its reverse complement; "insert" is plasmid in inverted orientation | ✅ `foldback` (the 500 bp before the insert are the reverse complement of the 500 bp after it, ≥ 10 shared 11-mers on one diagonal); `foldback_read` (caution: the read covers one plasmid region on both strands, its calls are one molecule twice) |
| C4 | Untrimmed low-quality or adapter ends | Soft clips at read ends | 🟡 read-end clips need resemblance to a both-flanked insert to be kept (*Unannotated* rule), else *Contamination* |
| C5 | Homopolymer/indel errors at junctions | DR length ±1, boundary jitter | 🟡 per-read DR; per-site consensus planned |
| C6 | Basecaller read splitting/joining | Read IDs containing parent IDs; truncated/merged reads | 🟡 read IDs with `_` are handled; joined reads behave like C1 (flagged); split reads (one continues where the other stops) count as one molecule in the site support |
| C7 | Duplex / complement reads | Template and complement of one molecule: opposite strands, identical ends | ✅ counted as one molecule (`min_molecules`, `single_molecule_support`) |

## D. Analysis and biology that blur boundaries

| # | Source | Signature | Status |
|---|--------|-----------|--------|
| D1 | Repeats inside the plasmid | Multi-mapping flanks; ambiguous `ref_pos` | ✅ flank MAPQ reported (`left/right_mapq`, `min_flank_mapq`); caution flag `low_mapq` (< 10) |
| D2 | **Element termini resemble the target site** (e.g. VD2 ISPs1 `TAGG…CCTA`, DR `CCTA`) | Boundary and DR length trade off; curated vs consensus lengths differ by a few bp | ⬜ per-site junction consensus + multi-copy genomic comparison, report boundary confidence |
| D3 | Chance flank extension | DR 1–2 bp longer than the family's typical TSD | 🟡 IS bounds (`--is_fasta`) clip flanks; novel elements +1 bp typical |
| D4 | Linear reference for a circular plasmid | Inserts near the origin carry plasmid sequence | ✅ fixed (B09) |
| D5 | Nested / composite elements | Inner element hit covers part of the insert | 🟡 one row per insert; inner elements not reported separately |
| D6 | Deletions next to insertions (adjacent deletions, excision footprints, IS26/Tn3 resolution) | Negative overlap (`target_site_deletion`) | ✅ reported as signed repeat length |
| D7 | Elements that leave no TSD (IS91, IS110, IS200/IS605 families, integron cassette recombination) | Blunt junction, length 0 | ✅ reported as `no_duplication`, not as an error |
| D8 | Insertion-site hotspots | Many independent events at one site | ✅ per-site table `<sample>/best_alignment_sites.tsv`: reads, reads per strand, distinct read start/end pairs, minimum number of molecules; caution flag `single_molecule_support` |

## Flags in the outputs

`insert_finder` writes the evidence per call to
`<sample>/insertion_from_bam_artefacts.tsv` (flank coordinates, strands,
MAPQ, flank gap, distance to the cut site, adapter hits, fold-back score,
plasmid self-match, flags). `blast_annot` copies the flags into the last
column `Artefact` of `best_alignment.tab`; with `--artefact_action exclude`
calls with an **artefact** flag are relabelled *Contamination*
(`excluded:<Note>:<flags>`). **Caution** flags (`one_sided`,
`one_sided_cut_site_is`, `cut_in_insert`, `foldback_read`, `low_mapq`,
`single_molecule_support`, `possible_bleed_through`) never change a Note.
`is_stat` reports `#Artefact_<flag>` counts, `dr_finder` writes flagged calls
under Subject ID `Artefact` (separate DR logo), `novel_elements` does not use
them. Thresholds and their calibration: docstring of `src/artefacts.py`.

## Benchmarking

`simulate/simulate_traposome.py --chimera_rate --ligation_chimera_rate
--overhang --dimer_rate --foldback_rate` adds artefact reads;
`simulate/evaluate.py` reports, per category, the reads called before and
after the flags (section `artefacts`). Fresh seeds, 2 × 2000 reads + 10 %
chimeras, 5 % ligation chimeras (+4 overhang), 3 % dimers, 5 % fold-backs:

| | uniform 2 % sub | ONT r9.4.1 |
|---|---|---|
| false insertion calls (annotation) before → after flags | 310 → 0 | 323 → 0 |
| true insertion reads with an artefact flag | 0 / 1102 | 0 / 1060 |
| precision / recall, `flag` → `exclude` | 0.78 / 0.944 → 1.00 / 0.944 | 0.77 / 0.937 → 1.00 / 0.937 |

Barcode bleed-through (A3) and 3′-overhang enzymes are not simulated yet.
