"""Characterise novel (non-database) inserts as candidate mobile elements.

Two kinds of insert are characterised (column ``source``):

* ``genome`` - inserts labelled ``Genome`` (steps 1-8 below);
* ``reads``  - inserts labelled ``Unannotated``: flanked by plasmid but found
  neither in the curated IS database nor in the genome (section UNANNOTATED
  INSERTS).

Inserts that best_alignment.tab labels ``Genome`` match the host genome but
not the curated IS database. If the same genomic interval jumped into the trap
plasmid in several reads (or samples), it may be an uncatalogued mobile
element. For every such interval this module reports:

* the element boundaries on the host genome (from the read alignments, then
  refined with terminal inverted repeats when they are found),
* terminal inverted repeats (TIRs) at the two ends,
* the direct repeat (DR, target-site duplication) flanking the genomic copy,
* the plasmid insertion sites and the samples it was seen in,
* the longest ORF inside the element (ATG to stop, standard code; a proxy for
  coding capacity such as a transposase - no annotation is implied),
* a verdict (see VERDICT RULE).

Algorithm
---------
1. Inserts: rows of ``<sample>/best_alignment.tab`` with Note == Genome. The
   insert is ``reads[read_id][start:end]`` of ``<sample>/<sample>.fasta``,
   coordinates from the Query ID ``{read_id}_{start}_{end}`` (0-based
   half-open, original read orientation, see e_insert_finder_from_bam).
2. Inserts are mapped to the genome with ``minimap2 -c -x map-ont``. The
   inserts are raw nanopore read segments (several % errors), which the
   assembly presets (asm5/asm10: < 1-5 % divergence) are not tuned for;
   ``-c`` computes base-level alignments so the alignment ends are exact
   rather than chaining estimates. Per insert the best primary alignment
   (most matching bases) is used.
3. Clustering: intervals on the same contig with reciprocal overlap >=
   ``min_overlap`` (default 0.5) form one cluster (pooled over samples; longest
   inserts seed the clusters). Clusters lying >= 90 % inside a larger cluster
   are absorbed into it (fragmented reads give partial inserts). Clusters with
   fewer than ``min_reads`` reads are dropped.
4. Read-based boundaries: the start (end) is the mode of the alignment
   starts (ends) inside the densest ``end_window`` (10) bp window - full
   length inserts agree, partial ones scatter. Spread = median absolute
   deviation of all reads from the estimate. In simulations these estimates
   are exact for most elements: the plasmid flanks differ from the genomic
   flanks, so alignments stop at the element ends. They overshoot by 1 bp
   when a plasmid flank base happens to match the genome - for every read of
   that insertion site.
5. TIR (find_tir): the genome +-``flank_search`` (60) bp around the estimated
   start is locally aligned (Biopython PairwiseAligner; match 2, mismatch -3,
   gap open -5, gap extend -2, close to blastn's nucleotide scoring, which
   tolerates ~10-20 % mismatches) with the reverse complement of the genome
   around the estimated end. The alignment is then trimmed: the maximal
   stretch under a strict score (+1 match, ``trim_mismatch`` -2 per mismatch,
   -3 per gap column; positive only above ~67 % identity) is kept, and both
   ends must start with ``edge_run`` (2) consecutive matches. Without this,
   chance pairing of flank bases extended TIRs outward by up to ~20 bp in
   tests. A TIR pair is accepted when it is >= ``min_tir_len`` (10) bp, >=
   ``min_tir_identity`` (0.75) identical and the local alignment scores >=
   ``min_tir_score`` (22). Calibration (1000 random E. coli K-12 intervals of
   1.2-3 kb): 2.1 % report a TIR at score >= 22 (>= 24: 0.7 %, >= 28: 0 %).
   Synthetic TIRs with 10 % mismatches are found in ~65 % (15 bp) and 100 %
   (>= 20 bp) of cases at 22; at 24 only 30 % of 15 bp TIRs would be found.
   All four TIRs of the curated test IS elements score >= 22 (one exactly 22).
6. Final boundaries and DR (find_dr): the genomic DR is the longest k in
   2..15 with genome[start-k:start] == genome[end:end+k] (one mismatch allowed
   for k >= 8; one mismatch in a shorter repeat is too likely by chance).
   The read estimate is kept unless
   a) a TIR was found: its outer edges, shifted symmetrically by up to
      ``dr_shift`` (3) bp, are scored by their exact DR; they replace the read
      estimate only if their DR is >= ``tir_override_gain`` (3) bp longer
      (a longer-by-2 DR arises by chance 1 in 16: the true DR plus one random
      pair on each side), or
   b) no TIR agrees with the estimate and other ends within 3 bp are each
      supported by >= 2 reads (and >= 20 % of the modal count): the read-
      supported combination with the longest exact DR wins (ties: more reads,
      then the estimate).
   ``boundary_source`` = reads | reads+TIR (TIR edges equal the read ends) |
   reads+DR (case b) | TIR (case a). The TIR is then re-searched within
   +-``tir_anchor`` (3) bp of the final ends so its length is measured from
   them (``tir_*_offset`` = TIR outer edge distance from the element end).
   TIR-only accuracy (estimates +-5 bp off, 600 synthetic elements with 15-40
   bp TIRs, 10 % mismatches outside the outer 2 bp): boundaries exact in 81 %
   and within 3 bp in 91 % of found TIRs, DR exact in 83 %.
7. Longest ORF on both strands of the element (ATG start, standard code,
   complete ORFs with stop codon only).
8. Context columns: ``genome_copies`` (genome alignments covering >= 90 % of
   the element, minimap2 asm10, own locus included), ``similar_clusters``
   (other clusters whose elements align over >= 80 % of both: reads of a
   multi-copy element are split between its copies), ``known_hit`` (best hit
   in ``is_fasta`` when given: name:identity%:element coverage%).

VERDICT RULE
------------
* ``mobile_element_candidate``: TIR found AND (genomic DR >= ``min_dr_len``
  (3 bp) OR longest ORF >= ``min_orf_aa`` (200 aa)).
* ``partial_evidence``: TIR found without DR/ORF support, or no TIR but both
  a DR >= ``min_dr_len`` and a long ORF.
  min_dr_len = 3: a 2 bp DR occurs by chance at ~1 in 16 loci, 3 bp at ~1 in
  64. min_orf_aa = 200: the longest ORFs of the curated test IS elements
  are 292-402 aa (Tn3_ISEc63: 1002 aa), while uniformly random 2 kb
  sequence gave a median longest ORF of 81 aa (max 169 in 200 trials).
* ``genome_segment``: anything else. A gene-dense genome segment often
  contains a long ORF (71 % of random 1.2-3 kb E. coli intervals have one
  >= 200 aa), so an ORF alone is not evidence. On 1000 such random intervals
  the verdicts were 96.9 % genome_segment, 2.0 % partial_evidence and 1.1 %
  mobile_element_candidate.

UNANNOTATED INSERTS (source = reads)
------------------------------------
No genome copy exists, so the element is reconstructed from the reads.
R1. Clustering: all-vs-all ``minimap2 -x ava-ont`` of the insert sequences;
    inserts aligned over >= ``min_overlap`` of both are linked and clusters
    are the connected components (see cluster_by_similarity: a chimeric
    insert cannot join two elements); a smaller component with a member
    covered >= 90 % by a member of a larger one is absorbed as partial.
    Clusters with < ``min_reads`` members are dropped.
R2. Consensus: the medoid (member with most matching bases to the other
    members) and up to ``max_consensus_reads`` - 1 members closest to it,
    oriented on the medoid, are aligned with MAFFT (--auto); a column-wise
    majority vote (terminal gaps of a read = no vote; a column needs votes
    from >= half of the reads; gap columns dropped) gives the consensus.
    racon/medaka are not used (no extra dependency). ``read_identity`` =
    median identity of the members to the consensus (minimap2 -c).
R3. TIR (find_tir) within ``consensus_flank`` (20) bp of the consensus
    ends; when found its outer edges trim the consensus to the element
    (boundary_source ``consensus+TIR``, else ``consensus``). The consensus
    may keep one copy of the plasmid TSD at an end when no TIR is found.
    Longest ORF as in step 7. Coordinates (start/end/orf_*) are on the
    consensus; contig, est_*, *_spread and the genomic DR columns are NA.
R4. Verdict as in VERDICT RULE, but DR evidence is the plasmid target-site
    duplication: the modal DR-finder repeat length of the cluster's inserts
    (plasmid_dr_len) >= ``min_dr_len`` in >= 50 % of its inserts with a
    DR-finder result; ``no_element_evidence`` replaces ``genome_segment``.
Plasmid DR columns (plasmid_dr_len, plasmid_dr_seq, plasmid_dr_support,
plasmid_dr_counts 'len:n;...') are filled for both sources from
``<data_dir>/insertions/<sample>_insertions.csv`` (run dr_finder first; NA
otherwise); for genome clusters they are informational only.

Outputs (``out_dir``, default ``<data_dir>/novel_elements/``):
  novel_elements.tsv    one row per cluster (columns: TSV_COLUMNS); genome
                        clusters ``novel_NNN`` first, then read clusters
                        ``unannotated_NNN``
  novel_elements.fasta  element sequences: from the genome (refined
                        boundaries) or the read consensus (TIR-trimmed)
  novel_elements.log    log of the run
All genome coordinates are 0-based half-open.

Requires minimap2 and (for read clusters) mafft in PATH, and Biopython.
"""

import csv
import logging
import math
import os
import shutil
import statistics
import subprocess
import tempfile
from collections import Counter

from Bio import SeqIO
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}

ANNOT_NAME = "best_alignment.tab"
INSERT_TAB_NAME = "insertion_from_bam.tab"

TSV_COLUMNS = [
    "cluster_id", "source", "contig", "start", "end", "length", "boundary_source",
    "n_reads", "n_partial_reads", "n_samples", "samples", "insert_strands",
    "est_start", "est_end", "start_spread", "end_spread",
    "consensus_reads", "consensus_length", "read_identity",
    "tir_found", "tir_left_len", "tir_right_len", "tir_identity", "tir_score",
    "tir_left_offset", "tir_right_offset", "tir_left_seq", "tir_right_seq",
    "dr_len", "dr_mismatches", "dr_left_seq", "dr_right_seq",
    "plasmid_dr_len", "plasmid_dr_seq", "plasmid_dr_support", "plasmid_dr_counts",
    "longest_orf_aa", "orf_start", "orf_end", "orf_strand",
    "genome_copies", "similar_clusters", "known_hit",
    "n_plasmid_sites", "plasmid_sites", "verdict",
]
NA = "NA"   # not applicable (e.g. genome coordinates / genomic DR of a read consensus)

# Element-characterisation defaults (see the module docstring for the reasons).
CHARACTERISE_DEFAULTS = dict(flank_search=60, max_tir=50, min_tir_len=10, min_tir_identity=0.75,
                             min_tir_score=22, edge_run=2, trim_mismatch=-2, dr_shift=3,
                             tir_override_gain=3, tir_anchor=3, min_dr_len=3, min_orf_aa=200)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


#------------------------------------------------------------------ TIR -----

def _aligner():
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score = 2
    al.mismatch_score = -3
    al.open_gap_score = -5
    al.extend_gap_score = -2
    return al


_ALIGNER = None


def find_tir(seq, start=0, end=None, flank=60, max_tir=50, min_len=10,
             min_identity=0.75, min_score=22, edge_run=2, trim_mismatch=-2):
    """Find a terminal inverted repeat pair around the ends of seq[start:end].

    The left end region ``seq[start-flank : start+flank+max_tir]`` is locally
    aligned with the reverse complement of the right end region
    ``seq[end-flank-max_tir : end+flank]`` (both clipped to the sequence and to
    the midpoint of the element). Use ``seq`` = a whole contig with estimated
    boundaries, or an element sequence with start=0, end=len(seq).

    A local alignment readily runs a few bp past the true TIR ends through
    chance matches (and across single mismatches). Both ends of the alignment
    are therefore trimmed back until they start with ``edge_run`` consecutive
    matching pairs.

    Returns None or a dict with
      start, end      outer edges of the TIR pair = element boundaries
                      (0-based half-open on seq)
      left, right     (start, end) of the left / right TIR on seq
      left_len, right_len, identity (matches / alignment columns),
      score (local alignment score before edge trimming), matches, columns,
      left_seq, right_seq  (right_seq as on the + strand of seq)
    """
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = _aligner()
    n = len(seq)
    if end is None:
        end = n
    mid = (start + end) // 2
    lw0 = max(0, start - flank)
    lw1 = min(n, start + flank + max_tir, mid)
    rw0 = max(0, end - flank - max_tir, mid)
    rw1 = min(n, end + flank)
    if lw1 - lw0 < min_len or rw1 - rw0 < min_len:
        return None
    left = seq[lw0:lw1].upper()
    right_rc = revcomp(seq[rw0:rw1].upper())
    try:
        aln = next(iter(_ALIGNER.align(left, right_rc)))
    except StopIteration:
        return None
    score = float(aln.score)
    if score < min_score:
        return None
    # aligned pairs (i on left window, j on reverse-complemented right window)
    pairs = [(int(a0) + k, int(b0) + k)
             for (a0, a1), (b0, _b1) in zip(*aln.aligned) for k in range(int(a1) - int(a0))]
    if not pairs:
        return None
    is_match = [left[i] == right_rc[j] for i, j in pairs]

    def run_at(x, step):
        """edge_run consecutive, gap-free matching pairs starting at x."""
        for r in range(edge_run):
            y = x + step * r
            if y < 0 or y >= len(pairs) or not is_match[y]:
                return False
            if r and (abs(pairs[y][0] - pairs[y - step][0]) != 1 or abs(pairs[y][1] - pairs[y - step][1]) != 1):
                return False
        return True

    # 1) keep the maximal-scoring stretch under a strict score (+1 match,
    #    trim_mismatch (-2) per mismatch, -3 per gap column): positive only
    #    above ~67 % identity, so low-identity or gapped chance extensions
    #    into the flanks are dropped.
    best_sum, best_rng, cur, cur_start = 0, None, 0, 0
    for x, ((i, j), m) in enumerate(zip(pairs, is_match)):
        gap = 0 if x == 0 else (i - pairs[x - 1][0] - 1) + (j - pairs[x - 1][1] - 1)
        val = (1 if m else trim_mismatch) - 3 * gap
        if cur <= 0:
            cur, cur_start = val, x
        else:
            cur += val
        if cur > best_sum:
            best_sum, best_rng = cur, (cur_start, x)
    if best_rng is None:
        return None
    lo, hi = best_rng
    # 2) both ends must start with edge_run gap-free matches
    f = next((x for x in range(lo, hi + 1) if run_at(x, 1)), None)
    b = next((x for x in range(hi, lo - 1, -1) if run_at(x, -1)), None)
    if f is None or b is None or b < f:
        return None
    (i0, j0), (ib, jb) = pairs[f], pairs[b]
    i1, j1 = ib + 1, jb + 1
    npairs = b - f + 1
    matches = sum(is_match[f:b + 1])
    # alignment columns = aligned pairs + gap columns in either sequence
    columns = (i1 - i0) + (j1 - j0) - npairs
    identity = matches / columns if columns else 0.0
    left_len, right_len = i1 - i0, j1 - j0
    tir_start = lw0 + i0
    tir_end = rw1 - j0
    ok = (min(left_len, right_len) >= min_len and identity >= min_identity
          and abs(tir_start - start) <= flank and abs(tir_end - end) <= flank)
    if not ok:
        return None
    return {
        "start": tir_start, "end": tir_end,
        "left": (tir_start, lw0 + i1), "right": (rw1 - j1, tir_end),
        "left_len": left_len, "right_len": right_len,
        "identity": round(identity, 3), "score": score,
        "matches": matches, "columns": columns,
        "left_seq": seq[tir_start:lw0 + i1].upper(), "right_seq": seq[rw1 - j1:tir_end].upper(),
    }


#------------------------------------------------------------------- DR -----

def find_dr(left, right, min_len=2, max_len=15, mm_from=8):
    """Longest direct repeat between the sequence just left of an element and
    the sequence just right of it.

    ``left`` ends at the element start, ``right`` starts at the element end
    (pass at least ``max_len`` bp of each). For each k from max_len down to
    min_len, ``left[-k:]`` is compared with ``right[:k]``; k >= ``mm_from``
    may have one mismatch, shorter repeats must be exact.

    Returns (k, mismatches, left_copy, right_copy) or (0, 0, "", "").
    """
    left, right = left.upper(), right.upper()
    for k in range(min(max_len, len(left), len(right)), min_len - 1, -1):
        a, b = left[len(left) - k:], right[:k]
        mm = sum(1 for x, y in zip(a, b) if x != y or x == "N")
        if mm == 0 or (k >= mm_from and mm <= 1):
            return k, mm, a, b
    return 0, 0, "", ""


def genome_dr(contig_seq, start, end, max_len=15, **kw):
    return find_dr(contig_seq[max(0, start - max_len):start], contig_seq[end:end + max_len],
                   max_len=max_len, **kw)


#------------------------------------------------------------------ ORF -----

def longest_orf(seq):
    """Longest complete ORF (ATG ... stop) on either strand of seq.

    Returns dict(aa, start, end, strand) with start/end 0-based half-open on
    seq (stop codon included; for strand '-' the ORF is revcomp(seq[start:end]))
    or None.
    """
    seq = seq.upper()
    n = len(seq)
    best = None
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            open_at = None
            for i in range(frame, n - 2, 3):
                codon = s[i:i + 3]
                if open_at is None:
                    if codon == "ATG":
                        open_at = i
                elif codon in STOPS:
                    length = i + 3 - open_at
                    if best is None or length > best[0]:
                        best = (length, open_at, i + 3, strand)
                    open_at = None
    if best is None:
        return None
    length, a, b, strand = best
    if strand == "-":
        a, b = n - b, n - a
    return {"aa": length // 3 - 1, "start": a, "end": b, "strand": strand}


#---------------------------------------------------------------- input -----

def parse_query_id(query_id):
    """'{read_id}_{start}_{end}' -> (read_id, start, end); read IDs may contain '_'."""
    read_id, start, end = query_id.rsplit("_", 2)
    return read_id, int(start), int(end)


def discover_samples(data_dir):
    """[(sample, best_alignment.tab, reads fasta)] for every usable sample folder."""
    found = []
    for name in sorted(os.listdir(data_dir)):
        folder = os.path.join(data_dir, name)
        if not os.path.isdir(folder):
            continue
        annot = os.path.join(folder, ANNOT_NAME)
        fasta = os.path.join(folder, f"{name}.fasta")
        ins = os.path.join(folder, INSERT_TAB_NAME)
        if os.path.isfile(annot) and os.path.isfile(fasta):
            if not os.path.isfile(ins):
                logger.warning("%s: no %s (only best_alignment.tab is needed; continuing)", name, INSERT_TAB_NAME)
            found.append((name, annot, fasta))
        elif os.path.isfile(annot) or os.path.isfile(ins):
            logger.warning("%s: missing %s or %s.fasta, skipping", name, ANNOT_NAME, name)
    return found


def read_genome_inserts(samples):
    """Collect the sequences of inserts labelled Genome (see read_inserts)."""
    return read_inserts(samples, ("Genome",))


def read_inserts(samples, notes=("Genome",)):
    """Collect the sequences of inserts whose Note is in ``notes``.

    ``samples``: [(sample, annotation tab, reads fasta)].
    Returns a list of dicts (sample, query_id, read_id, qstart, qend, ref_pos,
    note, seq).
    """
    notes = set(notes)
    label = "/".join(sorted(notes))
    inserts = []
    for sample, annot, fasta in samples:
        wanted = []
        with open(annot, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("Note") not in notes:
                    continue
                try:
                    read_id, qs, qe = parse_query_id(row["Query ID"])
                except ValueError:
                    logger.warning("%s: cannot parse Query ID %r", sample, row["Query ID"])
                    continue
                wanted.append((row["Query ID"], read_id, qs, qe, row.get("ref_pos", ""), row["Note"]))
        if not wanted:
            logger.info("%s: no %s inserts", sample, label)
            continue
        reads = SeqIO.index(fasta, "fasta")
        try:
            n = 0
            for qid, read_id, qs, qe, ref_pos, note in wanted:
                if read_id not in reads:
                    logger.warning("%s: read %s not in %s", sample, read_id, fasta)
                    continue
                seq = str(reads[read_id].seq[qs:qe]).upper()
                if seq:
                    inserts.append({"sample": sample, "query_id": qid, "read_id": read_id,
                                    "qstart": qs, "qend": qe, "ref_pos": ref_pos, "note": note,
                                    "seq": seq})
                    n += 1
        finally:
            reads.close()
        logger.info("%s: %d %s inserts", sample, n, label)
    return inserts


def read_plasmid_drs(dr_dir, sample_names):
    """Plasmid target-site duplications from the DR finder
    (``<dr_dir>/<sample>_insertions.csv``):
    {(sample, read_id, qstart, qend): (repeat length, DR sequence)}.
    Rows without a repeat length are skipped; missing files are ignored."""
    drs = {}
    if not dr_dir or not os.path.isdir(dr_dir):
        return drs
    for sample in sample_names:
        path = os.path.join(dr_dir, f"{sample}_insertions.csv")
        if not os.path.isfile(path):
            continue
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    key = (sample, row["Read ID"], int(float(row["Query Start CSV"])),
                           int(float(row["Query End CSV"])))
                    length = int(float(row["Repeat Length"]))
                except (KeyError, TypeError, ValueError):
                    continue
                seq = (row.get("Overlap Sequence") or "").strip()
                drs[key] = (length, "" if seq in ("", "N/A", "nan") else seq.upper())
    return drs


def plasmid_dr_summary(payloads, drs):
    """Modal plasmid DR of a cluster's inserts: dict(plasmid_dr_len,
    plasmid_dr_seq, plasmid_dr_support, plasmid_dr_counts). plasmid_dr_support
    = fraction of inserts with a DR-finder repeat length that have the modal
    length. All values are NA without DR-finder output."""
    vals = [drs.get((p["sample"], p["read_id"], p["qstart"], p["qend"])) for p in payloads]
    vals = [v for v in vals if v is not None]
    if not vals:
        return {"plasmid_dr_len": NA, "plasmid_dr_seq": NA, "plasmid_dr_support": NA,
                "plasmid_dr_counts": NA}
    counts = Counter(v[0] for v in vals)
    mode, n = max(counts.items(), key=lambda x: (x[1], x[0]))
    seqs = Counter(v[1] for v in vals if v[0] == mode and v[1])
    return {"plasmid_dr_len": mode, "plasmid_dr_seq": seqs.most_common(1)[0][0] if seqs else "",
            "plasmid_dr_support": round(n / len(vals), 3),
            "plasmid_dr_counts": ";".join(f"{k}:{c}" for k, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])))}


#-------------------------------------------------------------- mapping -----

def _run(cmd, **kw):
    logger.debug("Running: %s", " ".join(cmd))
    try:
        return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, **kw)
    except FileNotFoundError:
        raise RuntimeError(f"'{cmd[0]}' not found in PATH") from None
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{cmd[0]} failed (exit {e.returncode}): {e.stderr.strip()[-2000:]}") from None


def parse_paf(path):
    """Yield PAF records as dicts (only the fields used here)."""
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = dict(t.split(":", 1) for t in f[12:] if t.count(":") >= 2)
            yield {"qname": f[0], "qlen": int(f[1]), "qstart": int(f[2]), "qend": int(f[3]),
                   "strand": f[4], "tname": f[5], "tlen": int(f[6]), "tstart": int(f[7]),
                   "tend": int(f[8]), "matches": int(f[9]), "alen": int(f[10]), "mapq": int(f[11]),
                   "primary": tags.get("tp", "A:P") == "A:P"}


def map_inserts(inserts, genome_fasta, work_dir, threads=2):
    """Map inserts with minimap2; set insert['hit'] to the best primary alignment."""
    query = os.path.join(work_dir, "inserts.fasta")
    paf = os.path.join(work_dir, "inserts.paf")
    with open(query, "w") as fh:
        for i, ins in enumerate(inserts):
            fh.write(f">i{i}\n{ins['seq']}\n")
    _run(["minimap2", "-c", "-x", "map-ont", "--secondary=no", "-t", str(max(1, int(threads))),
          "-o", paf, os.path.abspath(genome_fasta), query])
    best = {}
    for rec in parse_paf(paf):
        if not rec["primary"]:
            continue
        cur = best.get(rec["qname"])
        if cur is None or rec["matches"] > cur["matches"]:
            best[rec["qname"]] = rec
    for i, ins in enumerate(inserts):
        ins["hit"] = best.get(f"i{i}")
    return inserts


#----------------------------------------------------------- clustering -----

def reciprocal_overlap(a0, a1, b0, b1):
    ov = min(a1, b1) - max(a0, b0)
    if ov <= 0:
        return 0.0
    return min(ov / (a1 - a0), ov / (b1 - b0))


def cluster_intervals(items, min_overlap=0.5, contain=0.9):
    """Greedy clustering of genome intervals.

    ``items``: list of (contig, start, end, payload). Longest intervals seed
    clusters; an interval joins the cluster (same contig) whose seed interval
    it overlaps reciprocally >= min_overlap (best overlap wins). Afterwards a
    cluster whose span lies >= ``contain`` inside a larger cluster's seed is
    absorbed (its members are flagged partial).
    Returns a list of clusters: dict(contig, seed=(s, e), members=[...], partial=[...]).
    """
    clusters = []
    for contig, s, e, payload in sorted(items, key=lambda x: -(x[2] - x[1])):
        best, best_ov = None, 0.0
        for c in clusters:
            if c["contig"] != contig:
                continue
            ov = reciprocal_overlap(s, e, *c["seed"])
            if ov >= min_overlap and ov > best_ov:
                best, best_ov = c, ov
        if best is None:
            clusters.append({"contig": contig, "seed": (s, e), "members": [(s, e, payload)], "partial": []})
        else:
            best["members"].append((s, e, payload))
    # absorb clusters of partial inserts into the enclosing cluster
    clusters.sort(key=lambda c: -(c["seed"][1] - c["seed"][0]))
    kept = []
    for c in clusters:
        lo = min(m[0] for m in c["members"])
        hi = max(m[1] for m in c["members"])
        host = None
        for k in kept:
            if k["contig"] != c["contig"]:
                continue
            ov = min(hi, k["seed"][1]) - max(lo, k["seed"][0])
            if ov >= contain * (hi - lo):
                host = k
                break
        if host is None:
            kept.append(c)
        else:
            host["partial"].extend(c["members"] + c["partial"])
    return kept


def robust_end(values, window=10):
    """Mode (or, if the mode is tied, the median) of the values inside the
    densest window of `window` bp, and the median absolute deviation of all
    values from it."""
    vals = sorted(values)
    best_lo, best_n, j = 0, 0, 0
    for i in range(len(vals)):
        while vals[i] - vals[j] > window:
            j += 1
        if i - j + 1 > best_n:
            best_n, best_lo = i - j + 1, j
    core = vals[best_lo:best_lo + best_n]
    counts = Counter(core).most_common()
    top = [v for v, n in counts if n == counts[0][1]]
    est = int(round(statistics.median(top if len(top) == 1 else core)))
    spread = statistics.median(abs(v - est) for v in vals)
    return est, spread


def boundary_candidates(values, est, radius=3, min_support=2, min_fraction=0.2):
    """Alternative values for one element end: values within `radius` bp of
    the estimate observed in >= max(min_support, min_fraction * mode count)
    reads. Returns {value: read count}, always including `est`."""
    counts = Counter(v for v in values if abs(v - est) <= radius)
    top = max(counts.values()) if counts else 0
    cands = {v: n for v, n in counts.items() if n >= max(min_support, min_fraction * top)}
    cands.setdefault(est, counts.get(est, 0))
    return cands


def _site_summary(positions, tolerance=20):
    """'pos:count;...' sorted by count, and the number of distinct sites
    (positions within `tolerance` bp are one site)."""
    counts = Counter(positions)
    text = ";".join(f"{p}:{n}" for p, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])))
    ints = sorted(int(p) for p in counts if str(p).lstrip("-").isdigit())
    n_sites, last = 0, None
    for p in ints:
        if last is None or p - last > tolerance:
            n_sites += 1
        last = p
    return text, n_sites


#------------------------------------------------------- characterise -----

def _exact_dr(contig_seq, start, end):
    return genome_dr(contig_seq, start, end, mm_from=99)[0]


def characterise(contig_seq, est_start, est_end, flank_search=60, max_tir=50, min_tir_len=10,
                 min_tir_identity=0.75, min_tir_score=22, edge_run=2, trim_mismatch=-2,
                 dr_shift=3, tir_override_gain=3, tir_anchor=3, min_dr_len=3, min_orf_aa=200,
                 candidates=None):
    """Boundaries, TIR, DR, ORF and verdict for one element candidate.

    ``est_start``/``est_end``: boundaries estimated from the read alignments.
    ``candidates``: optional ({start: reads}, {end: reads}) of alternative
    read-supported ends (see boundary_candidates); the combination with the
    longest exact genomic DR is used, ties going to the most reads and then to
    the estimate.
    Returns (row dict, element sequence). The boundary rules are described
    in step 6 of the module docstring.
    """
    tir_kw = dict(max_tir=max_tir, min_len=min_tir_len, min_identity=min_tir_identity,
                  min_score=min_tir_score, edge_run=edge_run, trim_mismatch=trim_mismatch)
    start, end, source = est_start, est_end, "reads"
    wide = find_tir(contig_seq, est_start, est_end, flank=flank_search, **tir_kw)
    if wide:
        # TIR edges can be off by a few bp (chance pairing of flank bases, or a
        # TIR mismatch near the edge): try symmetric shifts t in
        # [-dr_shift, dr_shift] (t > 0 = inward) and keep the longest EXACT
        # genomic DR, ties to the smallest |t|. Exact DRs are used because with
        # one mismatch allowed, shift t+1 at length k+2 contains the true DR of
        # shift t plus one random pair on each side.
        best = None
        for t in sorted(range(-dr_shift, dr_shift + 1), key=abs):
            ts, te = wide["start"] + t, wide["end"] - t
            if te - ts <= 2 * min_tir_len:
                continue
            k = _exact_dr(contig_seq, ts, te)
            if best is None or k > best[2]:
                best = (ts, te, k)
        if best is not None:
            ts, te, k_tir = best
            if (ts, te) == (est_start, est_end):
                source = "reads+TIR"
            elif k_tir >= _exact_dr(contig_seq, est_start, est_end) + tir_override_gain:
                start, end, source = ts, te, "TIR"
    if source == "reads" and candidates:
        # Several ends are supported by reads (e.g. all reads of one insertion
        # site overshoot by 1 bp because a plasmid flank base matches the
        # genome). The genomic DR - the trace of the transposition - decides.
        cs, ce = candidates
        best = None
        for s0, ns in cs.items():
            for e0, ne in ce.items():
                if e0 - s0 <= 2 * min_tir_len:
                    continue
                key = (_exact_dr(contig_seq, s0, e0), ns + ne, (s0, e0) == (est_start, est_end))
                if best is None or key > best[0]:
                    best = (key, s0, e0)
        if best is not None and (best[1], best[2]) != (start, end):
            start, end, source = best[1], best[2], "reads+DR"
    tir = find_tir(contig_seq, start, end, flank=tir_anchor, **tir_kw) or wide
    dr = genome_dr(contig_seq, start, end)
    element = contig_seq[start:end]
    orf = longest_orf(element)
    orf_aa = orf["aa"] if orf else 0
    has_dr = dr[0] >= min_dr_len
    has_orf = orf_aa >= min_orf_aa
    if tir and (has_dr or has_orf):
        verdict = "mobile_element_candidate"
    elif tir or (has_dr and has_orf):
        verdict = "partial_evidence"
    else:
        verdict = "genome_segment"
    if tir:
        tir_cols = {"tir_left_len": tir["left_len"], "tir_right_len": tir["right_len"],
                    "tir_left_offset": tir["start"] - start, "tir_right_offset": end - tir["end"],
                    "tir_identity": tir["identity"], "tir_score": tir["score"],
                    "tir_left_seq": tir["left_seq"], "tir_right_seq": tir["right_seq"]}
    else:
        tir_cols = {"tir_left_len": 0, "tir_right_len": 0, "tir_left_offset": "",
                    "tir_right_offset": "", "tir_identity": "", "tir_score": "",
                    "tir_left_seq": "", "tir_right_seq": ""}
    row = {
        "start": start, "end": end, "length": end - start, "boundary_source": source,
        "tir_found": bool(tir), **tir_cols,
        "dr_len": dr[0], "dr_mismatches": dr[1], "dr_left_seq": dr[2], "dr_right_seq": dr[3],
        "longest_orf_aa": orf_aa,
        "orf_start": start + orf["start"] if orf else "",
        "orf_end": start + orf["end"] if orf else "",
        "orf_strand": orf["strand"] if orf else "",
        "verdict": verdict,
    }
    return row, element


def count_genome_copies(fasta, genome_fasta, work_dir, threads=2, min_cov=0.9):
    """Number of genome alignments covering >= min_cov of each element."""
    paf = os.path.join(work_dir, "copies.paf")
    _run(["minimap2", "-c", "-x", "asm10", "-N", "1000", "-p", "0.5", "-t", str(max(1, int(threads))),
          "-o", paf, os.path.abspath(genome_fasta), os.path.abspath(fasta)])
    counts = Counter()
    for rec in parse_paf(paf):
        if (rec["qend"] - rec["qstart"]) >= min_cov * rec["qlen"]:
            counts[rec["qname"]] += 1
    return counts


def similar_clusters(fasta, work_dir, threads=2, min_cov=0.8):
    """{cluster: [other clusters]} for elements whose sequences align over
    >= min_cov of both (e.g. copies of one multi-copy element whose reads
    were split between genomic copies)."""
    paf = os.path.join(work_dir, "self.paf")
    _run(["minimap2", "-c", "-x", "asm10", "-N", "100", "-p", "0.5", "-t", str(max(1, int(threads))),
          "-o", paf, os.path.abspath(fasta), os.path.abspath(fasta)])
    sim = {}
    for rec in parse_paf(paf):
        if rec["qname"] == rec["tname"]:
            continue
        if (rec["qend"] - rec["qstart"] >= min_cov * rec["qlen"]
                and rec["tend"] - rec["tstart"] >= min_cov * rec["tlen"]):
            sim.setdefault(rec["qname"], set()).add(rec["tname"])
            sim.setdefault(rec["tname"], set()).add(rec["qname"])
    return {k: sorted(v) for k, v in sim.items()}


def known_element_hits(fasta, is_fasta, work_dir, threads=2):
    """Best hit of each element in a FASTA of known elements (e.g. the curated
    IS FASTA): {cluster: 'name:identity%:element_coverage%'}."""
    paf = os.path.join(work_dir, "known.paf")
    _run(["minimap2", "-c", "-x", "asm20", "-t", str(max(1, int(threads))),
          "-o", paf, os.path.abspath(is_fasta), os.path.abspath(fasta)])
    best = {}
    for rec in parse_paf(paf):
        if rec["qname"] not in best or rec["matches"] > best[rec["qname"]]["matches"]:
            best[rec["qname"]] = rec
    return {q: f"{r['tname']}:{100 * r['matches'] / r['alen']:.1f}:"
               f"{100 * (r['qend'] - r['qstart']) / r['qlen']:.1f}" for q, r in best.items()}


#----------------------------------------------- Unannotated: read-based -----

def all_vs_all(seqs, work_dir, threads=2, name="ava"):
    """minimap2 -x ava-ont over ``seqs`` (list of str). Returns
    {(i, j): dict(strand, cov={i: .., j: ..}, matches, alen)} with i < j, the
    best record (most matching bases) per pair; cov = aligned fraction of each
    sequence."""
    fa = os.path.join(work_dir, f"{name}.fasta")
    paf = os.path.join(work_dir, f"{name}.paf")
    with open(fa, "w") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f">{i}\n{seq}\n")
    _run(["minimap2", "-x", "ava-ont", "-t", str(max(1, int(threads))), "-o", paf, fa, fa])
    pairs = {}
    for rec in parse_paf(paf):
        a, b = int(rec["qname"]), int(rec["tname"])
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key in pairs and pairs[key]["matches"] >= rec["matches"]:
            continue
        pairs[key] = {"strand": rec["strand"], "matches": rec["matches"], "alen": rec["alen"],
                      "cov": {a: (rec["qend"] - rec["qstart"]) / rec["qlen"],
                              b: (rec["tend"] - rec["tstart"]) / rec["tlen"]}}
    return pairs


def _pair(pairs, i, j):
    return pairs.get((min(i, j), max(i, j)))


def cluster_by_similarity(lengths, pairs, min_overlap=0.5, contain=0.9):
    """Cluster sequences without a reference (cf. cluster_intervals).

    Sequences are linked when their alignment covers >= ``min_overlap`` of
    BOTH (reciprocal coverage); clusters are the connected components. Single
    linkage tolerates the pairs that ava-ont misses or under-reports at
    nanopore error rates (a seed-only greedy scheme split one element into up
    to 11 clusters on simulated r9.4.1 reads). It cannot chain two elements
    through a chimeric insert A+B: its coverage by A is a/(a+b) and by B
    b/(a+b), both >= 0.5 only when a == b (it joins the element of its
    longer part). A component with fewer members than another component is
    absorbed into it as partial when one of its members is covered >=
    ``contain`` by a typical member of the larger one (length <= 1.1 x the
    median member length, so a chimera containing an element cannot absorb
    that element's cluster) - fragmented inserts. Returns [dict(seed=longest
    member, members=[...], partial=[...])], largest first.
    """
    n = len(lengths)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, j), p in pairs.items():
        if min(p["cov"][i], p["cov"][j]) >= min_overlap:
            parent[find(i)] = find(j)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    clusters = [{"seed": max(m, key=lambda k: (lengths[k], -k)), "members": sorted(m), "partial": []}
                for m in comps.values()]
    clusters.sort(key=lambda c: (-len(c["members"]), -lengths[c["seed"]], c["seed"]))
    where = {}
    typical_max = []
    for ci, c in enumerate(clusters):
        for i in c["members"]:
            where[i] = ci
        typical_max.append(1.1 * statistics.median(lengths[i] for i in c["members"]))
    host_of = {}
    for (i, j), p in pairs.items():
        for small, big in ((i, j), (j, i)):
            cs, cb = where[small], where[big]
            if cs == cb or len(clusters[cb]["members"]) <= len(clusters[cs]["members"]):
                continue
            if lengths[big] > typical_max[cb]:
                continue
            if p["cov"][small] >= contain:
                cur = host_of.get(cs)
                if cur is None or cb < cur:      # prefer the largest host
                    host_of[cs] = cb
    kept = []
    for ci, c in enumerate(clusters):
        h = ci
        while h in host_of:                       # follow absorption chains to the final host
            h = host_of[h]
        if h == ci:
            kept.append(c)
        else:
            clusters[h]["partial"].extend(c["members"] + c["partial"])
    return kept


def medoid(members, pairs):
    """Member with the most matching bases to the other members (ties: lowest index)."""
    members = set(members)
    score = Counter()
    for (i, j), p in pairs.items():
        if i in members and j in members:
            score[i] += p["matches"]
            score[j] += p["matches"]
    return max(sorted(members), key=lambda i: score[i])


def majority_consensus(aligned, min_cov_fraction=0.5, base_fraction=0.5):
    """Column-wise vote over an MSA (list of equal-length strings).

    Leading/trailing gaps of a sequence are 'no data', not votes for a gap, so
    fragmented members do not shorten the consensus. A column needs votes from
    >= ``min_cov_fraction`` of the sequences (at least 2 when there are >= 2).
    The column is kept when >= ``base_fraction`` of its votes are bases (0.5:
    plain majority, a base beating a gap on ties); its base is the most
    frequent one. Returns the consensus string.
    """
    if not aligned:
        return ""
    if len(aligned) == 1:
        return aligned[0].replace("-", "").upper()
    spans = []
    for a in aligned:
        lo = len(a) - len(a.lstrip("-"))
        spans.append((lo, len(a.rstrip("-"))))
    need = max(2, math.ceil(min_cov_fraction * len(aligned)))
    out = []
    for col in range(len(aligned[0])):
        votes = Counter(a[col].upper() for a, (lo, hi) in zip(aligned, spans) if lo <= col < hi)
        if sum(votes.values()) < need:
            continue
        total = sum(votes.values())
        votes.pop("-", 0)
        votes.pop("N", None)
        if not votes:
            continue
        base, _n = max(votes.items(), key=lambda x: (x[1], x[0]))
        if sum(votes.values()) >= base_fraction * total:
            out.append(base)
    return "".join(out)


def _parse_cs(cs):
    """minimap2 short cs string -> list of (op, value): (':', n), ('*', 'ac'),
    ('+', 'seq'), ('-', 'seq')."""
    ops, i, n = [], 0, len(cs)
    while i < n:
        op = cs[i]
        j = i + 1
        if op == ":":
            while j < n and cs[j].isdigit():
                j += 1
            ops.append((":", int(cs[i + 1:j])))
        elif op == "*":
            j = i + 3
            ops.append(("*", cs[i + 1:j]))
        else:
            while j < n and cs[j].isalpha():
                j += 1
            ops.append((op, cs[i + 1:j]))
        i = j
    return ops


def pileup_polish(draft, seqs, work_dir, threads=2, base_fraction=0.5, name="polish"):
    """One round of polishing: ``seqs`` are mapped to ``draft`` (minimap2 -c
    --cs, map-ont, best alignment per read) and every draft position is
    re-voted: the base is kept when >= ``base_fraction`` of the reads
    covering it show a base there (most frequent base wins), an insertion
    before it is added when > 1 - ``base_fraction`` of the reads spanning it
    carry one (most frequent inserted string). Indels are left-aligned in
    repeats so equivalent placements vote together. Positions no read covers
    stay as in the draft. Indels are counted against one sequence, so the
    misplaced indels an MSA can leave in a column-wise consensus are
    corrected. Returns the polished sequence."""
    if not seqs or not draft:
        return draft
    dfa = os.path.join(work_dir, f"{name}_draft.fasta")
    qfa = os.path.join(work_dir, f"{name}_reads.fasta")
    paf = os.path.join(work_dir, f"{name}.paf")
    with open(dfa, "w") as fh:
        fh.write(f">draft\n{draft}\n")
    with open(qfa, "w") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f">r{i}\n{seq}\n")
    _run(["minimap2", "-c", "--cs", "-x", "map-ont", "--secondary=no", "-t", str(max(1, int(threads))),
          "-o", paf, dfa, qfa])
    best = {}
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or "tp:A:P" not in f[12:]:
                continue
            cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), None)
            if cs is None:
                continue
            if f[0] not in best or int(f[9]) > best[f[0]][0]:
                best[f[0]] = (int(f[9]), int(f[7]), int(f[8]), cs)
    n = len(draft)
    base_votes = [Counter() for _ in range(n)]
    ins_votes = [Counter() for _ in range(n + 1)]
    cover = [0] * (n + 1)          # reads spanning the gap before position p (p-1 and p aligned)
    for _m, ts, te, cs in best.values():
        t = ts
        for op, val in _parse_cs(cs):
            if op == ":":
                for k in range(val):
                    base_votes[t + k][draft[t + k]] += 1
                t += val
            elif op == "*":
                base_votes[t][val[1].upper()] += 1
                t += 1
            elif op == "-":
                # left-align: equivalent placements of a deletion in a repeat get one vote position
                d0, ln = t, len(val)
                while d0 > ts and draft[d0 - 1] == draft[d0 + ln - 1]:
                    d0 -= 1
                for k in range(ln):
                    base_votes[d0 + k]["-"] += 1
                t += ln
            elif op == "+":
                i0, seq = t, val.upper()
                while i0 > ts + 1 and draft[i0 - 1] == seq[-1]:      # left-align the insertion
                    seq = seq[-1] + seq[:-1]
                    i0 -= 1
                ins_votes[i0][seq] += 1
        for p in range(ts + 1, te):
            cover[p] += 1
    out = []
    for p in range(n):
        iv = ins_votes[p]
        if cover[p] and sum(iv.values()) > (1 - base_fraction) * cover[p]:
            out.append(iv.most_common(1)[0][0])
        bv = base_votes[p]
        total = sum(bv.values())
        if not total:
            out.append(draft[p])
            continue
        gaps = bv.pop("-", 0)
        if bv and total - gaps >= base_fraction * total:
            out.append(max(bv.items(), key=lambda x: (x[1], x[0]))[0])
    return "".join(out)


def build_consensus(seqs, work_dir, threads=2, name="cons", base_fraction=0.4, polish_rounds=2):
    """Consensus of oriented sequences (first = medoid): MAFFT MSA ->
    majority_consensus (draft) -> ``polish_rounds`` rounds of pileup_polish.
    Without MAFFT in PATH the medoid is the draft (warning)."""
    if len(seqs) == 1:
        return seqs[0]
    if shutil.which("mafft") is None:
        logger.warning("mafft not in PATH: the medoid insert is the draft consensus")
        draft = seqs[0]
    else:
        draft = _msa_consensus(seqs, work_dir, threads, name, base_fraction)
    for r in range(polish_rounds):
        polished = pileup_polish(draft, seqs, work_dir, threads, base_fraction, f"{name}_p{r}")
        if polished == draft:
            break
        draft = polished
    return draft


def _msa_consensus(seqs, work_dir, threads, name, base_fraction):
    fa = os.path.join(work_dir, f"{name}.fasta")
    with open(fa, "w") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f">m{i}\n{seq}\n")
    res = _run(["mafft", "--auto", "--nuc", "--quiet", "--thread", str(max(1, int(threads))), fa])
    aligned = {}
    name_ = None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            name_ = line[1:].split()[0]
            aligned[name_] = []
        elif name_ is not None:
            aligned[name_].append(line.strip())
    msa = ["".join(aligned[f"m{i}"]) for i in range(len(seqs)) if f"m{i}" in aligned]
    return majority_consensus(msa, base_fraction=base_fraction) or seqs[0]


def map_identity(queries, ref_fasta, work_dir, threads=2):
    """Identity (matches / alignment columns, minimap2 -c map-ont) of each
    query (name, seq, target name) to its own target in ``ref_fasta``:
    {name: identity}; queries without an alignment to their target are absent."""
    fa = os.path.join(work_dir, "members.fasta")
    paf = os.path.join(work_dir, "members.paf")
    with open(fa, "w") as fh:
        for qname, seq, _t in queries:
            fh.write(f">{qname}\n{seq}\n")
    _run(["minimap2", "-c", "-x", "map-ont", "-N", "50", "-t", str(max(1, int(threads))),
          "-o", paf, os.path.abspath(ref_fasta), fa])
    target = {q: t for q, _s, t in queries}
    best = {}
    for rec in parse_paf(paf):
        if target.get(rec["qname"]) != rec["tname"]:
            continue
        if rec["qname"] not in best or rec["matches"] > best[rec["qname"]]["matches"]:
            best[rec["qname"]] = rec
    return {q: r["matches"] / r["alen"] for q, r in best.items() if r["alen"]}


def characterise_consensus(cons, flank=20, max_tir=50, min_tir_len=10, min_tir_identity=0.75,
                           min_tir_score=22, edge_run=2, trim_mismatch=-2, min_dr_len=3,
                           min_orf_aa=200, plasmid_dr=None, min_dr_support=0.5, **_ignored):
    """TIR, boundaries, ORF and verdict for a read consensus (no genome copy).

    The consensus spans the insert as cut out of the reads, i.e. the element
    +- a few bp. TIRs are searched within ``flank`` bp of both ends; when one
    is found its outer edges become the element ends (boundary_source
    ``consensus+TIR``), otherwise the whole consensus is the element
    (``consensus``). There is no genomic DR (NA). The verdict uses the plasmid
    target-site duplication instead (``plasmid_dr``: plasmid_dr_summary of
    the cluster): DR evidence = modal plasmid DR >= min_dr_len in >=
    ``min_dr_support`` of the inserts with a DR-finder result.
    Returns (row dict, element sequence).
    """
    tir_kw = dict(max_tir=max_tir, min_len=min_tir_len, min_identity=min_tir_identity,
                  min_score=min_tir_score, edge_run=edge_run, trim_mismatch=trim_mismatch)
    n = len(cons)
    tir = find_tir(cons, 0, n, flank=flank, **tir_kw)
    if tir:
        start, end, source = tir["start"], tir["end"], "consensus+TIR"
    else:
        start, end, source = 0, n, "consensus"
    element = cons[start:end]
    orf = longest_orf(element)
    orf_aa = orf["aa"] if orf else 0
    pdr = plasmid_dr or {}
    has_dr = (isinstance(pdr.get("plasmid_dr_len"), int) and pdr["plasmid_dr_len"] >= min_dr_len
              and (pdr.get("plasmid_dr_support") or 0) >= min_dr_support)
    has_orf = orf_aa >= min_orf_aa
    if tir and (has_dr or has_orf):
        verdict = "mobile_element_candidate"
    elif tir or (has_dr and has_orf):
        verdict = "partial_evidence"
    else:
        verdict = "no_element_evidence"
    if tir:
        tir_cols = {"tir_left_len": tir["left_len"], "tir_right_len": tir["right_len"],
                    "tir_left_offset": 0, "tir_right_offset": 0,
                    "tir_identity": tir["identity"], "tir_score": tir["score"],
                    "tir_left_seq": tir["left_seq"], "tir_right_seq": tir["right_seq"]}
    else:
        tir_cols = {"tir_left_len": 0, "tir_right_len": 0, "tir_left_offset": "",
                    "tir_right_offset": "", "tir_identity": "", "tir_score": "",
                    "tir_left_seq": "", "tir_right_seq": ""}
    row = {
        "start": start, "end": end, "length": end - start, "boundary_source": source,
        "tir_found": bool(tir), **tir_cols,
        "dr_len": NA, "dr_mismatches": NA, "dr_left_seq": NA, "dr_right_seq": NA,
        "longest_orf_aa": orf_aa,
        "orf_start": orf["start"] if orf else "",
        "orf_end": orf["end"] if orf else "",
        "orf_strand": orf["strand"] if orf else "",
        "verdict": verdict,
    }
    return row, element


def read_clusters(inserts, work_dir, threads=2, min_reads=2, min_overlap=0.5, site_tolerance=20,
                  max_consensus_reads=30, drs=None, consensus_flank=20, cparams=None,
                  base_fraction=0.4, polish_rounds=2):
    """Cluster inserts without genome hit by all-vs-all similarity and
    characterise a consensus per cluster. Returns (rows, fasta records)."""
    cparams = cparams or {}
    pairs = all_vs_all([x["seq"] for x in inserts], work_dir, threads, "ava")
    lengths = [len(x["seq"]) for x in inserts]
    clusters = [c for c in cluster_by_similarity(lengths, pairs, min_overlap)
                if len(c["members"]) >= min_reads]
    rows, records, idq = [], [], []
    for k, c in enumerate(clusters, 1):
        cid = f"unannotated_{k:03d}"
        med = medoid(c["members"], pairs)
        # orient members on the medoid; members without a direct alignment are not used
        oriented = [(med, "+")]
        others = []
        for i in c["members"]:
            if i == med:
                continue
            p = _pair(pairs, i, med)
            if p is not None:
                others.append((-p["matches"], i, p["strand"]))
        oriented += [(i, s) for _m, i, s in sorted(others)][:max(0, max_consensus_reads - 1)]
        seqs = [inserts[i]["seq"] if s == "+" else revcomp(inserts[i]["seq"]) for i, s in oriented]
        cons = build_consensus(seqs, work_dir, threads, f"cons_{k}", base_fraction, polish_rounds)
        payloads = [inserts[i] for i in c["members"] + c["partial"]]
        strands = Counter()
        for i in c["members"] + c["partial"]:
            p = _pair(pairs, i, med)
            strands["+" if i == med else (p["strand"] if p else "?")] += 1
        pdr = plasmid_dr_summary(payloads, drs or {})
        char, element = characterise_consensus(cons, flank=consensus_flank, plasmid_dr=pdr, **cparams)
        sample_counts = Counter(p["sample"] for p in payloads)
        sites, n_sites = _site_summary([p["ref_pos"] for p in payloads], site_tolerance)
        row = {"cluster_id": cid, "source": "reads", "contig": NA,
               "n_reads": len(c["members"]), "n_partial_reads": len(c["partial"]),
               "n_samples": len(sample_counts),
               "samples": ";".join(f"{s}:{n}" for s, n in sorted(sample_counts.items())),
               "insert_strands": ";".join(f"{s}:{n}" for s, n in sorted(strands.items())),
               "est_start": NA, "est_end": NA, "start_spread": NA, "end_spread": NA,
               "consensus_reads": len(oriented), "consensus_length": len(cons),
               "n_plasmid_sites": n_sites, "plasmid_sites": sites, "genome_copies": "",
               "similar_clusters": "", "known_hit": "", **pdr, **char}
        rows.append(row)
        records.append((cid, f"consensus len={char['length']} reads={len(c['members'])} "
                             f"verdict={char['verdict']}", element))
        idq += [(f"{cid}|{i}", inserts[i]["seq"], cid) for i in c["members"]]
    if records:
        fa = os.path.join(work_dir, "consensus.fasta")
        with open(fa, "w") as fh:
            for cid, _d, seq in records:
                fh.write(f">{cid}\n{seq}\n")
        ident = map_identity(idq, fa, work_dir, threads)
        for row in rows:
            vals = [v for q, v in ident.items() if q.split("|")[0] == row["cluster_id"]]
            row["read_identity"] = round(statistics.median(vals), 4) if vals else NA
    return rows, records


#----------------------------------------------------------------- main -----

def main_novel_elements(data_dir=None, genome_fasta=None, out_dir=None, min_reads=2, threads=2,
                        flank_search=60, samples=None, min_overlap=0.5, end_window=10,
                        site_tolerance=20, is_fasta=None, keep_temp=False, unannotated=True,
                        max_consensus_reads=30, consensus_flank=20, dr_dir=None,
                        consensus_base_fraction=0.4, polish_rounds=2, **characterise_params):
    """Cluster and characterise the Genome- and Unannotated-labelled inserts.

    Inputs: ``data_dir`` with one folder per sample (best_alignment.tab and
    <sample>.fasta), or an explicit ``samples`` list of
    (sample, best_alignment.tab, reads fasta). ``genome_fasta`` must be the
    genome used for annotation.

    Writes novel_elements.tsv / .fasta / .log in out_dir (default
    <data_dir>/novel_elements) and returns the list of row dicts: first the
    genome clusters (``source`` genome, ``novel_NNN``), then the read clusters
    (``source`` reads, ``unannotated_NNN``; skipped with unannotated=False).

    ``is_fasta`` (optional, e.g. the curated IS FASTA): each element is aligned
    to it and the best hit is reported in ``known_hit`` - useful because an
    insert of a known IS can still be labelled Genome (e.g. when plasmid
    overhang lowers its IS coverage below the annotation threshold).

    ``min_reads`` counts reads whose insert overlaps the cluster reciprocally
    (absorbed partial inserts are reported in n_partial_reads but not counted).

    Read clusters: ``max_consensus_reads`` members closest to the medoid
    build the consensus; ``consensus_flank`` bp at each consensus end are
    searched for TIRs. Plasmid DRs come from the DR finder tables in
    ``dr_dir`` (default <data_dir>/insertions; NA when absent).

    Further keyword arguments override CHARACTERISE_DEFAULTS (max_tir,
    min_tir_len, min_tir_identity, min_tir_score, edge_run, trim_mismatch,
    dr_shift, tir_override_gain, tir_anchor, min_dr_len, min_orf_aa).
    """
    unknown = set(characterise_params) - set(CHARACTERISE_DEFAULTS)
    if unknown:
        raise TypeError(f"unknown parameter(s): {sorted(unknown)}")
    cparams = {**CHARACTERISE_DEFAULTS, **characterise_params, "flank_search": flank_search}
    if genome_fasta is None or not os.path.isfile(genome_fasta):
        raise FileNotFoundError(f"genome FASTA not found: {genome_fasta}")
    if samples is None:
        if data_dir is None or not os.path.isdir(data_dir):
            raise FileNotFoundError(f"data directory not found: {data_dir}")
        samples = discover_samples(data_dir)
    if out_dir is None:
        if data_dir is None:
            raise ValueError("out_dir is required when data_dir is not given")
        out_dir = os.path.join(data_dir, "novel_elements")
    if dr_dir is None and data_dir is not None:
        dr_dir = os.path.join(data_dir, "insertions")
    os.makedirs(out_dir, exist_ok=True)
    tsv_path = os.path.join(out_dir, "novel_elements.tsv")
    fasta_path = os.path.join(out_dir, "novel_elements.fasta")
    log_path = os.path.join(out_dir, "novel_elements.log")

    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    prev_level = logger.level
    if logger.getEffectiveLevel() > logging.INFO:
        logger.setLevel(logging.INFO)
    work = tempfile.mkdtemp(prefix="novel_", dir=out_dir)
    try:
        params = dict(min_reads=min_reads, min_overlap=min_overlap, end_window=end_window,
                      unannotated=unannotated, max_consensus_reads=max_consensus_reads,
                      consensus_flank=consensus_flank, consensus_base_fraction=consensus_base_fraction,
                      polish_rounds=polish_rounds, **cparams)
        logger.info("novel_elements: genome %s, %d sample(s), parameters %s", genome_fasta, len(samples), params)
        drs = read_plasmid_drs(dr_dir, [x[0] for x in samples])
        logger.info("plasmid DRs: %d DR-finder rows from %s", len(drs), dr_dir)
        inserts = read_inserts(samples, ("Genome",))
        logger.info("%d Genome inserts in total", len(inserts))
        rows, records = [], []
        if inserts:
            map_inserts(inserts, genome_fasta, work, threads)
            mapped = [x for x in inserts if x["hit"] is not None]
            logger.info("%d / %d inserts mapped to the genome", len(mapped), len(inserts))
            items = [(x["hit"]["tname"], x["hit"]["tstart"], x["hit"]["tend"], x) for x in mapped]
            clusters = cluster_intervals(items, min_overlap)
            genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(genome_fasta, "fasta")}
            clusters = [c for c in clusters if len(c["members"]) >= min_reads]
            clusters.sort(key=lambda c: (c["contig"], c["seed"][0]))
            for k, c in enumerate(clusters, 1):
                cid = f"novel_{k:03d}"
                starts = [m[0] for m in c["members"]]
                ends = [m[1] for m in c["members"]]
                est_s, sp_s = robust_end(starts, end_window)
                est_e, sp_e = robust_end(ends, end_window)
                cands = (boundary_candidates(starts, est_s), boundary_candidates(ends, est_e))
                payloads = [m[2] for m in c["members"]] + [m[2] for m in c["partial"]]
                char, element = characterise(genome[c["contig"]], est_s, est_e, candidates=cands,
                                             **cparams)
                sample_counts = Counter(p["sample"] for p in payloads)
                strands = Counter(p["hit"]["strand"] for p in payloads)
                sites, n_sites = _site_summary([p["ref_pos"] for p in payloads], site_tolerance)
                row = {"cluster_id": cid, "source": "genome", "contig": c["contig"],
                       "n_reads": len(c["members"]), "n_partial_reads": len(c["partial"]),
                       "n_samples": len(sample_counts),
                       "samples": ";".join(f"{s}:{n}" for s, n in sorted(sample_counts.items())),
                       "insert_strands": ";".join(f"{s}:{n}" for s, n in sorted(strands.items())),
                       "est_start": est_s, "est_end": est_e, "start_spread": sp_s, "end_spread": sp_e,
                       "consensus_reads": NA, "consensus_length": NA, "read_identity": NA,
                       "n_plasmid_sites": n_sites, "plasmid_sites": sites, "genome_copies": "",
                       "similar_clusters": "", "known_hit": "",
                       **plasmid_dr_summary(payloads, drs), **char}
                rows.append(row)
                records.append((cid, f"{c['contig']}:{char['start']}-{char['end']} "
                                     f"len={char['length']} verdict={char['verdict']}", element))
        if unannotated:
            un = read_inserts(samples, ("Unannotated",))
            logger.info("%d Unannotated inserts in total", len(un))
            if un:
                urows, urecs = read_clusters(un, work, threads, min_reads, min_overlap, site_tolerance,
                                             max_consensus_reads, drs, consensus_flank, cparams,
                                             consensus_base_fraction, polish_rounds)
                logger.info("%d Unannotated cluster(s) with >= %d reads", len(urows), min_reads)
                rows += urows
                records += urecs
        with open(fasta_path, "w") as fh:
            for cid, desc, seq in records:
                fh.write(f">{cid} {desc}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
        if records:
            copies = count_genome_copies(fasta_path, genome_fasta, work, threads)
            similar = similar_clusters(fasta_path, work, threads)
            known = known_element_hits(fasta_path, is_fasta, work, threads) if is_fasta else {}
            for row in rows:
                cid = row["cluster_id"]
                row["genome_copies"] = copies.get(cid, 0)
                row["similar_clusters"] = ";".join(similar.get(cid, []))
                row["known_hit"] = known.get(cid, "")
        with open(tsv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)
        verdicts = Counter((r["source"], r["verdict"]) for r in rows)
        logger.info("%d cluster(s) with >= %d reads: %s", len(rows), min_reads,
                    {f"{a}:{b}": n for (a, b), n in sorted(verdicts.items())})
        logger.info("Outputs: %s, %s", tsv_path, fasta_path)
        return rows
    finally:
        if not keep_temp:
            shutil.rmtree(work, ignore_errors=True)
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(prev_level)


def scan_known_tirs(fasta, flank=0, **kw):
    """Run find_tir on the ends of every record of a FASTA of known elements
    (e.g. test-data/IS_curated.fasta). With flank=0 the TIR must start at the
    annotated ends; a few bp of flank tolerates slightly trimmed annotations.
    Returns [(name, length, tir dict or None)]."""
    out = []
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        out.append((rec.id, len(seq), find_tir(seq, 0, len(seq), flank=max(flank, 0), **kw)))
    return out


def _cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Novel mobile element characterisation (stand-alone entry).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("-d", "--data_dir", required=True)
    r.add_argument("-g", "--genome", required=True)
    r.add_argument("-o", "--out_dir")
    r.add_argument("--min_reads", type=int, default=2)
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--flank_search", type=int, default=60)
    r.add_argument("--is_fasta")
    r.add_argument("--no_unannotated", action="store_true")
    r.add_argument("--max_consensus_reads", type=int, default=30)
    k = sub.add_parser("known")
    k.add_argument("fasta")
    k.add_argument("--flank", type=int, default=5)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if a.cmd == "run":
        main_novel_elements(a.data_dir, a.genome, a.out_dir, a.min_reads, a.threads, a.flank_search,
                            is_fasta=a.is_fasta, unannotated=not a.no_unannotated,
                            max_consensus_reads=a.max_consensus_reads)
    else:
        print("name\tlength\ttir_left\ttir_right\tleft_len\tright_len\tidentity\tscore\tleft_seq\tright_seq")
        for name, length, t in scan_known_tirs(a.fasta, a.flank):
            if t is None:
                print(f"{name}\t{length}\tnone")
            else:
                print(f"{name}\t{length}\t{t['left'][0]}-{t['left'][1]}\t{t['right'][0]}-{t['right'][1]}\t"
                      f"{t['left_len']}\t{t['right_len']}\t{t['identity']}\t{t['score']}\t"
                      f"{t['left_seq']}\t{t['right_seq']}")


if __name__ == "__main__":
    _cli()
