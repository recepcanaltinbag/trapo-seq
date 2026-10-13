#!/usr/bin/env python3
"""Per-insertion artefact evidence and flags (docs/ARTEFACTS.md).

A genuine simple insertion has plasmid on BOTH sides of the insert, on the
SAME strand, and the two plasmid flanks are ADJACENT on the circular plasmid
(they overlap by the target-site duplication, TSD). It is seen in independent
molecules. ONT chimeras, ligation products at the restriction site, fold-back
reads, internal adapters and plasmid multimers violate at least one of these.
This module measures, for every call of insertion_from_bam.tab, the evidence
and writes ``<tab stem>_artefacts.tsv`` next to it (the 8-column tab itself is
unchanged); it also infers the linearisation (restriction) site from the read
ends and writes ``<tab stem>_cut_site.tsv``.

Evidence per call (columns of the artefacts table, ARTEFACT_COLUMNS)
--------------------------------------------------------------------
* ``flanks``: both / left_only / right_only / none - plasmid alignment pieces
  ending at the insert start / starting at the insert end on the read
  (read-left and read-right flank, original read orientation), within
  FLANK_TOLERANCE bp. ``left_*`` / ``right_*``: the flank's plasmid interval
  (0-based half-open), strand, length on the read and MAPQ.
* ``plasmid_left_end`` (A) and ``plasmid_right_start`` (B): the flanks in
  PLASMID orientation - A is the exclusive end of the plasmid sequence before
  the insert, B the first plasmid base after it (for a minus-strand read the
  read-right flank is the plasmid-left one).
* ``flank_gap`` = A - B, signed and circular (in (-L/2, L/2]): the same sign
  convention as dr_finder's Repeat Length. A genuine insertion has
  flank_gap = TSD length (0-14 bp, plus chance extension of a flank into the
  insert); negative = plasmid missing between the flanks.
* ``strand_consistent``: both flanks on the same plasmid strand.
* ``min_flank_mapq``.
* ``junction``: A, else B (one-sided calls: the coordinate of the one flank);
  ``dist_to_cut``: circular distance from A / B to the nearest molecule end
  of the linearised plasmid (cut site, see below).
* ``reaches_read_end``: the insert starts/ends within READ_END_TOLERANCE bp
  of a read end.
* ``adapter_hits``: adapter/barcode motifs found in insert +- ADAPTER_FLANK bp
  (not within ADAPTER_END_MARGIN bp of a read end, where adapters are
  expected), as name:read_pos:edit_distance.
* ``foldback_kmers``: number of k-mers (k = FOLDBACK_K) of the read-left
  window (FOLDBACK_WINDOW bp before the insert) whose reverse complement occurs
  on one diagonal band of the read-right window: the read folds back on
  itself across the call (read = X + junk + revcomp(X)).
* ``plasmid_cover``: fraction of the insert interior (insert minus
  PLASMID_SELF_EDGE bp at each end) covered by k-mers (k = PLASMID_SELF_K)
  that occur in the plasmid (both strands, circular). Needs the plasmid FASTA.
* ``read_strand``, ``read_plasmid_start``, ``read_plasmid_end``: strand of the
  read relative to the plasmid and the plasmid coordinates of the molecule's
  two ends (first/last aligned base in plasmid orientation), used for the
  per-site independence check in g_blast_annot.
* ``flags``: comma-separated flags (or "-"), artefact flags first.

Flags
-----
ARTEFACT flags (evidence that the call is not a transposition into the
plasmid; g_blast_annot --artefact_action exclude relabels these calls
Contamination; see ARTEFACT_FLAGS):

* ``nonadjacent_flanks`` (C1): both flanks on one strand, but flank_gap >
  MAX_TSD_OVERLAP (40) bp or < MIN_FLANK_GAP (-100) bp. Known TSDs are 0-14
  bp, but the aligner extends a flank into the element when the element end
  resembles the target (D2) or through sequencing errors: on simulated true
  insertions (uniform 2 % substitutions and ONT r9.4.1 errors, 2 x 1000
  reads each, calibration seeds) flank_gap ranged from -47 to +30 bp (+30 = IS186B, 11 bp TSD plus
  19 bp extension). Errors at a junction make the aligner stop a flank early
  (negative values). A chimera joins two unrelated plasmid positions, so its
  flank_gap is uniform over the plasmid: the window [-100, +40] keeps
  ~2 % of 7 kb-plasmid chimeras. Real adjacent deletions (D6, IS26-type
  footprints) larger than 100 bp are also flagged: check them by hand.
  For IN calls (CIGAR insertion) the aligner puts the TSD copy inside the
  insertion; the overlap is recovered from the insert ends (in_call_overlap,
  exact matches, so it is a lower bound on noisy reads).
* ``strand_flip`` (C1/C3): the two flanks lie on opposite plasmid strands
  (inversions next to an insertion are possible biology but rare).
* ``unaligned_flank`` (C1): the alignment identity of a flank over its
  FLANK_CHECK_LEN (100) aligned read bases next to the insert (matches /
  columns, indels <= 50 bp counted) is < FLANK_MIN_IDENTITY (0.65): the aligner
  bridged a chimera junction by aligning unrelated DNA to a short plasmid
  stretch (seen with random foreign DNA between two plasmid pieces: several
  IN calls with "adjacent" flanks). Needs the plasmid FASTA. Evidence columns
  ``left_flank_identity`` / ``right_flank_identity``. Simulated true
  insertions: >= 0.89 (uniform errors) and >= 0.75 (ONT r9.4.1); bridged
  junctions: 0.42-0.61. Also raised (without the plasmid FASTA too) for a
  flank piece shorter than SANDWICH_MAX_FLANK (50) bp with > 50 bp of
  unaligned sequence on its other side as well (not a read end, not another
  piece, e.g. at the origin split of the circular mapping): a "flank" between two
  large insertions is a bridged junction, not evidence of a plasmid flank
  (a short stub at a read end - B1 - is not flagged).
* ``cut_site_junction`` (B2): a two-sided call with a junction within
  ``|end_offset| + CUT_SITE_SLACK`` (20) bp of the cut site (simulated
  ligation junctions lay up to 18 bp from it on ONT r9.4.1 reads). An insertion inside
  the recognition site destroys it (the molecule would not be linearised
  there) and one just next to it leaves a plasmid stub too short to align, so
  a two-sided call AT the cut site is a ligation of molecule ends.
* ``overhang_dr`` (B2): in addition, flank_gap equals the end signature of the
  enzyme (+-1 bp): the molecule-end offset (+k for a 5' overhang filled in
  by end-repair, -k for a 3' overhang removed by it, 0 blunt), or the length
  of the recognition site (sticky-end ligation of two fragments cut by the
  same enzyme recreates the site on both sides).
* ``one_sided_cut_site`` (B1/B2): one-sided call whose flank ends at the cut
  site and whose insert reaches the read end: foreign/genomic DNA ligated to
  the end of the linear plasmid (B2) or an insertion next to the cut site
  (B1, plasmid stub too short to align). g_blast_annot treats it as an
  artefact unless the insert is IS_DB with >= B1_EXPLAINED (90) % explained
  (B1: then the clip is the element plus a short stub).
* ``internal_adapter`` (C2): an adapter motif inside the read at the call.
* ``foldback`` (C3): foldback_kmers >= FOLDBACK_MIN_KMERS (10), spread over
  >= FOLDBACK_MIN_SPAN (100) bp. Random sequence gives ~0.1 shared 11-mers in
  two 500 bp windows; a fold-back read with 8 % errors per copy ~80.
* ``plasmid_self`` (A4): plasmid_cover >= PLASMID_SELF_MIN_COVER (0.5): the
  insert is plasmid sequence (dimer/concatemer junction or a misaligned
  plasmid copy).

CAUTION flags (context, never relabel):

* ``one_sided``: one-sided call away from the cut site (fragmented molecule
  ending inside an insert, or a chimera; g_blast_annot's contamination rule
  already handles one-sided calls without IS/genome hit).
* ``cut_in_insert`` (B3): one-sided call whose insert reaches the read end
  and whose terminal RECOGNITION_SLACK + ADAPTER_END_MARGIN (160) read bases
  contain a recognition site of the enzyme: the enzyme cut the element, which
  explains the missing flank (needs the enzyme, given or guessed). A 6 bp site
  occurs by chance in ~4 % of 160 bp windows, so this is a hint, not proof.
* ``foldback_read``: the read aligns to the same plasmid region on both
  strands (>= FOLDBACK_READ_OVERLAP bp): duplex/fold-back read, its calls are
  one molecule seen twice.
* ``low_mapq`` (D1): a flank with MAPQ < LOW_MAPQ (10): repeat in the plasmid.
* ``single_molecule_support`` and ``possible_bleed_through`` are site-level
  and added by g_blast_annot (see there).

Cut site (ARTEFACTS B)
----------------------
Reads of unfragmented plasmid molecules (aligned plasmid >= CUT_MIN_COVER of
the plasmid, <= CUT_MAX_END_CLIP bp unaligned at each read end) start and end
at the linearisation site. For each such read the plasmid coordinates of the
molecule start (first base, plasmid orientation) and end (exclusive) are
collected; the start and end with most reads within +-CUT_WINDOW bp are the
inferred molecule start (``cut_site`` = first base of the linear molecule)
and end. ``end_offset`` = end - start: +k for a 5' overhang of k bases filled
in by end-repair, -k for a 3' overhang removed by it, 0 for a blunt cut.
A site is reported when >= CUT_MIN_SUPPORT of the eligible reads (and >= 5
reads) start within +-CUT_WINDOW bp of it; otherwise (random fragmentation,
transposase kits) no cut site is used. With the plasmid sequence the enzyme
whose recognition site explains start and end (+-CUT_GUESS_TOL = 3 bp each)
is reported (``enzyme_guess``). The observed ends can lie a few bp outside the
enzyme's cut: with untrimmed adapters/barcodes the alignment extends into the
adapter by chance-matching bases (VD2 data: XbaI start 3239 / filled-in end
3243 observed as 3237 / 3246: reads start ...cagcac|CTCTAGAGG and end
...TCTAG|AGGtgctg..., the native-barcode flank CAGCACCT and its reverse
complement AGGTGCTG contribute CT / AGG that continue the plasmid). Both the
observed and the enzyme offset count as end signatures for ``overhang_dr``.
User input overrides: ``enzyme`` (ENZYMES table) or
``cut_site`` (0-based first base of the linear molecule) with ``overhang``
(5' overhang sequence; negative offset for 3' overhangs is given with a
leading '-', e.g. '-TGCA').

Adapters
--------
ADAPTERS holds the motif list; it is user-configurable (``adapter_fasta``
replaces it). Only the ONT ligation-kit Y-adapter top-strand prefix (as in
Porechop's adapter table) is built in; kit-specific barcode sequences are not,
because they differ between kits - supply them as a FASTA.

Command line:
    artefacts.py BAM TAB [--plasmid PLASMID.fasta] [--enzyme NAME |
        --cut_site POS --overhang SEQ] [--adapter_fasta FASTA]
"""

import argparse
import csv
import glob
import logging
import os
import sys
from collections import Counter, defaultdict

logger = logging.getLogger(__name__)

ARTEFACTS_SUFFIX = "_artefacts.tsv"
CUT_SITE_SUFFIX = "_cut_site.tsv"

# ------------------------------------------------------------- thresholds
FLANK_TOLERANCE = 20          # bp between an insert boundary and its flank piece on the read
READ_END_TOLERANCE = 20       # = g_blast_annot.READ_END_TOLERANCE
LARGE_INSERTION = 50          # CIGAR insertions longer than this split alignment pieces
MAX_TSD_OVERLAP = 40          # flank overlap above this: not a TSD
MIN_FLANK_GAP = -100          # flank gap below this: plasmid missing between the flanks
CUT_SITE_SLACK = 20           # bp added to |overhang| for the cut-site window
LOW_MAPQ = 10
ADAPTER_FLANK = 100           # bp on each side of the insert scanned for adapters
ADAPTER_END_MARGIN = 150      # adapter hits this close to a read end are ignored (expected there)
ADAPTER_SEED = 8
ADAPTER_MAX_ERROR = 0.18      # max edit distance / motif length (5 of 28)
FOLDBACK_K = 11
FOLDBACK_WINDOW = 500
FOLDBACK_MIN_WINDOW = 100
FOLDBACK_BAND = 30
FOLDBACK_MIN_KMERS = 10
FOLDBACK_MIN_SPAN = 100
FOLDBACK_READ_OVERLAP = 200
PLASMID_SELF_K = 13
PLASMID_SELF_EDGE = 100
PLASMID_SELF_MIN_LEN = 100    # interior bp needed to judge
PLASMID_SELF_MIN_COVER = 0.5
RECOGNITION_SLACK = 10
FLANK_CHECK_LEN = 100         # aligned read bp of a flank next to the insert checked
FLANK_MIN_IDENTITY = 0.65     # below this the "flank" is not plasmid sequence
SANDWICH_MAX_FLANK = 50       # a flank this short with unaligned sequence on both sides is no evidence
B1_EXPLAINED = 90.0           # used by g_blast_annot (one_sided_cut_site exemption)
# cut-site inference
CUT_MIN_COVER = 0.8
CUT_MAX_END_CLIP = 200
CUT_WINDOW = 3
CUT_MIN_SUPPORT = 0.2
CUT_MIN_READS = 5
CUT_GUESS_TOL = 3             # bp between observed molecule start/end and an enzyme's cut
PLASMID_FASTA_MAX_BYTES = 5_000_000

ARTEFACT_FLAGS = ("nonadjacent_flanks", "strand_flip", "unaligned_flank", "cut_site_junction",
                  "overhang_dr", "one_sided_cut_site", "internal_adapter", "foldback", "plasmid_self")
CAUTION_FLAGS = ("one_sided", "cut_in_insert", "foldback_read", "low_mapq",
                 "single_molecule_support", "possible_bleed_through")
ALL_FLAGS = ARTEFACT_FLAGS + CAUTION_FLAGS

# Restriction enzymes: recognition site (palindromic, so one strand is enough),
# top-strand cut and bottom-strand cut as offsets from the site start in
# top-strand coordinates. Overhang = site[top:bottom] (5') or site[bottom:top] (3').
ENZYMES = {
    "EcoRI": ("GAATTC", 1, 5), "BamHI": ("GGATCC", 1, 5), "HindIII": ("AAGCTT", 1, 5),
    "XbaI": ("TCTAGA", 1, 5), "SpeI": ("ACTAGT", 1, 5), "NotI": ("GCGGCCGC", 2, 6),
    "PstI": ("CTGCAG", 5, 1), "SalI": ("GTCGAC", 1, 5), "KpnI": ("GGTACC", 5, 1),
    "SacI": ("GAGCTC", 5, 1), "XhoI": ("CTCGAG", 1, 5), "NcoI": ("CCATGG", 1, 5),
    "NdeI": ("CATATG", 2, 4), "EcoRV": ("GATATC", 3, 3), "SmaI": ("CCCGGG", 3, 3),
    "BglII": ("AGATCT", 1, 5), "NheI": ("GCTAGC", 1, 5), "SphI": ("GCATGC", 5, 1),
    "ApaI": ("GGGCCC", 5, 1), "ClaI": ("ATCGAT", 2, 4),
}

# Adapter / barcode motifs searched inside reads (user-configurable: pass a
# FASTA with adapter_fasta to replace this list). Built in: the top-strand
# prefix of the ONT ligation-kit Y-adapter as listed in Porechop's adapter
# table. Barcodes and adapters of other kits (rapid, native barcoding, newer
# chemistries) differ - give them as a FASTA rather than relying on this list.
ADAPTERS = {
    "ONT_ligation_Y_top": "AATGTACTTCGTTCAGTTACGTATTGCT",
}

ARTEFACT_COLUMNS = [
    "Query ID", "read_id", "query_start", "query_end", "insert_type", "read_len",
    "flanks", "left_ref_start", "left_ref_end", "left_strand", "left_len", "left_mapq",
    "right_ref_start", "right_ref_end", "right_strand", "right_len", "right_mapq",
    "plasmid_left_end", "plasmid_right_start", "flank_gap", "strand_consistent",
    "min_flank_mapq", "junction", "cut_site", "dist_to_cut", "reaches_read_end",
    "left_flank_identity", "right_flank_identity",
    "adapter_hits", "foldback_kmers", "plasmid_cover", "read_strand",
    "read_plasmid_start", "read_plasmid_end", "flags",
]
CUT_SITE_KEYS = ["contig", "plasmid_len", "source", "cut_site", "cut_end", "end_offset",
                 "overhang", "enzyme", "enzyme_guess", "enzyme_offset", "recognition_sites", "eligible_reads",
                 "start_support", "end_support", "start_support_fraction", "end_support_fraction"]

M, I, D, N, S, H, P, EQ, X = range(9)
QUERY_OPS = {M, I, S, H, EQ, X}
REF_OPS = {M, D, N, EQ, X}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def circ_signed(d, L):
    """d modulo L in (-L/2, L/2]."""
    if not L:
        return d
    d %= L
    if d > L // 2:
        d -= L
    return d


def circ_dist(a, b, L):
    return abs(circ_signed(a - b, L))


# ------------------------------------------------------------------ FASTA

def read_fasta(path):
    """[(name, sequence upper-case)] in file order."""
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def find_plasmid_fasta(bam_path, contig, length):
    """A FASTA next to the BAM or in its parent folder (the data dir) whose
    first record has the BAM reference's name and length, or None."""
    bam_dir = os.path.dirname(os.path.abspath(bam_path))
    for d in (bam_dir, os.path.dirname(bam_dir)):
        cands = []
        for ext in ("*.fasta", "*.fa", "*.fna", "*.fas"):
            cands += glob.glob(os.path.join(glob.escape(d), ext))
        for path in sorted(cands):
            try:
                if os.path.getsize(path) > PLASMID_FASTA_MAX_BYTES:
                    continue
                recs = read_fasta(path)
            except (OSError, UnicodeDecodeError):
                continue
            for name, seq in recs[:5]:
                if name == contig and len(seq) == length:
                    return path
    return None


def read_adapters(adapter_fasta=None):
    if not adapter_fasta:
        return dict(ADAPTERS)
    recs = read_fasta(adapter_fasta)
    if not recs:
        raise ValueError(f"no sequences in adapter FASTA {adapter_fasta}")
    return {n: s for n, s in recs if s}


# --------------------------------------------------------------- enzymes

def enzyme_info(name):
    """dict(name, site, top, bottom, offset, overhang (signed text)) for a table enzyme."""
    key = {k.lower(): k for k in ENZYMES}.get(str(name).lower())
    if key is None:
        raise ValueError(f"unknown enzyme {name!r}; known: {', '.join(sorted(ENZYMES))} "
                         f"(or give --cut_site and --overhang)")
    site, top, bottom = ENZYMES[key]
    off = bottom - top
    if off > 0:
        oh = site[top:bottom]
    elif off < 0:
        oh = "-" + site[bottom:top]
    else:
        oh = ""
    return {"name": key, "site": site, "top": top, "bottom": bottom, "offset": off, "overhang": oh}


def parse_overhang(overhang):
    """'AATT' -> +4, '-TGCA' -> -4, '' / None -> 0."""
    if not overhang:
        return 0
    s = str(overhang).strip().upper()
    if s.startswith("-"):
        return -len(s[1:])
    return len(s)


def recognition_sites(plasmid_seq, site):
    """0-based start positions of ``site`` on the circular plasmid (+ strand;
    the table enzymes are palindromic)."""
    L = len(plasmid_seq)
    ext = plasmid_seq + plasmid_seq[:len(site) - 1]
    out, i = [], ext.find(site)
    while i != -1:
        if i < L:
            out.append(i)
        i = ext.find(site, i + 1)
    return out


# ---------------------------------------------------------------- reads

class Piece:
    """One aligned piece of a read (original read orientation)."""
    __slots__ = ("q_start", "q_end", "ref_start", "ref_end", "reverse", "mapq", "contig", "record")

    def __init__(self, q_start, q_end, ref_start, ref_end, reverse, mapq, contig, record=None):
        self.q_start, self.q_end = q_start, q_end
        self.ref_start, self.ref_end = ref_start, ref_end
        self.reverse, self.mapq, self.contig = reverse, mapq, contig
        self.record = record          # (reference_start, cigartuples, total read length)

    @property
    def ref_at_q_start(self):
        """Plasmid coordinate at the read-left end (exclusive end for a minus piece)."""
        return self.ref_end if self.reverse else self.ref_start

    @property
    def ref_at_q_end(self):
        return self.ref_start if self.reverse else self.ref_end


def _read_length(rec):
    return sum(n for op, n in rec.cigartuples if op in QUERY_OPS)


def record_pieces(rec, contig, max_ins=LARGE_INSERTION, keep_cigar=True):
    """Aligned pieces of a record, split at clips and CIGAR insertions > max_ins.
    keep_cigar: attach the record's CIGAR (needed for flank identity only)."""
    total = _read_length(rec)
    out, q, ref = [], 0, rec.reference_start
    q0 = r0 = None
    for op, n in rec.cigartuples:
        split = op in (S, H) or (op == I and n > max_ins)
        if split:
            if q0 is not None and q > q0:
                out.append((q0, q, r0, ref))
            q0 = None
        elif q0 is None and op in QUERY_OPS:
            q0, r0 = q, ref
        if op in QUERY_OPS:
            q += n
        if op in REF_OPS:
            ref += n
    if q0 is not None and q > q0:
        out.append((q0, q, r0, ref))
    pieces = []
    record = (rec.reference_start, tuple(rec.cigartuples), total) if keep_cigar else None
    for s, e, rs, re in out:
        if rec.is_reverse:
            s, e = total - e, total - s
        pieces.append(Piece(s, e, rs, re, rec.is_reverse, rec.mapping_quality, contig, record))
    return pieces


def window_identity(piece, lo, hi, read_seq, plasmid_seq):
    """Alignment identity (matches / aligned columns incl. indels) of the
    piece's record over the read interval [lo, hi) (original orientation)."""
    if piece.record is None or not read_seq or not plasmid_seq:
        return None
    ref, cigar, total = piece.record
    if piece.reverse:
        lo, hi = total - hi, total - lo
    L = len(plasmid_seq)
    q = 0
    match = cols = 0
    for op, n in cigar:
        if op in (M, EQ, X):
            a, b = max(q, lo), min(q + n, hi)
            for k in range(a, b):
                qi = k if not piece.reverse else total - 1 - k
                base = read_seq[qi] if not piece.reverse else read_seq[qi].translate(COMPLEMENT)
                cols += 1
                match += base == plasmid_seq[(ref + k - q) % L]
            q += n
            ref += n
        elif op == I:
            if lo <= q < hi and n <= LARGE_INSERTION:
                cols += min(n, hi - q)
            q += n
        elif op in (S, H):
            q += n
        elif op in (D, N):
            if lo < q < hi and n <= LARGE_INSERTION:
                cols += n
            ref += n
        if q >= hi:
            break
    return match / cols if cols >= 30 else None


def load_reads(bam_path, names=None, detail_names=None):
    """{read: dict(pieces, length, seq (original orientation or None))} from
    primary + supplementary records (secondary ignored). Also returns
    {contig: length}. Only reads in ``names`` (default all) are loaded; the
    sequence and CIGARs are kept only for reads in ``detail_names`` (default:
    all loaded reads) so that a large BAM does not have to fit in memory."""
    import pysam
    reads = {}
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or not rec.cigartuples:
                continue
            if names is not None and rec.query_name not in names:
                continue
            detail = detail_names is None or rec.query_name in detail_names
            info = reads.setdefault(rec.query_name, {"pieces": [], "length": 0, "seq": None})
            info["pieces"].extend(record_pieces(rec, rec.reference_name, keep_cigar=detail))
            total = _read_length(rec)
            info["length"] = max(info["length"], total)
            if (detail and not rec.is_supplementary and rec.query_sequence
                    and len(rec.query_sequence) == total):
                seq = rec.query_sequence.upper()
                info["seq"] = revcomp(seq) if rec.is_reverse else seq
    for info in reads.values():
        info["pieces"].sort(key=lambda p: (p.q_start, p.q_end))
    return reads, lengths


def molecule_ends(info, L):
    """(start, end, strand, clip_start, clip_end) of a read in plasmid
    orientation: start = plasmid coordinate of the molecule's first base,
    end = exclusive end; strand '+' if the read starts with a plus piece."""
    pcs = info["pieces"]
    if not pcs:
        return None
    first, last = pcs[0], max(pcs, key=lambda p: p.q_end)
    if not first.reverse:
        start = first.ref_start
        end = last.ref_end if not last.reverse else None
    else:
        end = first.ref_end
        start = last.ref_start if last.reverse else None
    if start is None or end is None:        # read starts and ends on different strands
        return None
    strand = "-" if first.reverse else "+"
    return start % L, end % L if L else end, strand, first.q_start, info["length"] - last.q_end


def covered_ref(pieces, L):
    """Plasmid bases covered by the pieces (union, circular)."""
    ivs = sorted((p.ref_start, p.ref_end) for p in pieces)
    tot, cs, ce = 0, None, None
    for s, e in ivs:
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if ce is not None:
        tot += ce - cs
    return min(tot, L) if L else tot


# --------------------------------------------------------- cut-site inference

def _circular_mode(values, L, window=CUT_WINDOW):
    """(position, count within +-window) of the densest point (ties: lowest position)."""
    if not values:
        return None, 0
    cnt = Counter(v % L for v in values)
    keys = sorted(cnt)
    best, best_n = None, -1
    for k in keys:
        n = sum(cnt.get((k + d) % L, 0) for d in range(-window, window + 1))
        if n > best_n or (n == best_n and cnt[k] > cnt.get(best, 0)):
            best, best_n = k, n
    # refine: the most common exact value within the window
    sub = [((best + d) % L, cnt.get((best + d) % L, 0)) for d in range(-window, window + 1)]
    best = max(sub, key=lambda x: (x[1], -abs(circ_signed(x[0] - best, L))))[0]
    n = sum(cnt.get((best + d) % L, 0) for d in range(-window, window + 1))
    return best, n


def infer_cut_site(reads, lengths):
    """Infer the dominant linearisation site from read ends (see module docstring).
    Returns a dict (CUT_SITE_KEYS subset) or None if no reads."""
    by_contig = defaultdict(list)
    for info in reads.values():
        if not info["pieces"]:
            continue
        contig = Counter(p.contig for p in info["pieces"]).most_common(1)[0][0]
        by_contig[contig].append(info)
    if not by_contig:
        return None
    contig = max(by_contig, key=lambda c: len(by_contig[c]))
    L = lengths.get(contig)
    starts, ends = [], []
    for info in by_contig[contig]:
        pcs = [p for p in info["pieces"] if p.contig == contig]
        if covered_ref(pcs, L) < CUT_MIN_COVER * L:
            continue
        me = molecule_ends(dict(info, pieces=pcs), L)
        if me is None:
            continue
        start, end, _, c0, c1 = me
        if c0 > CUT_MAX_END_CLIP or c1 > CUT_MAX_END_CLIP:
            continue
        starts.append(start)
        ends.append(end)
    n = len(starts)
    res = {"contig": contig, "plasmid_len": L, "eligible_reads": n, "source": "none",
           "cut_site": None, "cut_end": None, "end_offset": None, "start_support": 0,
           "end_support": 0, "start_support_fraction": None, "end_support_fraction": None}
    if not n:
        return res
    s, ns = _circular_mode(starts, L)
    e, ne = _circular_mode(ends, L)
    res.update(start_support=ns, end_support=ne, start_support_fraction=round(ns / n, 3),
               end_support_fraction=round(ne / n, 3))
    if ns >= CUT_MIN_READS and ns / n >= CUT_MIN_SUPPORT:
        res.update(source="inferred", cut_site=s)
        if ne >= CUT_MIN_READS and ne / n >= CUT_MIN_SUPPORT:
            res.update(cut_end=e, end_offset=circ_signed(e - s, L))
        else:
            res.update(cut_end=s, end_offset=0)
    return res


def resolve_cut_site(inferred, plasmid_seq=None, enzyme=None, cut_site=None, overhang=None):
    """Combine the inferred site with user input / the enzyme table.
    Returns a dict with cut_site, cut_end, end_offset, overhang, enzyme,
    enzyme_guess, recognition_sites, source (+ the inference statistics)."""
    res = dict(inferred or {})
    L = res.get("plasmid_len") or (len(plasmid_seq) if plasmid_seq else None)
    res.setdefault("plasmid_len", L)
    res.update(enzyme="", enzyme_guess="", recognition_sites="", overhang="")
    site_len = None
    if enzyme:
        info = enzyme_info(enzyme)
        res["enzyme"] = info["name"]
        res["overhang"] = info["overhang"]
        res["enzyme_offset"] = info["offset"]
        site_len = len(info["site"])
        if plasmid_seq:
            sites = recognition_sites(plasmid_seq, info["site"])
            res["recognition_sites"] = ",".join(map(str, sites))
            if not sites:
                logger.warning("enzyme %s has no recognition site in the plasmid", info["name"])
            else:
                starts = [(s + info["top"]) % L for s in sites]
                if res.get("cut_site") is not None:
                    best = min(starts, key=lambda x: circ_dist(x, res["cut_site"], L))
                    if circ_dist(best, res["cut_site"], L) > 20:
                        logger.warning("enzyme %s: no recognition site within 20 bp of the inferred "
                                       "cut site %s; using its site at %d", info["name"],
                                       res["cut_site"], best)
                elif len(starts) > 1:
                    logger.warning("enzyme %s cuts the plasmid %d times (%s) and no dominant read-end "
                                   "pile-up was found; using the first site", info["name"], len(starts),
                                   res["recognition_sites"])
                    best = starts[0]
                else:
                    best = starts[0]
                res.update(cut_site=best, cut_end=(best + info["offset"]) % L,
                           end_offset=info["offset"], source="enzyme")
        else:
            res["end_offset"] = info["offset"]
            if res.get("cut_site") is not None:
                res["cut_end"] = (res["cut_site"] + info["offset"]) % L
                res["source"] = "inferred+enzyme"
    if cut_site is not None:
        off = parse_overhang(overhang)
        res.update(cut_site=int(cut_site) % L, end_offset=off, cut_end=(int(cut_site) + off) % L,
                   overhang=(overhang or "").upper(), source="user")
    if not enzyme and cut_site is None and plasmid_seq and res.get("cut_site") is not None:
        # which table enzyme explains the observed molecule start and end?
        guesses = []
        for name in sorted(ENZYMES):
            info = enzyme_info(name)
            for s in recognition_sites(plasmid_seq, info["site"]):
                d0 = circ_dist((s + info["top"]) % L, res["cut_site"], L)
                d1 = circ_dist((s + info["bottom"]) % L, res["cut_end"], L)
                if d0 <= CUT_GUESS_TOL and d1 <= CUT_GUESS_TOL:
                    guesses.append((d0 + d1, name, s, info))
        if guesses:
            guesses.sort()
            _, name, s, info = guesses[0]
            res["enzyme_guess"] = ",".join(sorted({g[1] for g in guesses}))
            res["overhang"] = info["overhang"]
            res["enzyme_offset"] = info["offset"]
            res["recognition_sites"] = ",".join(map(str, recognition_sites(plasmid_seq, info["site"])))
            site_len = len(info["site"])
    res["site_len"] = site_len
    return res


# --------------------------------------------------------------- per call

def call_flanks(pieces, qs, qe, tol=FLANK_TOLERANCE):
    """(read-left flank piece, read-right flank piece), either may be None."""
    left = [p for p in pieces if abs(p.q_end - qs) <= tol and p.q_start < qs]
    right = [p for p in pieces if abs(p.q_start - qe) <= tol and p.q_end > qe]
    lp = min(left, key=lambda p: (abs(p.q_end - qs), -(p.q_end - p.q_start))) if left else None
    rp = min(right, key=lambda p: (abs(p.q_start - qe), -(p.q_end - p.q_start))) if right else None
    return lp, rp


def plasmid_junctions(lp, rp, L):
    """(A, B): exclusive end of the plasmid-left flank and start of the
    plasmid-right flank (None where missing)."""
    A = B = None
    for p, side in ((lp, "L"), (rp, "R")):
        if p is None:
            continue
        # + strand: read-left flank is plasmid-left; - strand: swapped
        plasmid_left = (side == "L") != p.reverse
        if plasmid_left:
            A = p.ref_end % L
        else:
            B = p.ref_start % L
    return A, B


def _edit_distance_semiglobal(pattern, text):
    """Min edit distance of ``pattern`` against any substring of ``text`` and the
    end position (exclusive) in text."""
    m = len(pattern)
    prev = list(range(m + 1))
    best, best_end = prev[m], 0
    for j, tc in enumerate(text, 1):
        cur = [0] * (m + 1)
        for i in range(1, m + 1):
            cost = 0 if pattern[i - 1] == tc else 1
            cur[i] = min(prev[i - 1] + cost, prev[i] + 1, cur[i - 1] + 1)
        if cur[m] < best:
            best, best_end = cur[m], j
        prev = cur
    return best, best_end


class AdapterScanner:
    """Seed-and-verify search of adapter motifs (both strands)."""

    def __init__(self, adapters, k=ADAPTER_SEED, max_error=ADAPTER_MAX_ERROR):
        self.k, self.max_error = k, max_error
        self.motifs = []
        self.seeds = defaultdict(list)
        for name, seq in adapters.items():
            seq = seq.upper()
            if len(seq) < k:
                continue
            for strand, s in (("+", seq), ("-", revcomp(seq))):
                idx = len(self.motifs)
                self.motifs.append((name, strand, s))
                for o in range(len(s) - k + 1):
                    self.seeds[s[o:o + k]].append((idx, o))

    def scan(self, seq, lo, hi):
        """Hits [(name, read_pos, edit_distance)] in seq[lo:hi]."""
        if not self.motifs or hi - lo < self.k:
            return []
        k, seeds = self.k, self.seeds
        cand = Counter()
        for p in range(lo, hi - k + 1):
            for idx, o in seeds.get(seq[p:p + k], ()):
                cand[(idx, (p - o) // 4)] += 1
        hits, seen = [], []
        for (idx, b), n in cand.items():
            if n < 2:
                continue
            name, strand, motif = self.motifs[idx]
            m = len(motif)
            a = max(lo, b * 4 - m // 4 - 4)
            z = min(hi, b * 4 + m + m // 4 + 4)
            dist, end = _edit_distance_semiglobal(motif, seq[a:z])
            if dist <= self.max_error * m:
                pos = a + max(0, end - m)
                if all(not (nm == name and abs(pp - pos) < m) for nm, pp, _ in seen):
                    seen.append((name, pos, dist))
        hits.extend(sorted(seen, key=lambda h: h[1]))
        return hits


def foldback_score(seq, qs, qe, window=FOLDBACK_WINDOW, k=FOLDBACK_K, band=FOLDBACK_BAND):
    """(matches, span): k-mers of the read-left window whose reverse complement
    lies in the read-right window on one diagonal band (+-band), and the read
    span they cover. None if a window is shorter than FOLDBACK_MIN_WINDOW."""
    left = seq[max(0, qs - window):qs]
    right = seq[qe:qe + window]
    if len(left) < FOLDBACK_MIN_WINDOW or len(right) < FOLDBACK_MIN_WINDOW:
        return None
    rc = revcomp(right)
    index = defaultdict(list)
    for j in range(len(rc) - k + 1):
        index[rc[j:j + k]].append(j)
    diag = defaultdict(list)
    for i in range(len(left) - k + 1):
        for j in index.get(left[i:i + k], ()):
            diag[(i - j) // band].append(i)
    best, span = 0, 0
    for d in diag:
        pts = diag[d] + diag.get(d + 1, [])
        if len(pts) > best:
            best, span = len(pts), max(pts) - min(pts) + k
    return best, span


class PlasmidKmers:
    def __init__(self, plasmid_seq, k=PLASMID_SELF_K):
        self.k = k
        ext = plasmid_seq + plasmid_seq[:k - 1]
        kms = {ext[i:i + k] for i in range(len(plasmid_seq))}
        self.kmers = kms | {revcomp(x) for x in kms}

    def cover(self, insert, edge=PLASMID_SELF_EDGE, min_len=PLASMID_SELF_MIN_LEN):
        """Fraction of the insert interior covered by plasmid k-mers (None if too short)."""
        core = insert[edge:len(insert) - edge] if len(insert) > 2 * edge else ""
        if len(core) < min_len:
            return None
        k, kms = self.k, self.kmers
        covered = [False] * len(core)
        last = -1
        for i in range(len(core) - k + 1):
            if core[i:i + k] in kms:
                for j in range(max(i, last + 1), i + k):
                    covered[j] = True
                last = i + k - 1
        return sum(covered) / len(core)


def read_foldback(pieces, L, min_overlap=FOLDBACK_READ_OVERLAP):
    """The read aligns to the same plasmid region on both strands."""
    plus = [(p.ref_start, p.ref_end) for p in pieces if not p.reverse]
    minus = [(p.ref_start, p.ref_end) for p in pieces if p.reverse]
    for a0, a1 in plus:
        for b0, b1 in minus:
            if min(a1, b1) - max(a0, b0) >= min_overlap:
                return True
    return False


def in_call_overlap(ins_p, plasmid_seq, X, max_ext=60):
    """IN calls: the aligner puts the TSD copy inside the CIGAR insertion, so
    both flanks end at the same coordinate X. The flank overlap is recovered as
    the number of insert bases that continue the plasmid on either side:
    (bases of ins_p's start equal to plasmid[X:], bases of its end equal to
    plasmid[:X]). ins_p: insert in plasmid orientation."""
    L = len(plasmid_seq)
    ext = plasmid_seq + plasmid_seq
    X %= L
    k_left = 0
    while k_left < min(max_ext, len(ins_p)) and ins_p[k_left] == ext[X + k_left]:
        k_left += 1
    k_right = 0
    while k_right < min(max_ext, len(ins_p)) and ins_p[-1 - k_right] == ext[X + L - 1 - k_right]:
        k_right += 1
    return k_left, k_right


def evaluate_call(call, info, cut, L, scanner=None, pk=None, recog=None, plasmid_seq=None):
    """Evidence dict (ARTEFACT_COLUMNS) and flags for one insertion call.
    call: (read_id, query_start, query_end, insert_type)."""
    read_id, qs, qe, itype = call
    row = {c: "" for c in ARTEFACT_COLUMNS}
    row.update({"Query ID": f"{read_id}_{qs}_{qe}", "read_id": read_id, "query_start": qs,
                "query_end": qe, "insert_type": itype})
    flags = []
    if info is None:
        row.update(flanks="none", flags="-")
        return row
    pieces, rlen, seq = info["pieces"], info["length"], info["seq"]
    row["read_len"] = rlen
    lp, rp = call_flanks(pieces, qs, qe)
    row["flanks"] = "both" if lp and rp else "left_only" if lp else "right_only" if rp else "none"
    mapqs = []
    for p, side in ((lp, "left"), (rp, "right")):
        if p is None:
            continue
        row.update({f"{side}_ref_start": p.ref_start, f"{side}_ref_end": p.ref_end,
                    f"{side}_strand": "-" if p.reverse else "+", f"{side}_len": p.q_end - p.q_start,
                    f"{side}_mapq": p.mapq})
        mapqs.append(p.mapq)
    A, B = plasmid_junctions(lp, rp, L)
    row["plasmid_left_end"] = "" if A is None else A
    row["plasmid_right_start"] = "" if B is None else B
    if mapqs:
        row["min_flank_mapq"] = min(mapqs)
        if min(mapqs) < LOW_MAPQ:
            flags.append("low_mapq")
    reaches_end = qs <= READ_END_TOLERANCE or qe >= rlen - READ_END_TOLERANCE
    row["reaches_read_end"] = reaches_end
    two_sided = lp is not None and rp is not None
    if (two_sided and itype == "IN" and plasmid_seq and seq and lp.reverse == rp.reverse
            and A is not None and B is not None and circ_dist(A, B, L) <= 1):
        ins_p = seq[qs:qe] if not lp.reverse else revcomp(seq[qs:qe])
        k_left, k_right = in_call_overlap(ins_p, plasmid_seq, A)
        A, B = (A + k_left) % L, (B - k_right) % L
        row["plasmid_left_end"], row["plasmid_right_start"] = A, B
    if two_sided:
        consistent = lp.reverse == rp.reverse
        row["strand_consistent"] = consistent
        if consistent:
            gap = circ_signed(A - B, L)
            row["flank_gap"] = gap
            if gap > MAX_TSD_OVERLAP or gap < MIN_FLANK_GAP:
                flags.append("nonadjacent_flanks")
        else:
            flags.append("strand_flip")
    junctions = [x for x in (A, B) if x is not None]
    row["junction"] = junctions[0] if junctions else ""
    # cut site
    if cut and cut.get("cut_site") is not None and junctions:
        ends = [cut["cut_site"], cut["cut_end"] if cut.get("cut_end") is not None else cut["cut_site"]]
        d = min(circ_dist(j, e, L) for j in junctions for e in ends)
        row["cut_site"] = cut["cut_site"]
        row["dist_to_cut"] = d
        win = abs(cut.get("end_offset") or 0) + CUT_SITE_SLACK
        if d <= win:
            if two_sided:
                flags.append("cut_site_junction")
                gap = row["flank_gap"]
                sigs = {cut.get("end_offset") or 0}
                if cut.get("enzyme_offset") is not None:
                    sigs.add(cut["enzyme_offset"])
                if cut.get("site_len"):
                    sigs.add(cut["site_len"])
                if gap != "" and any(abs(gap - s) <= 1 for s in sigs):
                    flags.append("overhang_dr")
            elif reaches_end:
                flags.append("one_sided_cut_site")
    if not two_sided and (lp is not None or rp is not None):
        if "one_sided_cut_site" not in flags:
            flags.append("one_sided")
            if recog and seq and reaches_end:
                # B3: the read end inside the insert lies at a recognition site
                end = 0 if qs <= READ_END_TOLERANCE else rlen
                for site in recog:
                    # untrimmed adapter/barcode may precede the molecule end
                    span = len(site) + RECOGNITION_SLACK + ADAPTER_END_MARGIN
                    if end == 0:
                        hit = seq[:span].find(site) != -1
                    else:
                        hit = seq[max(0, rlen - span):].find(site) != -1
                    if hit:
                        flags.append("cut_in_insert")
                        break
    # sequence-based evidence
    if seq:
        if scanner is not None:
            lo = max(ADAPTER_END_MARGIN, qs - ADAPTER_FLANK)
            hi = min(rlen - ADAPTER_END_MARGIN, qe + ADAPTER_FLANK)
            hits = scanner.scan(seq, lo, hi) if hi > lo else []
            row["adapter_hits"] = ";".join(f"{n}:{p}:{d}" for n, p, d in hits) or "-"
            if hits:
                flags.append("internal_adapter")
        fb = foldback_score(seq, qs, qe)
        if fb is not None:
            row["foldback_kmers"] = fb[0]
            if fb[0] >= FOLDBACK_MIN_KMERS and fb[1] >= FOLDBACK_MIN_SPAN:
                flags.append("foldback")
        if plasmid_seq:
            weak = []
            for side, piece in (("left", lp), ("right", rp)):
                if piece is None:
                    continue
                if side == "left":
                    lo, hi = max(piece.q_start, piece.q_end - FLANK_CHECK_LEN), piece.q_end
                else:
                    lo, hi = piece.q_start, min(piece.q_end, piece.q_start + FLANK_CHECK_LEN)
                ident = window_identity(piece, lo, hi, seq, plasmid_seq)
                if ident is not None:
                    row[f"{side}_flank_identity"] = round(ident, 3)
                    weak.append(ident < FLANK_MIN_IDENTITY)
            if any(weak):
                flags.append("unaligned_flank")
        if pk is not None:
            cov = pk.cover(seq[qs:qe])
            if cov is not None:
                row["plasmid_cover"] = round(cov, 3)
                if cov >= PLASMID_SELF_MIN_COVER:
                    flags.append("plasmid_self")
    for side, piece in (("left", lp), ("right", rp)):
        if piece is None or piece.q_end - piece.q_start >= SANDWICH_MAX_FLANK:
            continue
        if side == "left":
            far = piece.q_start
            inside = far > READ_END_TOLERANCE
            abutting = any(o is not piece and abs(o.q_end - far) <= LARGE_INSERTION for o in pieces)
        else:
            far = piece.q_end
            inside = far < rlen - READ_END_TOLERANCE
            abutting = any(o is not piece and abs(o.q_start - far) <= LARGE_INSERTION for o in pieces)
        if inside and not abutting and "unaligned_flank" not in flags:
            flags.append("unaligned_flank")
    if read_foldback(pieces, L) and "foldback" not in flags:
        flags.append("foldback_read")
    me = molecule_ends(info, L)
    if me is not None:
        row.update(read_plasmid_start=me[0], read_plasmid_end=me[1], read_strand=me[2])
    elif pieces:
        row["read_strand"] = Counter("-" if p.reverse else "+" for p in pieces).most_common(1)[0][0]
    flags = [f for f in ALL_FLAGS if f in flags]
    row["flags"] = ",".join(flags) or "-"
    return row


def split_flags(value):
    """'a,b' -> (artefact flags, caution flags); '-'/'' -> ([], [])."""
    fl = [f for f in str(value or "").replace(";", ",").split(",") if f and f != "-"]
    return [f for f in fl if f in ARTEFACT_FLAGS], [f for f in fl if f not in ARTEFACT_FLAGS]


# ------------------------------------------------------------------ main

def artefacts_path_for(tab_file):
    return os.path.splitext(tab_file)[0] + ARTEFACTS_SUFFIX


def cut_site_path_for(tab_file):
    return os.path.splitext(tab_file)[0] + CUT_SITE_SUFFIX


def _atomic_write(path, write):
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)), f".tmp_{os.getpid()}_{os.path.basename(path)}")
    try:
        with open(tmp, "w", newline="") as fh:
            write(fh)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def write_rows(rows, path, columns=ARTEFACT_COLUMNS):
    def w(fh):
        wr = csv.DictWriter(fh, columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    _atomic_write(path, w)


def write_cut_site(cut, path):
    def w(fh):
        fh.write("key\tvalue\n")
        for k in CUT_SITE_KEYS:
            v = cut.get(k)
            fh.write(f"{k}\t{'' if v is None else v}\n")
    _atomic_write(path, w)


def read_artefacts(path):
    """{Query ID: row dict} of an artefacts table."""
    with open(path, newline="") as fh:
        return {r["Query ID"]: r for r in csv.DictReader(fh, delimiter="\t")}


def flag_calls(insertions, bamfile, plasmid_fasta=None, enzyme=None, cut_site=None, overhang=None,
               adapter_fasta=None):
    """insertions: iterable of (read_id, query_start, query_end, insert_type).
    Returns (rows, cut_site_dict)."""
    insertions = list(insertions)
    reads, lengths = load_reads(bamfile, detail_names={c[0] for c in insertions})
    cut = infer_cut_site(reads, lengths) or {}
    contig = cut.get("contig") or (next(iter(lengths)) if lengths else None)
    L = lengths.get(contig) if contig else None
    plasmid_seq = None
    if plasmid_fasta is None and contig is not None:
        plasmid_fasta = find_plasmid_fasta(bamfile, contig, L)
    if plasmid_fasta:
        recs = dict(read_fasta(plasmid_fasta))
        plasmid_seq = recs.get(contig) or (next(iter(recs.values())) if recs else None)
        if plasmid_seq and L and len(plasmid_seq) != L:
            logger.warning("%s: plasmid length %d differs from the BAM reference (%d); sequence "
                           "checks skipped", plasmid_fasta, len(plasmid_seq), L)
            plasmid_seq = None
    cut = resolve_cut_site(cut, plasmid_seq, enzyme, cut_site, overhang)
    cut["plasmid_fasta"] = plasmid_fasta or ""
    scanner = AdapterScanner(read_adapters(adapter_fasta))
    pk = PlasmidKmers(plasmid_seq) if plasmid_seq else None
    recog = None
    name = cut.get("enzyme") or (cut.get("enzyme_guess") or "").split(",")[0]
    if name:
        recog = [enzyme_info(name)["site"]]
    rows = [evaluate_call(c, reads.get(c[0]), cut, L, scanner, pk, recog, plasmid_seq) for c in insertions]
    return rows, cut


def describe_cut(cut):
    if cut.get("cut_site") is None:
        return (f"no dominant linearisation site ({cut.get('eligible_reads', 0)} unfragmented plasmid "
                f"reads; best start support {cut.get('start_support_fraction')})")
    txt = (f"linearisation site {cut['cut_site']} (molecule end {cut['cut_end']}, end offset "
           f"{cut['end_offset']:+d}; source {cut['source']}")
    if cut.get("eligible_reads"):
        txt += (f"; {cut['start_support']}/{cut['eligible_reads']} unfragmented reads start and "
                f"{cut['end_support']} end there")
    if cut.get("enzyme"):
        txt += f"; enzyme {cut['enzyme']}"
    if cut.get("enzyme_guess"):
        txt += f"; consistent with {cut['enzyme_guess']}"
    if cut.get("overhang"):
        txt += f"; overhang {cut['overhang']}"
    return txt + ")"


def flag_for_tab(tab_file, insertions, bamfile, plasmid_fasta=None, enzyme=None, cut_site=None,
                 overhang=None, adapter_fasta=None):
    """Write <tab stem>_artefacts.tsv and <tab stem>_cut_site.tsv next to tab_file.

    ``insertions``: Insertion records of e_insert_finder_from_bam (or tuples
    starting with read_id, query_start, query_end, ..., insert_type at index 6).
    Returns the rows."""
    calls = []
    for ins in insertions:
        if hasattr(ins, "insert_type"):
            calls.append((ins.read_id, ins.query_start, ins.query_end, ins.insert_type))
        else:
            calls.append((ins[0], ins[1], ins[2], ins[6] if len(ins) > 6 else ins[3]))
    rows, cut = flag_calls(calls, bamfile, plasmid_fasta, enzyme, cut_site, overhang, adapter_fasta)
    write_rows(rows, artefacts_path_for(tab_file))
    write_cut_site(cut, cut_site_path_for(tab_file))
    logger.info("%s: %s", tab_file, describe_cut(cut))
    counts = Counter(f for r in rows for f in split_flags(r["flags"])[0])
    n_art = sum(bool(split_flags(r["flags"])[0]) for r in rows)
    if n_art:
        logger.warning("%s: %d of %d call(s) carry artefact flags (%s); see %s", tab_file, n_art,
                       len(rows), ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())),
                       artefacts_path_for(tab_file))
    else:
        logger.info("%s: no artefact flags on %d call(s)", tab_file, len(rows))
    return rows


# ----------------------------------------------------- sites (support/independence)

SITE_TOLERANCE = 20
SAME_MOLECULE_TOL = 20
SITE_COLUMNS = ["site_id", "label", "position", "n_calls", "n_reads", "n_plus", "n_minus",
                "n_distinct_ends", "min_molecules", "n_artefact_flagged", "single_molecule_support",
                "query_ids"]


def _to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def cluster_positions(items, L=None, tol=SITE_TOLERANCE):
    """items: [(pos, payload)] -> clusters (lists of items), chaining positions
    within tol (first/last cluster merged across the origin when circular)."""
    pts = sorted(items, key=lambda x: x[0])
    cs, cur = [], []
    for it in pts:
        if cur and it[0] - cur[-1][0] > tol:
            cs.append(cur)
            cur = []
        cur.append(it)
    if cur:
        cs.append(cur)
    if L and len(cs) > 1 and cs[0][0][0] + L - cs[-1][-1][0] <= tol:
        cs[0] = cs.pop() + cs[0]
    return cs


def could_be_same_molecule(a, b, L=None, tol=SAME_MOLECULE_TOL):
    """Two supporting reads that may come from one molecule: opposite strands
    with the same molecule ends (duplex template/complement, or a fold-back
    split into two reads), or one strand with one read starting where the
    other ends (basecaller split). Several calls in one read are merged before
    (site_table counts reads)."""
    if a["read_id"] == b["read_id"]:
        return True
    sa, ea = _to_int(a.get("start")), _to_int(a.get("end"))
    sb, eb = _to_int(b.get("start")), _to_int(b.get("end"))
    if None in (sa, ea, sb, eb):
        return a.get("strand") != b.get("strand")
    near = (lambda x, y: circ_dist(x, y, L) <= tol) if L else (lambda x, y: abs(x - y) <= tol)
    if a.get("strand") != b.get("strand"):
        return near(sa, sb) and near(ea, eb)
    # basecaller split: one read continues where the other stops (and they do
    # not both span the whole linear molecule)
    if near(sa, sb) or near(ea, eb):
        return False
    return near(ea, sb) or near(eb, sa)


def site_table(calls, L=None, tol=SITE_TOLERANCE):
    """calls: dicts with query_id, read_id, label, pos, strand, start, end,
    artefact (bool). Returns (site rows, {query_id: site_id})."""
    by_label = defaultdict(list)
    for c in calls:
        if c.get("pos") is None:
            continue
        by_label[c["label"]].append((c["pos"], c))
    rows, member = [], {}
    for label in sorted(by_label, key=str):
        for cl in cluster_positions(by_label[label], L, tol):
            cs = [c for _, c in cl]
            good = [c for c in cs if not c.get("artefact")]
            reads = {}
            for c in good:
                reads.setdefault(c["read_id"], c)
            rl = list(reads.values())
            n_plus = sum(c.get("strand") == "+" for c in rl)
            n_minus = sum(c.get("strand") == "-" for c in rl)
            # one molecule gives at most one read per strand (template/complement,
            # duplex) - or two pieces of one read split by the basecaller
            min_molecules = max(n_plus, n_minus, (len(rl) + 1) // 2)
            if len(rl) == 2 and could_be_same_molecule(rl[0], rl[1], L):
                min_molecules = 1
            elif len(rl) == 2:
                min_molecules = 2
            sid = f"site_{len(rows) + 1:04d}"
            ends = {(_to_int(c.get("start")), _to_int(c.get("end"))) for c in rl}
            pos = sorted(p for p, _ in cl)[len(cl) // 2]
            rows.append({"site_id": sid, "label": label, "position": pos, "n_calls": len(cs),
                         "n_reads": len(rl), "n_plus": sum(c.get("strand") == "+" for c in rl),
                         "n_minus": sum(c.get("strand") == "-" for c in rl),
                         "n_distinct_ends": len(_distinct(ends, L)),
                         "min_molecules": min_molecules,
                         "n_artefact_flagged": len(cs) - len(good),
                         "single_molecule_support": bool(rl) and min_molecules <= 1,
                         "query_ids": ",".join(c["query_id"] for c in cs)})
            for c in cs:
                member[c["query_id"]] = sid
    return rows, member


def _distinct(ends, L, tol=SAME_MOLECULE_TOL):
    out = []
    for s, e in sorted(ends, key=lambda x: (x[0] is None, x)):
        if s is None or e is None:
            out.append((s, e))
            continue
        if not any(o[0] is not None and circ_dist(o[0], s, L) <= tol and circ_dist(o[1], e, L) <= tol
                   for o in out):
            out.append((s, e))
    return out


def cross_barcode_sites(per_sample, L=None, tol=SITE_TOLERANCE, bleed_pct=5.0):
    """per_sample: {sample: [site rows]} -> rows of sites_across_barcodes.tsv.
    A site (label + position within tol) of sample X is a possible
    bleed-through when another sample has > 100/bleed_pct times its reads."""
    items = []
    for sample, rows in per_sample.items():
        for r in rows:
            if r["n_reads"]:
                items.append((r["position"], dict(r, sample=sample)))
    by_label = defaultdict(list)
    for pos, r in items:
        by_label[r["label"]].append((pos, r))
    out = []
    samples = sorted(per_sample)
    for label in sorted(by_label, key=str):
        for cl in cluster_positions(by_label[label], L, tol):
            counts = Counter()
            for _, r in cl:
                counts[r["sample"]] += r["n_reads"]
            top = max(counts.values())
            pos = sorted(p for p, _ in cl)[len(cl) // 2]
            bleed = sorted(s for s, n in counts.items() if n * 100.0 < bleed_pct * top)
            row = {"label": label, "position": pos, "n_samples": len(counts),
                   "total_reads": sum(counts.values()),
                   "possible_bleed_through": ",".join(bleed) or "-",
                   "site_ids": ",".join(f"{r['sample']}:{r['site_id']}" for _, r in cl)}
            for s in samples:
                row[f"reads_{s}"] = counts.get(s, 0)
            out.append(row)
    return out


# ----------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("tab", help="insertion_from_bam.tab of the same BAM")
    ap.add_argument("--plasmid", default=None, help="plasmid FASTA (default: found next to the BAM)")
    ap.add_argument("--enzyme", default=None, help=f"restriction enzyme ({', '.join(sorted(ENZYMES))})")
    ap.add_argument("--cut_site", type=int, default=None, help="0-based first base of the linear molecule")
    ap.add_argument("--overhang", default=None, help="5' overhang sequence (prefix '-' for a 3' overhang)")
    ap.add_argument("--adapter_fasta", default=None, help="adapter/barcode motifs (replaces the built-in list)")
    args = ap.parse_args(argv)
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    calls = []
    with open(args.tab) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7 and f[1].isdigit():
                calls.append((f[0], int(f[1]), int(f[2]), f[6]))
    rows, cut = flag_calls(calls, args.bam, args.plasmid, args.enzyme, args.cut_site, args.overhang,
                           args.adapter_fasta)
    write_rows(rows, artefacts_path_for(args.tab))
    write_cut_site(cut, cut_site_path_for(args.tab))
    print(describe_cut(cut), file=sys.stderr)


if __name__ == "__main__":
    main()
