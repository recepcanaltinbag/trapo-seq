"""Find large insertions in long reads mapped to the plasmid.

Two kinds of insertion are reported:

* ``IN`` - a CIGAR insertion (``I``) longer than the threshold inside one
  alignment record.
* ``SC`` - a stretch of the read longer than the threshold that is covered by
  none of the read's alignments (primary + supplementary). This catches
  inserts that the aligner soft-clipped or that split the read into several
  alignments (minimap2, bwa-mem).

Secondary alignments are ignored everywhere. Supplementary alignments are used
for the coverage (SC) logic and for CIGAR insertions.

Output (TSV, no header), one row per insertion:
    read_id, query_start, query_end, insertion_length, ref_pos, quality,
    insert_type, is_reverse

* ``query_start``/``query_end``: 0-based half-open interval on the read in its
  original orientation (as written by ``samtools fastq``), so
  ``seq[query_start:query_end]`` is the inserted sequence.
* ``ref_pos``: 0-based reference position next to the insertion. For ``IN``
  it is the position of the first reference base after the insertion; for
  ``SC`` it is the reference coordinate of the boundary of the alignment that
  flanks the gap. The plasmid is circular, so a boundary at the contig end
  (coordinate == contig length) is reported as 0.

Circular plasmid (AUDIT B09): the BAM is expected to come from the
circular-aware mapping (src/circular_mapping.py, run by the ``map`` step), in
which plasmid sequence next to the origin is aligned and an alignment that
crosses the origin is split there into two adjacent records. The gap logic
below then sees only the inserted element. With a plain linear mapping the
contig-edge rule in ``flank_of_gap`` (AUDIT B08) still picks the right flank,
but the insert may then carry plasmid sequence from next to the origin.
* ``quality``: integer part of the mean base quality of the primary record.
* ``is_reverse``: strand of the alignment the insertion was found in / flanked by.

Plasmid-genome shared segments (GitHub issue #5): when the ``map`` step wrote
``plasmid_genome_homology.tsv`` (next to the BAM or in its parent folder) and
the genome BAM ``07_*.bam`` is next to the plasmid BAM, every call is checked
and ``<output stem>_homology_flags.tsv`` is written next to the output (see
src/plasmid_genome_homology.py). A call is flagged when each plasmid flank of
the insert ends inside a shared segment and the read's genome alignment runs
contiguously across that flank/insert junction (a genomic read that only looks
plasmid-flanked). insertion_from_bam.tab itself is unchanged; g_blast_annot
uses the flags.

Artefact evidence (docs/ARTEFACTS.md, src/artefacts.py): every call is also
checked for chimera/ligation/fold-back/adapter/plasmid-multimer signatures
(flank adjacency and strands, distance to the restriction cut site inferred
from the read ends, ...). ``<output stem>_artefacts.tsv`` and
``<output stem>_cut_site.tsv`` are written next to the output; g_blast_annot,
is_stat, dr_finder and novel_elements use them.

Recep Can Altınbağ, 22 10 2024, v0.2
"""

import logging
import os
from collections import OrderedDict, namedtuple

import pysam

from . import artefacts
from . import plasmid_genome_homology as homology

logger = logging.getLogger(__name__)

OUTPUT_NAME = "insertion_from_bam.tab"

# CIGAR operation codes (SAM spec)
M, I, D, N, S, H, P, EQ, X = range(9)
QUERY_OPS = {M, I, S, H, EQ, X}       # ops that advance along the read (incl. hard clips)
REF_OPS = {M, D, N, EQ, X}            # ops that advance along the reference
MATCH_OPS = {M, EQ, X}

Insertion = namedtuple(
    "Insertion",
    "read_id query_start query_end length ref_pos quality insert_type is_reverse",
)

# One alignment record, with the query interval in original read orientation.
Segment = namedtuple("Segment", "q_start q_end ref_start ref_end is_reverse ref_len")


def _clip(cigar, at_start):
    """Total soft+hard clip length at the start or end of a CIGAR."""
    ops = cigar if at_start else reversed(cigar)
    total = 0
    for op, length in ops:
        if op in (S, H):
            total += length
        else:
            break
    return total


def read_length(read):
    """Full read length including hard clips."""
    return sum(length for op, length in read.cigartuples if op in QUERY_OPS)


def segment_of(read, ref_len=None):
    """The aligned query interval of a record in original read orientation."""
    cigar = read.cigartuples
    total = read_length(read)
    aligned = sum(length for op, length in cigar if op in (M, I, EQ, X))
    left = _clip(cigar, True)
    if read.is_reverse:
        q_start = total - left - aligned
    else:
        q_start = left
    return Segment(q_start, q_start + aligned, read.reference_start,
                   read.reference_end, read.is_reverse, ref_len)


def find_large_insertions(read, insertion_threshold):
    """CIGAR insertions longer than the threshold in one alignment record.

    Returns a list of (query_start, query_end, length, ref_pos) with the query
    interval in original read orientation (half-open).
    """
    total = read_length(read)
    q = 0
    ref = read.reference_start
    found = []
    for op, length in read.cigartuples:
        if op == I and length > insertion_threshold:
            if read.is_reverse:
                start = total - q - length
            else:
                start = q
            found.append((start, start + length, length, ref))
        if op in QUERY_OPS:
            q += length
        if op in REF_OPS:
            ref += length
    return found


def uncovered_intervals(segments, length, insertion_threshold):
    """Maximal read intervals [start, end) not covered by any segment and
    longer than the threshold."""
    intervals = sorted((s.q_start, s.q_end) for s in segments)
    gaps = []
    pos = 0
    for start, end in intervals:
        if start > pos and start - pos > insertion_threshold:
            gaps.append((pos, start))
        pos = max(pos, end)
    if length - pos > insertion_threshold:
        gaps.append((pos, length))
    return gaps


# A flank boundary this close to a contig end is treated as "at the edge".
CONTIG_EDGE_TOLERANCE = 20


def _at_contig_edge(ref_coord, seg):
    if ref_coord <= CONTIG_EDGE_TOLERANCE:
        return True
    return seg.ref_len is not None and seg.ref_len - ref_coord <= CONTIG_EDGE_TOLERANCE


def _wrap(ref, seg):
    """Circular plasmid: the contig end (coordinate == length) is position 0."""
    if seg.ref_len and ref >= seg.ref_len:
        return ref - seg.ref_len
    return ref


def flank_of_gap(segments, gap_start, gap_end):
    """Pick the alignment flanking an uncovered read interval.

    Returns (ref_pos, is_reverse) of the flanking alignment boundary (a
    boundary at the contig end is returned as 0), or (None, None) if no
    segment exists.

    A segment ending at ``gap_start`` flanks the gap on the left, one starting
    at ``gap_end`` on the right. The left flank is preferred (legacy
    behaviour) unless its boundary lies on a contig edge (within CONTIG_EDGE_TOLERANCE bp of 0
    or of the contig length)
    and the right one does not: on a circular plasmid a flank that stops at
    the contig edge only reflects the linear reference, not the insertion
    site (e.g. an insertion a few hundred bp after the origin).
    """
    def ref_at_read_end(seg):     # reference coordinate at the segment's right end in read orientation
        return seg.ref_start if seg.is_reverse else seg.ref_end

    def ref_at_read_start(seg):
        return seg.ref_end if seg.is_reverse else seg.ref_start

    left = [s for s in segments if s.q_end == gap_start]
    right = [s for s in segments if s.q_start == gap_end]
    candidates = []
    if left:
        seg = max(left, key=lambda s: s.q_end - s.q_start)
        candidates.append((ref_at_read_end(seg), seg))
    if right:
        seg = max(right, key=lambda s: s.q_end - s.q_start)
        candidates.append((ref_at_read_start(seg), seg))
    if not candidates:
        # No adjacent segment (should not happen for a maximal gap); use the
        # nearest one in read coordinates.
        if not segments:
            return None, None
        seg = min(segments, key=lambda s: min(abs(s.q_end - gap_start), abs(s.q_start - gap_end)))
        if abs(seg.q_end - gap_start) <= abs(seg.q_start - gap_end):
            return _wrap(ref_at_read_end(seg), seg), seg.is_reverse
        return _wrap(ref_at_read_start(seg), seg), seg.is_reverse

    for ref, seg in candidates:
        if not _at_contig_edge(ref, seg):
            return _wrap(ref, seg), seg.is_reverse
    ref, seg = candidates[0]
    return _wrap(ref, seg), seg.is_reverse


def _mean_quality(read):
    qualities = read.query_qualities
    if qualities is None or len(qualities) == 0:
        return 0.0
    return sum(qualities) / len(qualities)


def find_insertions(bamfile, insertion_threshold=500):
    """Return the list of Insertion records for a BAM file."""
    by_read = OrderedDict()      # read_id -> dict(segments, inserts, length, quality, has_primary)
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or not read.cigartuples:
                continue
            info = by_read.setdefault(read.query_name, {
                "segments": [], "inserts": [], "length": 0,
                "quality": None, "has_primary": False,
            })
            ref_len = bam.get_reference_length(read.reference_name)
            seg = segment_of(read, ref_len)
            info["segments"].append(seg)
            info["length"] = max(info["length"], read_length(read))
            if not read.is_supplementary:
                info["quality"] = _mean_quality(read)
                info["has_primary"] = True
            elif info["quality"] is None:
                info["quality"] = _mean_quality(read)   # fallback until the primary is seen
            for start, end, length, ref in find_large_insertions(read, insertion_threshold):
                info["inserts"].append((start, end, length, ref, "IN", read.is_reverse))

    results = []
    for read_id, info in by_read.items():
        quality = int(info["quality"] or 0)
        rows = []
        for gap_start, gap_end in uncovered_intervals(info["segments"], info["length"], insertion_threshold):
            ref, is_reverse = flank_of_gap(info["segments"], gap_start, gap_end)
            rows.append((gap_start, gap_end, gap_end - gap_start, ref, "SC", is_reverse))
            logger.debug("%s: SC %d-%d ref %s reverse %s", read_id, gap_start, gap_end, ref, is_reverse)
        rows.extend(info["inserts"])
        for start, end, length, ref, insert_type, is_reverse in rows:
            if length > insertion_threshold:
                results.append(Insertion(read_id, start, end, length, ref, quality,
                                         insert_type, is_reverse))
    logger.info("%s: %d reads, %d insertions > %d bp", bamfile, len(by_read),
                len(results), insertion_threshold)
    return results


def write_insertions(insertions, output_file):
    """Write insertions as TSV (atomically: temp file, then rename)."""
    tmp = os.path.join(os.path.dirname(os.path.abspath(output_file)),
                       f".tmp_{os.getpid()}_{os.path.basename(output_file)}")
    try:
        with open(tmp, "w") as f:
            for ins in insertions:
                f.write(f"{ins.read_id}\t{ins.query_start}\t{ins.query_end}\t{ins.length}\t"
                        f"{ins.ref_pos}\t{ins.quality}\t{ins.insert_type}\t{ins.is_reverse}\n")
        os.replace(tmp, output_file)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main_insert_finder_from_bam(bamfile, output_file, insertion_threshold=500,
                                homology_tsv=None, genome_bam=None, plasmid_fasta=None,
                                enzyme=None, cut_site=None, overhang=None, adapter_fasta=None,
                                artefact_check=True):
    """Find insertions in ``bamfile`` and write them to ``output_file``.

    Also writes the plasmid-genome homology flags next to ``output_file`` when
    the homology table and the genome BAM are found (``homology_tsv`` /
    ``genome_bam``, default: looked up next to ``bamfile``, see module
    docstring), and the artefact evidence (src/artefacts.py; ``plasmid_fasta``
    default: a FASTA next to the BAM or in its parent folder whose record
    matches the BAM reference; ``enzyme`` or ``cut_site``/``overhang``
    override the cut site inferred from the read ends; ``adapter_fasta``
    replaces the built-in adapter motifs).

    Returns the list of Insertion records.
    """
    insertions = find_insertions(bamfile, insertion_threshold)
    write_insertions(insertions, output_file)
    logger.info("Insertions written to %s", output_file)
    homology.flag_for_tab(output_file, insertions, bamfile, homology_tsv, genome_bam)
    if artefact_check:
        artefacts.flag_for_tab(output_file, insertions, bamfile, plasmid_fasta, enzyme, cut_site,
                               overhang, adapter_fasta)
    return insertions
