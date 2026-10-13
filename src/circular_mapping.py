#!/usr/bin/env python3
"""Circular-aware read mapping to the trap plasmid (AUDIT B09).

The trap plasmid is circular, but minimap2 maps against a linear reference.
A read that crosses the reference origin is split into two alignments, and
when an insertion sits a few hundred bp from the origin, the short plasmid
stretch between the origin and the insertion is often not aligned at all: the
uncovered read interval (the "SC" insert) then contains that plasmid sequence
on top of the element, and its flank is the contig edge instead of the
insertion site.

A second effect has the same symptom: minimap2's long-join (map-ont default)
chains a short plasmid stretch between an insert and a read end (an insertion
a few hundred bp from the restriction site) with the alignment on the other
side of the insert; the alignment then stops at the insert and the stretch is
lost. ``--no-long-join`` aligns it as its own supplementary record.

Approach: the reads are mapped twice, to the plasmid as given ("linear"
mapping, minimap2 defaults, as before) and, with ``--no-long-join``, to the
plasmid rotated by half its length ("rotated" mapping, in which the original
origin lies in the middle of the reference so a stretch next to it is aligned
contiguously with its neighbours). Per read, the linear alignments are kept
unless the rotated ones align more of the read (primary + supplementary,
union of aligned query intervals; CIGAR insertions longer than
``LARGE_INSERTION`` do not count as aligned) by more than ``LARGE_GAIN`` bp,
or by more than ``MIN_GAIN`` bp when a linear alignment stops at a contig end
inside the read (it was cut by the linear origin). Small gains elsewhere are
chance extensions of a flank into the insert and keep the linear mapping, so
most reads (and their CIGAR insertions) are unchanged. Chosen rotated
alignments are converted back to the ORIGINAL plasmid coordinates; an
alignment that crosses the original origin is split there into two records
(the larger part keeps the primary flag, the other becomes supplementary;
insertion/deletion operations that touch the cut are turned into clipping or
dropped, so the pieces stay valid alignments).

So the merged BAM is a normal BAM against the original plasmid (same @SQ
names and lengths), every read keeps exactly one primary record with its full
sequence, and reads whose linear mapping was already complete are passed
through unchanged (byte-identical records).

Tags of converted records: unsplit records keep all tags; split pieces get a
recomputed NM, and the alignment-wide tags (AS, ms, de, s1, s2, cs, MD, cg)
are dropped. SA is rebuilt for every converted read.

Command line (used by scripts/01_data_prep_from_py.sh; needs only pysam):

    circular_mapping.py rotate PLASMID.fasta ROTATED.fasta
    circular_mapping.py merge PLASMID.fasta LINEAR.bam ROTATED.bam OUT.bam

LINEAR.bam and ROTATED.bam are minimap2 outputs for the same FASTQ, in input
(read) order, with unmapped records kept (they keep the two streams in step).
OUT.bam is unsorted and contains mapped reads only.
"""

import argparse
import sys

import pysam

# The rotated alignments are used if they align more than LARGE_GAIN extra
# read bases, or more than MIN_GAIN when a linear alignment stops within
# EDGE_TOLERANCE bp of a contig end inside the read. Ties and bp-level jitter
# keep the linear mapping.
MIN_GAIN = 10
LARGE_GAIN = 30
EDGE_TOLERANCE = 20
# CIGAR insertions longer than this are inserted elements, not aligned bases,
# when the coverage of the two mappings is compared.
LARGE_INSERTION = 50

M, I, D, N, S, H, P, EQ, X = range(9)
QUERY_OPS = {M, I, S, H, EQ, X}
REF_OPS = {M, D, N, EQ, X}
ALIGNED_OPS = {M, EQ, X}
FLAG_SUPPLEMENTARY = 0x800
DROP_TAGS_ON_SPLIT = ("AS", "ms", "de", "s1", "s2", "cs", "MD", "cg", "NM", "SA")


# ------------------------------------------------------------------ FASTA

def read_fasta(path):
    """[(name, sequence)] in file order (name = first word of the header)."""
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0] if line[1:].split() else "", []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def rotation_offset(length):
    """The rotated reference starts at original position ``length // 2``."""
    return length // 2


def write_rotated_fasta(plasmid_fasta, out_fasta):
    """Write every record rotated by half its length (same names)."""
    records = read_fasta(plasmid_fasta)
    if not records:
        raise ValueError(f"no sequences in {plasmid_fasta}")
    with open(out_fasta, "w") as out:
        for name, seq in records:
            s = rotation_offset(len(seq))
            rotated = seq[s:] + seq[:s]
            out.write(f">{name}\n")
            for i in range(0, len(rotated), 80):
                out.write(rotated[i:i + 80] + "\n")


# --------------------------------------------------------------- coverage

def read_length(rec):
    return sum(n for op, n in rec.cigartuples if op in QUERY_OPS)


def query_interval(rec):
    """Aligned query interval [start, end) in original read orientation."""
    cig = rec.cigartuples
    total = read_length(rec)
    left = 0
    for op, n in cig:
        if op in (S, H):
            left += n
        else:
            break
    aligned = sum(n for op, n in cig if op in (M, I, EQ, X))
    start = total - left - aligned if rec.is_reverse else left
    return start, start + aligned


def aligned_intervals(rec, max_ins=None):
    """Query intervals [start, end) in original read orientation that the
    record aligns to the reference. The record's span is split at CIGAR
    insertions longer than ``max_ins`` (default LARGE_INSERTION): an inserted
    element inside one alignment is not aligned sequence."""
    if max_ins is None:
        max_ins = LARGE_INSERTION
    total = read_length(rec)
    ivs, q, start = [], 0, None
    for op, n in rec.cigartuples:
        if op in (S, H):
            if start is not None:
                ivs.append((start, q))
                start = None
            q += n
        elif op == I and n > max_ins:
            if start is not None:
                ivs.append((start, q))
                start = None
            q += n
        elif op in QUERY_OPS:
            if start is None:
                start = q
            q += n
    if start is not None:
        ivs.append((start, q))
    if rec.is_reverse:
        ivs = [(total - e, total - s) for s, e in ivs]
    return ivs


def covered_bases(records):
    """Number of read bases aligned by the primary + supplementary records
    (union over records; large CIGAR insertions do not count)."""
    ivs = sorted(iv for r in records
                 if not r.is_unmapped and not r.is_secondary and r.cigartuples
                 for iv in aligned_intervals(r))
    total, cur_s, cur_e = 0, None, None
    for s, e in ivs:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


# ------------------------------------------------------------- conversion

def _split_cigar_at(cigar_core, ref_start, cut):
    """Split core CIGAR ops (no clips) of an alignment starting at
    ``ref_start`` at reference coordinate ``cut``.

    Returns (left_ops, right_ops) as lists of (op, len)."""
    left, right = [], []
    ref = ref_start
    for op, n in cigar_core:
        if op in REF_OPS:
            if ref + n <= cut:
                left.append((op, n))
            elif ref >= cut:
                right.append((op, n))
            else:
                k = cut - ref
                left.append((op, k))
                right.append((op, n - k))
            ref += n
        else:                       # I (or P): belongs to the side we are on
            (left if ref < cut else right).append((op, n))
    return left, right


def _trim_ends(ops):
    """Strip leading/trailing I/D/N/P ops. Returns (ops, lead_query, lead_ref,
    trail_query): query bases moved into the clips and reference bases skipped
    at the start."""
    lead_q = lead_r = trail_q = 0
    ops = list(ops)
    while ops and ops[0][0] not in ALIGNED_OPS:
        op, n = ops.pop(0)
        if op == I:
            lead_q += n
        elif op in (D, N):
            lead_r += n
    while ops and ops[-1][0] not in ALIGNED_OPS:
        op, n = ops.pop()
        if op == I:
            trail_q += n
    return ops, lead_q, lead_r, trail_q


def _nm(ops, read_seq, q0, ref_seq, r0):
    """Edit distance of core ops against the original reference.
    ``read_seq[q0]`` / ``ref_seq[r0]`` are the first bases of the piece."""
    nm, q, r = 0, q0, r0
    for op, n in ops:
        if op in ALIGNED_OPS:
            for k in range(n):
                if read_seq[q + k].upper() != ref_seq[(r + k) % len(ref_seq)].upper():
                    nm += 1
            q += n
            r += n
        elif op == I:
            nm += n
            q += n
        elif op in (D, N):
            if op == D:
                nm += n
            r += n
    return nm


def _new_record(template, header, flag, ref_start, core, q_start, q_end, total, seq_aln,
                qual_aln, soft, ref_seq):
    """Build one record. q_start/q_end: aligned query interval in alignment
    orientation (full read). seq_aln/qual_aln: full read in alignment
    orientation (or None)."""
    rec = pysam.AlignedSegment(header)
    rec.query_name = template.query_name
    rec.flag = flag
    rec.reference_id = template.reference_id
    rec.reference_start = ref_start
    rec.mapping_quality = template.mapping_quality
    clip = S if soft else H
    cigar = ([(clip, q_start)] if q_start else []) + core + \
            ([(clip, total - q_end)] if total - q_end else [])
    rec.cigartuples = cigar
    if seq_aln is not None:
        lo, hi = (0, total) if soft else (q_start, q_end)
        rec.query_sequence = seq_aln[lo:hi]
        if qual_aln is not None:
            rec.query_qualities = qual_aln[lo:hi]
    for tag, value, vtype in template.get_tags(with_value_type=True):
        if tag not in DROP_TAGS_ON_SPLIT:
            rec.set_tag(tag, value, value_type=vtype)
    if seq_aln is not None and ref_seq is not None:
        rec.set_tag("NM", _nm(core, seq_aln, q_start, ref_seq, ref_start), value_type="i")
    return rec


def convert_record(rec, header, ref_len, ref_seq=None):
    """Convert one record of the rotated mapping to original coordinates.

    Returns a list of 0, 1 or 2 records (2 when the alignment crosses the
    original origin)."""
    if rec.is_unmapped:
        return []
    s = rotation_offset(ref_len)
    cut = ref_len - s                 # rotated coordinate of original position 0
    b0, b1 = rec.reference_start, rec.reference_end
    if b1 <= cut or b0 >= cut:
        out = pysam.AlignedSegment.fromstring(rec.to_string(), header)
        out.reference_start = b0 + s if b1 <= cut else b0 - cut
        return [out]

    cig = rec.cigartuples
    total = read_length(rec)
    lead_clip = 0
    lead_hard = 0
    for op, n in cig:
        if op in (S, H):
            lead_clip += n
            if op == H:
                lead_hard += n
        else:
            break
    trail_clip = 0
    for op, n in reversed(cig):
        if op in (S, H):
            trail_clip += n
        else:
            break
    core = cig[_n_clip_ops(cig, True): len(cig) - _n_clip_ops(cig, False)]

    # Full read in alignment orientation (positions without SEQ are 'N').
    seq_aln = qual_aln = None
    if rec.query_sequence is not None:
        seq = rec.query_sequence
        quals = rec.query_qualities
        seq_aln = "N" * lead_hard + seq + "N" * (total - lead_hard - len(seq))
        if quals is not None:
            qual_aln = [0] * lead_hard + list(quals) + [0] * (total - lead_hard - len(seq))
    full_seq = rec.query_sequence is not None and len(rec.query_sequence) == total

    left_ops, right_ops = _split_cigar_at(core, b0, cut)
    pieces = []
    q = lead_clip
    left_q_len = sum(n for op, n in left_ops if op in QUERY_OPS)
    for ops, r_start, q_start in ((left_ops, b0 + s, q), (right_ops, None, q + left_q_len)):
        ops, lq, lr, tq = _trim_ends(ops)
        if not ops:
            continue
        q_len = sum(n for op, n in ops if op in QUERY_OPS)
        start = q_start + lq
        ref = (r_start if r_start is not None else 0) + lr
        pieces.append((ops, ref, start, start + q_len))
    if not pieces:
        return []

    def aligned(p):
        return sum(n for op, n in p[0] if op in ALIGNED_OPS)
    keep = max(range(len(pieces)), key=lambda i: aligned(pieces[i]))
    out = []
    for i, (ops, ref, qs, qe) in enumerate(pieces):
        flag = rec.flag
        if i != keep and not rec.is_secondary:
            flag |= FLAG_SUPPLEMENTARY
        # Soft clips (full SEQ) only where the record carried the full read and
        # this piece keeps the original role (the primary keeps its full SEQ).
        soft = full_seq and (i == keep or rec.is_supplementary or rec.is_secondary)
        out.append(_new_record(rec, header, flag, ref, ops, qs, qe, total,
                               seq_aln, qual_aln, soft, ref_seq))
    return out


def _n_clip_ops(cig, at_start):
    n = 0
    for op, _ in (cig if at_start else reversed(cig)):
        if op in (S, H):
            n += 1
        else:
            break
    return n


def _rebuild_sa(records, header):
    """Set SA on every primary/supplementary record of one read."""
    main = [r for r in records if not r.is_secondary]
    for r in main:
        others = [o for o in main if o is not r]
        if not others:
            if r.has_tag("SA"):
                r.set_tag("SA", None)
            continue
        parts = []
        for o in others:
            nm = o.get_tag("NM") if o.has_tag("NM") else 0
            parts.append(f"{header.get_reference_name(o.reference_id)},{o.reference_start + 1},"
                         f"{'-' if o.is_reverse else '+'},{o.cigarstring},{o.mapping_quality},{nm};")
        r.set_tag("SA", "".join(parts), value_type="Z")


def convert_read(records, header, ref_lens, ref_seqs):
    out = []
    for rec in records:
        name = header.get_reference_name(rec.reference_id) if rec.reference_id >= 0 else None
        if name is None:
            continue
        out.extend(convert_record(rec, header, ref_lens[name], ref_seqs.get(name)))
    _rebuild_sa(out, header)
    return out


# ------------------------------------------------------------------ merge

def _groups(bam):
    """Yield (read_name, [records]) for consecutive records of one read."""
    name, group = None, []
    for rec in bam:
        if rec.query_name != name and group:
            yield name, group
            group = []
        name = rec.query_name
        group.append(rec)
    if group:
        yield name, group


def edge_inside_read(records, tol=EDGE_TOLERANCE):
    """True if a primary/supplementary alignment stops within ``tol`` bp of a
    contig end at a read position that is neither a read end nor the start of
    another alignment of the read, i.e. the linear origin cut the alignment
    and the read continues unaligned (the B09 situation). A read that simply
    crosses the origin (two abutting alignments) does not count."""
    recs = [r for r in records if not r.is_unmapped and not r.is_secondary and r.cigartuples]
    intervals = [query_interval(r) for r in recs]
    for r, (q_start, q_end) in zip(recs, intervals):
        ref_len = r.header.get_reference_length(r.reference_name)
        total = read_length(r)
        # read-orientation boundary of the alignment at ref_start / ref_end
        q_at_ref_start, q_at_ref_end = (q_end, q_start) if r.is_reverse else (q_start, q_end)
        for at_edge, q in ((r.reference_start <= tol, q_at_ref_start),
                           (ref_len - r.reference_end <= tol, q_at_ref_end)):
            if not at_edge or not tol < q < total - tol:
                continue
            abutting = any(o is not r and (abs(oq_start - q) <= tol or abs(oq_end - q) <= tol)
                           for o, (oq_start, oq_end) in zip(recs, intervals))
            if not abutting:
                return True
    return False


def choose(linear, rotated):
    """'linear' or 'rotated' for one read's mapped records.

    The rotated mapping is used when it aligns more than LARGE_GAIN extra read
    bases, or more than MIN_GAIN extra bases while the linear mapping has an
    alignment cut at a contig end inside the read (origin crossing)."""
    lin = [r for r in linear if not r.is_unmapped]
    rot = [r for r in rotated if not r.is_unmapped]
    if not rot:
        return "linear"
    if not lin:
        return "rotated"
    gain = covered_bases(rot) - covered_bases(lin)
    if gain > LARGE_GAIN or (gain > MIN_GAIN and edge_inside_read(lin)):
        return "rotated"
    return "linear"


def merge(plasmid_fasta, linear_bam, rotated_bam, out_bam):
    """Write the per-read best of the two mappings in original coordinates.

    Returns a dict of counters."""
    ref_seqs = dict(read_fasta(plasmid_fasta))
    counts = {"reads": 0, "mapped": 0, "rotated_used": 0, "split_records": 0}
    with pysam.AlignmentFile(linear_bam, "rb", check_sq=False) as lin, \
            pysam.AlignmentFile(rotated_bam, "rb", check_sq=False) as rot:
        lin_lens = dict(zip(lin.references, lin.lengths))
        rot_lens = dict(zip(rot.references, rot.lengths))
        if lin_lens != rot_lens:
            raise ValueError(f"reference names/lengths differ between {linear_bam} and {rotated_bam}")
        for name, length in lin_lens.items():
            if name in ref_seqs and len(ref_seqs[name]) != length:
                raise ValueError(f"{name}: FASTA length {len(ref_seqs[name])} != BAM length {length}")
        hdr = lin.header.to_dict()
        hdr.setdefault("PG", []).append({
            "ID": "trapo-seq-circular", "PN": "circular_mapping.py",
            "DS": "per read: linear or rotated (by L/2) plasmid mapping, converted to original coordinates",
        })
        header = pysam.AlignmentHeader.from_dict(hdr)
        rot_iter = _groups(rot)
        with pysam.AlignmentFile(out_bam, "wb", header=header) as out:
            for name, lin_group in _groups(lin):
                rot_name, rot_group = next(rot_iter, (None, []))
                if rot_name != name:
                    raise ValueError(f"read order differs between the two mappings at {name} / {rot_name}")
                counts["reads"] += 1
                which = choose(lin_group, rot_group)
                if which == "linear":
                    recs = [r for r in lin_group if not r.is_unmapped]
                    for r in recs:
                        out.write(pysam.AlignedSegment.fromstring(r.to_string(), header))
                else:
                    recs = convert_read([r for r in rot_group if not r.is_unmapped],
                                        header, lin_lens, ref_seqs)
                    counts["rotated_used"] += 1
                    counts["split_records"] += len(recs) - len([r for r in rot_group if not r.is_unmapped])
                    for r in recs:
                        out.write(r)
                if recs:
                    counts["mapped"] += 1
            if next(rot_iter, None) is not None:
                raise ValueError("the rotated mapping has more reads than the linear one")
    return counts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rotate", help="write the plasmid rotated by half its length")
    r.add_argument("plasmid")
    r.add_argument("out")
    m = sub.add_parser("merge", help="merge linear and rotated mappings")
    m.add_argument("plasmid")
    m.add_argument("linear_bam")
    m.add_argument("rotated_bam")
    m.add_argument("out_bam")
    args = ap.parse_args(argv)
    if args.cmd == "rotate":
        write_rotated_fasta(args.plasmid, args.out)
    else:
        c = merge(args.plasmid, args.linear_bam, args.rotated_bam, args.out_bam)
        print(f"  circular merge: {c['mapped']}/{c['reads']} reads mapped, "
              f"{c['rotated_used']} taken from the rotated mapping "
              f"({c['split_records']} extra records from splits at the origin)", file=sys.stderr)


if __name__ == "__main__":
    main()
