#!/usr/bin/env python3
"""Plasmid-genome shared segments (GitHub issue #5).

When the trap plasmid and the host genome share a segment (a resistance gene,
an origin/replicon fragment, an IS copy), three things go wrong:

1. A genomic DNA read that covers the shared segment plus genome sequence on
   either side maps to the plasmid through the shared segment. The genome
   sequence next to it is then reported as an insertion (an SC clip at a
   plasmid position on the edge of the shared segment), and BLAST labels it
   ``Genome``: a false "novel insertion".
2. Every plasmid read covering the shared segment also aligns to the genome
   copy, so the ``map`` step counted it as "Both plasmid and genome" and the
   transposition rate was inflated (on test-data/ a 230 bp shared segment made
   every plasmid read "genome-mapped": rate 1.00).
3. Real inserts that contain the shared sequence align partly to the plasmid
   (not handled here, see README).

This module provides

* ``detect``: plasmid vs genome with minimap2 (``-x asm20``, secondary
  alignments kept so every genome copy is listed; the plasmid is treated as
  circular). Segments of at least ``min_len`` bp and ``min_identity`` are
  written to ``<data_dir>/plasmid_genome_homology.tsv`` by the ``map`` step.
* ``mapstats``: mapping statistics that do not count a read as genome-mapped
  when its genome alignments only cover sequence that is aligned to the plasmid
  (the shared segment), and that drop genomic reads that reach the plasmid
  only through a shared segment (a genome alignment running from
  plasmid-aligned sequence into unaligned sequence). For a read whose genome
  hits lie within plasmid-aligned sequence, the plasmid-unaligned stretches
  (>= 100 bp) are mapped to the genome on their own: an insert inside a shared
  segment is otherwise hidden in a flank-[insertion]-flank genome alignment.
* ``flag_insertions``: per insertion call, whether it is a homology artefact:
  every plasmid flank of the insert ends inside a shared segment AND the read's
  genome alignment runs contiguously across that flank/insert junction, i.e.
  the read is genome DNA that only looks plasmid-flanked. A real insertion into
  (or next to) a shared segment is not flagged: its insert is not next to the
  shared segment in the genome. Written next to insertion_from_bam.tab as
  ``insertion_from_bam_homology_flags.tsv``; g_blast_annot relabels flagged
  calls ``Contamination`` (column ``Homology`` keeps the BLAST note).

Plasmid coordinates are 0-based half-open; ``plasmid_end`` may exceed
``plasmid_len`` for a segment that crosses the origin of the circular plasmid.

Command line (used by scripts/01_data_prep_from_py.sh):

    plasmid_genome_homology.py detect PLASMID.fasta GENOME.fasta OUT.tsv
        [--min-len 100] [--min-identity 0.8] [--warn-len 500] [--threads 1]
    plasmid_genome_homology.py mapstats PLASMID.bam GENOME.bam TOTAL OUT.txt
        [--genome GENOME.fasta] [--threads 1]
"""

import argparse
import csv
import glob
import logging
import os
import subprocess
import sys
import tempfile

logger = logging.getLogger(__name__)

HOMOLOGY_NAME = "plasmid_genome_homology.tsv"
FLAGS_SUFFIX = "_homology_flags.tsv"
HOMOLOGY_COLUMNS = ["plasmid", "plasmid_len", "plasmid_start", "plasmid_end", "genome_contig",
                    "genome_start", "genome_end", "strand", "identity", "length"]
FLAG_COLUMNS = ["Query ID", "read_id", "query_start", "query_end", "insert_type", "junctions",
                "junctions_in_shared_segment", "junctions_genome_contiguous", "shared_segments",
                "homology_flag"]

DEFAULT_MIN_LEN = 100          # shortest shared segment reported
DEFAULT_MIN_IDENTITY = 0.8     # matches / alignment block length (minimap2 PAF columns 10/11)
DEFAULT_WARN_LEN = 500         # = default insertion threshold of insert_finder
MINIMAP2_PRESET = "asm20"
ORIGIN_EXTENSION = 50000       # bp of the plasmid start appended to find segments crossing the origin

# CIGAR insertions longer than this split an alignment (an inserted element is
# not aligned sequence; same value as circular_mapping.LARGE_INSERTION).
LARGE_INSERTION = 50
# A genome alignment must extend this many bp on both sides of a
# plasmid/non-plasmid boundary to count as contiguous across it (capped by the
# length of the flank and of the insert).
JUNCTION_SPAN = 100
JUNCTION_SLACK = 10
# A read is genome-mapped (mapping stats) only with at least this many
# genome-aligned bases that are not aligned to the plasmid.
MIN_GENOME_ONLY_BASES = 50
# A plasmid flank "ends inside" a shared segment if its boundary lies within
# this many bp of the segment.
SEGMENT_TOLERANCE = 50
# Plasmid alignment pieces ending this close to an insert boundary flank it.
FLANK_TOLERANCE = 20

M, I, D, N, S, H, P, EQ, X = range(9)
QUERY_OPS = {M, I, S, H, EQ, X}
REF_OPS = {M, D, N, EQ, X}


# ----------------------------------------------------------------- FASTA

def read_fasta(path):
    """[(name, sequence)] in file order (name = first word of the header)."""
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


# ------------------------------------------------------------- detection

def _circular_ranges(start, end, length):
    """[start, end) with end possibly > length -> list of ranges in [0, length)."""
    if length <= 0 or end - start >= length:
        return [(0, length)]
    start %= length
    end = start + (end - start)
    if end <= length:
        return [(start, end)]
    return [(start, length), (0, end - length)]


def _contained(inner, outer):
    """Every range of ``inner`` lies inside one range of ``outer``."""
    return all(any(os_ <= s and e <= oe for os_, oe in outer) for s, e in inner)


def parse_paf_hits(lines, lengths, min_len, min_identity):
    """PAF lines of plasmid (query, possibly origin-extended) vs genome ->
    segment dicts in plasmid coordinates (deduplicated)."""
    hits = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qname, qstart, qend = f[0], int(f[2]), int(f[3])
        plen = lengths.get(qname)
        if plen is None:
            continue
        matches, block = int(f[9]), int(f[10])
        identity = matches / block if block else 0.0
        if qstart >= plen:                  # entirely in the appended copy of the start
            qstart, qend = qstart - plen, qend - plen
        length = qend - qstart
        if length < min_len or identity < min_identity:
            continue
        hits.append({"plasmid": qname, "plasmid_len": plen, "plasmid_start": qstart,
                     "plasmid_end": qend, "genome_contig": f[5], "genome_start": int(f[7]),
                     "genome_end": int(f[8]), "strand": f[4], "identity": round(identity, 4),
                     "length": length})
    # drop hits contained in a longer one (same genome locus): copies found
    # twice through the origin extension, or pieces of an origin-crossing hit
    hits.sort(key=lambda h: (-h["length"], -h["identity"], h["genome_contig"], h["genome_start"]))
    kept = []
    for h in hits:
        rng = _circular_ranges(h["plasmid_start"], h["plasmid_end"], h["plasmid_len"])
        dup = False
        for k in kept:
            if (k["plasmid"] == h["plasmid"] and k["genome_contig"] == h["genome_contig"]
                    and k["genome_start"] <= h["genome_start"] and h["genome_end"] <= k["genome_end"]
                    and _contained(rng, _circular_ranges(k["plasmid_start"], k["plasmid_end"],
                                                         k["plasmid_len"]))):
                dup = True
                break
        if not dup:
            kept.append(h)
    kept.sort(key=lambda h: (h["plasmid"], h["plasmid_start"], h["genome_contig"], h["genome_start"]))
    return kept


def detect_homology(plasmid_fasta, genome_fasta, min_len=DEFAULT_MIN_LEN,
                    min_identity=DEFAULT_MIN_IDENTITY, threads=1, preset=MINIMAP2_PRESET):
    """Shared segments between every plasmid record and the genome (minimap2)."""
    records = read_fasta(plasmid_fasta)
    if not records:
        raise ValueError(f"no sequences in {plasmid_fasta}")
    lengths = {name: len(seq) for name, seq in records}
    with tempfile.TemporaryDirectory(prefix="trapo_homology_") as tmp:
        query = os.path.join(tmp, "plasmid_extended.fasta")
        with open(query, "w") as fh:
            for name, seq in records:
                # circular plasmid: append its start so a segment crossing the origin is found whole
                fh.write(f">{name}\n{seq}{seq[:min(len(seq), ORIGIN_EXTENSION)]}\n")
        cmd = ["minimap2", "-c", "-x", preset, "--secondary=yes", "-N", "1000", "-p", "0.05",
               "-t", str(max(1, int(threads))), os.path.abspath(genome_fasta), query]
        try:
            res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True)
        except FileNotFoundError:
            raise RuntimeError("minimap2 not found in PATH") from None
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"minimap2 failed (exit {e.returncode}): {e.stderr.strip()[-2000:]}") from None
    return parse_paf_hits(res.stdout.splitlines(), lengths, min_len, min_identity)


def write_homology(segments, path):
    """Write the segment table atomically (header always written)."""
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)), f".tmp_{os.getpid()}_{os.path.basename(path)}")
    try:
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, HOMOLOGY_COLUMNS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for s in segments:
                w.writerow(s)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_homology(path):
    """Segments from a table written by write_homology ([] if no rows)."""
    segs = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for k in ("plasmid_len", "plasmid_start", "plasmid_end", "genome_start", "genome_end", "length"):
                row[k] = int(row[k])
            row["identity"] = float(row["identity"])
            segs.append(row)
    return segs


def homology_messages(segments, warn_len=DEFAULT_WARN_LEN):
    """(warnings, infos): human-readable summary lines for a segment table."""
    long_ = [s for s in segments if s["length"] >= warn_len]
    short = [s for s in segments if s["length"] < warn_len]
    warnings, infos = [], []

    def describe(ss):
        loci = {}
        for s in ss:
            key = (s["plasmid"], s["plasmid_start"], s["plasmid_end"])
            loci.setdefault(key, []).append(s)
        parts = []
        for (p, a, b), cs in list(loci.items())[:5]:
            parts.append(f"{p}:{a}-{b} ({b - a} bp, {len(cs)} genome cop{'y' if len(cs) == 1 else 'ies'}, "
                         f"max identity {max(c['identity'] for c in cs):.3f})")
        more = f" and {len(loci) - 5} more" if len(loci) > 5 else ""
        return "; ".join(parts) + more

    if long_:
        warnings.append(
            f"plasmid and genome share {len(long_)} segment(s) >= {warn_len} bp: {describe(long_)}. "
            f"Genomic reads through such a segment look like plasmid reads with an insertion; "
            f"these calls are flagged (insertion_from_bam_homology_flags.tsv, column Homology of "
            f"best_alignment.tab) and labelled Contamination.")
    if short:
        infos.append(
            f"plasmid and genome share {len(short)} segment(s) shorter than {warn_len} bp: "
            f"{describe(short)}. They cannot create insertion calls; plasmid reads that align to the "
            f"genome only through them are not counted as genome-mapped in 08_mapping_stats.txt.")
    if not segments:
        infos.append("no plasmid-genome shared segments found.")
    return warnings, infos


# ------------------------------------------------------- alignment pieces

def read_length(rec):
    return sum(n for op, n in rec.cigartuples if op in QUERY_OPS)


def aligned_pieces(rec, max_ins=LARGE_INSERTION):
    """Aligned pieces of a record in ORIGINAL read orientation:
    [(q_start, q_end, ref_at_q_start, ref_at_q_end)], split at clips and at
    CIGAR insertions longer than ``max_ins``. ``ref_at_q_*`` are the reference
    coordinates at the piece's read-orientation boundaries."""
    total = read_length(rec)
    out = []
    q, ref = 0, rec.reference_start
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
    if rec.is_reverse:
        return [(total - e, total - s, re, rs) for s, e, rs, re in out]
    return out


def _merge(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def _overlap(a, b, merged):
    """Bases of [a, b) covered by the merged intervals."""
    return sum(max(0, min(b, e) - max(a, s)) for s, e in merged)


def read_pieces(bam_path, names=None, max_ins=LARGE_INSERTION):
    """{read: [(q_start, q_end, ref_at_q_start, ref_at_q_end, contig)]} for the
    primary + supplementary records of a BAM (optionally only ``names``)."""
    import pysam
    out = {}
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or not rec.cigartuples:
                continue
            if names is not None and rec.query_name not in names:
                continue
            pcs = out.setdefault(rec.query_name, [])
            for p in aligned_pieces(rec, max_ins):
                pcs.append(p + (rec.reference_name,))
    return out


# --------------------------------------------------------- mapping stats

def classify_genome_hits(plasmid_ivs, genome_ivs, span=JUNCTION_SPAN, min_extra=MIN_GENOME_ONLY_BASES):
    """One read: 'genomic' (a genome alignment covers >= span plasmid-aligned
    AND >= span other read bases, i.e. it runs from the plasmid-aligned part
    into the rest of the read: genome DNA reaching the plasmid through a shared
    segment), 'homology_only' (fewer than min_extra genome-aligned bases lie
    outside the plasmid-aligned part) or 'genome'."""
    merged = _merge(plasmid_ivs)
    extra = 0
    for s, e in _merge(genome_ivs):
        extra += (e - s) - _overlap(s, e, merged)
    for s, e in genome_ivs:
        inside = _overlap(s, e, merged)
        if inside >= span and (e - s) - inside >= span:
            return "genomic"
    if extra < min_extra:
        return "homology_only"
    return "genome"


REVCOMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def unaligned_stretches(ivs, length, min_len=JUNCTION_SPAN):
    """Read intervals not covered by ``ivs`` and at least ``min_len`` long."""
    out, pos = [], 0
    for s, e in _merge(ivs):
        if s - pos >= min_len:
            out.append((pos, s))
        pos = max(pos, e)
    if length - pos >= min_len:
        out.append((pos, length))
    return out


def remap_stretches(plasmid_bam, plasmid_ivs, genome_fasta, threads=1, min_bases=MIN_GENOME_ONLY_BASES):
    """Reads (of ``plasmid_ivs``) whose plasmid-unaligned stretches align to
    the genome on their own (minimap2 map-ont) over >= min_bases read bases.

    Needed for reads whose genome alignment runs through a shared segment:
    for an insert inside the shared segment minimap2 aligns
    flank-[insertion]-flank to the genome copy of the segment and does not
    report the insert's own genome hit."""
    import pysam
    found = set()
    with tempfile.TemporaryDirectory(prefix="trapo_homology_") as tmp:
        fa = os.path.join(tmp, "stretches.fasta")
        n = 0
        with pysam.AlignmentFile(plasmid_bam, "rb", check_sq=False) as bam, open(fa, "w") as out:
            for rec in bam.fetch(until_eof=True):
                if (rec.query_name not in plasmid_ivs or rec.is_unmapped or rec.is_secondary
                        or rec.is_supplementary or not rec.query_sequence):
                    continue
                seq = rec.query_sequence
                if rec.is_reverse:
                    seq = seq.translate(REVCOMP)[::-1]
                for s, e in unaligned_stretches(plasmid_ivs[rec.query_name], len(seq)):
                    out.write(f">{rec.query_name}\t{s}\n{seq[s:e]}\n")
                    n += 1
        if not n:
            return found
        cmd = ["minimap2", "-c", "-x", "map-ont", "--secondary=no", "-t", str(max(1, int(threads))),
               os.path.abspath(genome_fasta), fa]
        try:
            res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except FileNotFoundError:
            raise RuntimeError("minimap2 not found in PATH") from None
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"minimap2 failed (exit {e.returncode}): {e.stderr.strip()[-2000:]}") from None
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 12 and int(f[3]) - int(f[2]) >= min_bases:
            found.add(f[0])
    return found


def mapping_stats(plasmid_bam, genome_bam, total=None, genome_fasta=None, threads=1):
    """Read counts for 08_mapping_stats.txt (see module docstring).

    With ``genome_fasta``, reads classified 'homology_only' whose
    plasmid-unaligned stretches (>= JUNCTION_SPAN bp) align to the genome on
    their own are counted as genome-mapped (inserts inside a shared segment)."""
    pl = read_pieces(plasmid_bam)
    gn = read_pieces(genome_bam, names=set(pl))
    counts = {"genome": 0, "homology_only": 0, "genomic": 0}
    hom_only = {}
    for name, pcs in gn.items():
        if not pcs:
            continue
        p_ivs = [(p[0], p[1]) for p in pl[name]]
        kind = classify_genome_hits(p_ivs, [(p[0], p[1]) for p in pcs])
        counts[kind] += 1
        if kind == "homology_only":
            hom_only[name] = p_ivs
    if genome_fasta and hom_only:
        rescued = remap_stretches(plasmid_bam, hom_only, genome_fasta, threads)
        counts["homology_only"] -= len(rescued)
        counts["genome"] += len(rescued)
    plasmid_raw = len(pl)
    both_raw = sum(counts.values())
    plasmid = plasmid_raw - counts["genomic"]
    both = counts["genome"]
    return {"total": total, "plasmid_raw": plasmid_raw, "both_raw": both_raw,
            "homology_only": counts["homology_only"], "genomic": counts["genomic"],
            "plasmid": plasmid, "both": both, "only_plasmid": plasmid - both}


def _rate(num, den):
    return f"{num / den:.2f}" if den > 0 else "NA"


def format_mapping_stats(st):
    lines = []
    if st.get("total") is not None:
        lines.append(f"Total reads: {st['total']}")
    lines += [f"Plasmid: {st['plasmid']}",
              f"Only Plasmid: {st['only_plasmid']}",
              f"Both plasmid and genome: {st['both']}",
              f"Transposition Rate: {_rate(st['both'], st['plasmid'])}",
              f"Genome hits only within plasmid-aligned sequence (not counted as genome-mapped): "
              f"{st['homology_only']}",
              f"Genomic reads reaching the plasmid through a shared segment (excluded): {st['genomic']}",
              f"Plasmid (uncorrected): {st['plasmid_raw']}",
              f"Both plasmid and genome (uncorrected): {st['both_raw']}",
              f"Transposition Rate (uncorrected): {_rate(st['both_raw'], st['plasmid_raw'])}"]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------ flagging

def flags_path_for(tab_file):
    """insertion_from_bam.tab -> insertion_from_bam_homology_flags.tsv (same folder)."""
    return os.path.splitext(tab_file)[0] + FLAGS_SUFFIX


def find_homology_table(*dirs):
    """First plasmid_genome_homology.tsv in the given folders, or None."""
    for d in dirs:
        p = os.path.join(d, HOMOLOGY_NAME)
        if os.path.isfile(p):
            return p
    return None


def find_genome_bam(folder):
    hits = sorted(glob.glob(os.path.join(glob.escape(folder), "07_*.bam")))
    return hits[0] if hits else None


def segments_at(contig, ref, segments, tol=SEGMENT_TOLERANCE):
    """Indices of segments on ``contig`` whose plasmid interval contains
    ``ref`` (+- tol, circular)."""
    out = []
    for i, s in enumerate(segments):
        if s["plasmid"] != contig:
            continue
        L = s["plasmid_len"]
        width = s["plasmid_end"] - s["plasmid_start"]
        if L <= 0 or width >= L:
            out.append(i)
            continue
        d_in = (ref - s["plasmid_start"]) % L
        if d_in <= width + tol or (s["plasmid_start"] - ref) % L <= tol:
            out.append(i)
    return out


def insertion_junctions(qs, qe, pieces, tol=FLANK_TOLERANCE):
    """Plasmid flanks of the insert [qs, qe): [(side, read_pos, ref, contig, flank_len)]."""
    out = []
    left = [p for p in pieces if abs(p[1] - qs) <= tol and p[0] < qs]
    right = [p for p in pieces if abs(p[0] - qe) <= tol and p[1] > qe]
    if left:
        p = min(left, key=lambda p: (abs(p[1] - qs), -(p[1] - p[0])))
        out.append(("left", qs, p[3], p[4], p[1] - p[0]))
    if right:
        p = min(right, key=lambda p: (abs(p[0] - qe), -(p[1] - p[0])))
        out.append(("right", qe, p[2], p[4], p[1] - p[0]))
    return out


def genome_contiguous(side, pos, flank_len, insert_len, genome_ivs, span=JUNCTION_SPAN,
                      slack=JUNCTION_SLACK):
    """A genome alignment covers the read from ``span`` bp inside the flank to
    ``span`` bp inside the insert across the junction at ``pos``."""
    a, b = min(span, flank_len), min(span, insert_len)
    lo, hi = (pos - a, pos + b) if side == "left" else (pos - b, pos + a)
    return any(s <= lo + slack and e >= hi - slack for s, e in genome_ivs)


def flag_insertions(insertions, plasmid_bam, genome_bam, segments):
    """insertions: iterable of (read_id, query_start, query_end, insert_type).
    Returns one dict per insertion (FLAG_COLUMNS)."""
    insertions = list(insertions)
    rows = []
    pl = gn = {}
    if segments and insertions:
        names = {r[0] for r in insertions}
        pl = read_pieces(plasmid_bam, names)
        gn = read_pieces(genome_bam, names)
    for read_id, qs, qe, itype in insertions:
        juncs = insertion_junctions(qs, qe, pl.get(read_id, [])) if segments else []
        g_ivs = [(p[0], p[1]) for p in gn.get(read_id, [])]
        in_seg = cont = 0
        seg_ids = set()
        for side, pos, ref, contig, flank_len in juncs:
            hit = segments_at(contig, ref, segments)
            if hit:
                in_seg += 1
                seg_ids.update(hit)
            if genome_contiguous(side, pos, flank_len, qe - qs, g_ivs):
                cont += 1
        flag = bool(juncs) and in_seg == len(juncs) and cont == len(juncs)
        segs = ",".join(f"{segments[i]['plasmid_start']}-{segments[i]['plasmid_end']}"
                        for i in sorted(seg_ids))
        rows.append({"Query ID": f"{read_id}_{qs}_{qe}", "read_id": read_id, "query_start": qs,
                     "query_end": qe, "insert_type": itype, "junctions": len(juncs),
                     "junctions_in_shared_segment": in_seg, "junctions_genome_contiguous": cont,
                     "shared_segments": segs or "-", "homology_flag": flag})
    return rows


def write_flags(rows, path):
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)), f".tmp_{os.getpid()}_{os.path.basename(path)}")
    try:
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, FLAG_COLUMNS, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_flagged_ids(path):
    """Query IDs flagged in a flags table."""
    with open(path, newline="") as fh:
        return {r["Query ID"] for r in csv.DictReader(fh, delimiter="\t")
                if str(r.get("homology_flag", "")).strip() == "True"}


def flag_for_tab(tab_file, insertions, plasmid_bam, homology_tsv=None, genome_bam=None):
    """Write the flags table for an insertion tab next to it.

    ``insertions``: Insertion records (read_id, query_start, query_end, ...,
    insert_type). The homology table is looked up next to the BAM and in its
    parent folder (the data dir), the genome BAM as 07_*.bam next to the
    plasmid BAM. Without either, no flags are written and a stale flags file is
    removed. Returns the flag rows or None."""
    out = flags_path_for(tab_file)
    bam_dir = os.path.dirname(os.path.abspath(plasmid_bam))
    homology_tsv = homology_tsv or find_homology_table(bam_dir, os.path.dirname(bam_dir))
    genome_bam = genome_bam or find_genome_bam(bam_dir)
    if not homology_tsv or not genome_bam:
        if os.path.exists(out):
            os.remove(out)
        logger.info("%s: no %s or genome BAM; insertion calls not checked for plasmid-genome homology "
                    "(run the map step)", tab_file, HOMOLOGY_NAME)
        return None
    segments = read_homology(homology_tsv)
    rows = flag_insertions([(i.read_id, i.query_start, i.query_end, i.insert_type) for i in insertions],
                           plasmid_bam, genome_bam, segments)
    write_flags(rows, out)
    n = sum(r["homology_flag"] for r in rows)
    if n:
        logger.warning("%s: %d of %d insertion call(s) are genomic reads that only look plasmid-flanked "
                       "through a plasmid-genome shared segment (flagged in %s)", tab_file, n, len(rows), out)
    else:
        logger.info("%s: no insertion call flagged for plasmid-genome homology", tab_file)
    return rows


# ----------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("detect", help="find plasmid-genome shared segments")
    d.add_argument("plasmid")
    d.add_argument("genome")
    d.add_argument("out")
    d.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN)
    d.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY)
    d.add_argument("--warn-len", type=int, default=DEFAULT_WARN_LEN)
    d.add_argument("--threads", type=int, default=1)
    m = sub.add_parser("mapstats", help="mapping statistics corrected for shared segments")
    m.add_argument("plasmid_bam")
    m.add_argument("genome_bam")
    m.add_argument("total", type=int)
    m.add_argument("out")
    m.add_argument("--genome", default=None, help="genome FASTA: re-map plasmid-unaligned stretches of "
                   "reads whose genome hits lie within plasmid-aligned sequence")
    m.add_argument("--threads", type=int, default=1)
    args = ap.parse_args(argv)
    if args.cmd == "detect":
        segs = detect_homology(args.plasmid, args.genome, args.min_len, args.min_identity, args.threads)
        write_homology(segs, args.out)
        warnings, infos = homology_messages(segs, args.warn_len)
        for w in warnings:
            print(f"WARNING: {w}", file=sys.stderr)
        for i in infos:
            print(f"  homology: {i}", file=sys.stderr)
        print(f"  homology: {len(segs)} shared segment(s) (>= {args.min_len} bp, identity >= "
              f"{args.min_identity}) written to {args.out}", file=sys.stderr)
    else:
        st = mapping_stats(args.plasmid_bam, args.genome_bam, args.total, args.genome, args.threads)
        text = format_mapping_stats(st)
        tmp = os.path.join(os.path.dirname(os.path.abspath(args.out)),
                           f".tmp_{os.getpid()}_{os.path.basename(args.out)}")
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, args.out)
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
