"""Direct-repeat (target-site duplication, TSD) finder.

For every insertion in ``<sample>/best_alignment.tab`` the read window
``read[start - open_gap : end + open_gap]`` is BLASTed against the trap
plasmid. The plasmid HSP whose read end is closest to the insertion start is
the *left flank*, the HSP whose read start is closest to the insertion end is
the *right flank*. Both flanks extend into the duplicated target site, so the
overlap of the two flanks on the plasmid is the direct repeat.

* Coordinates: insertion coordinates in the ``Query ID`` (``{read}_{start}_{end}``)
  are 0-based, half-open read coordinates. All reported plasmid coordinates are
  1-based and lie in ``1..len(plasmid)``.
* Circular plasmid: flanks are searched against the plasmid concatenated with
  itself, so flanks spanning the origin are found as one HSP; overlaps are
  computed modulo the plasmid length and DR sequences are extracted circularly.
* Strand: a read from the minus strand has its left flank on the plasmid-right
  side of the insertion; this is handled explicitly.
* If the annotation hit is an IS from the curated database (``Note == IS_DB``),
  the IS boundaries on the read are used to clip flank HSPs that extended into
  the IS by chance base matches (the IS database entry cannot extend past its
  own ends, so its boundaries are the more reliable ones), and to extend flanks
  that stop up to MAX_EXTEND bp short of the IS because of a sequencing error.
  IS hits that miss terminal IS bases are extended by the missing bases; pass
  ``is_fasta`` (IS lengths) to do this at both IS ends. Without IS boundaries
  (Genome / novel inserts) a flank can extend into the insert by chance base
  matches, so the DR length is then an upper bound (+1 bp in ~1/3 of cases).

Output: ``<input_dir>/insertions/<sample>_insertions.csv``, one row per
insertion. Columns (unchanged ones first):

``Read ID, Subject ID, Query Start CSV, Query End CSV`` - read and insertion
coordinates as in the input.
``Start Subject, End Subject`` - plasmid coordinates of the target site: first
base after the insertion on the plasmid-right side and last base before it on
the plasmid-left side. For a DR these are the first and last base of the
duplicated sequence.
``Insertion Point`` - 1-based plasmid coordinate of the first base of the
target site (= simulator ``ref_pos + 1``).
``Repeat Length`` - signed overlap of the two flanks: > 0 direct repeat,
0 no duplication (blunt insertion), < 0 target-site deletion of that many bases;
empty when it cannot be computed.
``Overlap Sequence`` - DR sequence on the plasmid + strand, ``N/A`` otherwise.
``Read Strand`` - ``+``/``-`` orientation of the read relative to the plasmid.
``Status`` - DR, no_duplication, target_site_deletion, repeat_too_long,
left_flank_only, right_flank_only, no_flank_hit, strand_conflict,
read_not_found, bad_query_id.
"""

import logging
import os
import shutil
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__name__)

OUTPUT_COLUMNS = [
    "Read ID", "Subject ID", "Query Start CSV", "Query End CSV",
    "Start Subject", "End Subject", "Insertion Point", "Repeat Length",
    "Overlap Sequence", "Read Strand", "Status",
]
BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq"]
SKIP_FOLDERS = {"insertions", "temp-dr-finder"}
MAX_EXTEND = 10   # bp a flank may be extended to reach a known IS boundary


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def parse_query_id(query_id):
    """``{read_id}_{start}_{end}`` -> (read_id, start, end); read IDs may contain '_'."""
    parts = str(query_id).rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        return parts[0], int(parts[1]), int(parts[2])
    except ValueError:
        return None


def read_plasmid(plasmid_fasta):
    for record in SeqIO.parse(plasmid_fasta, "fasta"):
        seq = str(record.seq).upper()
        if seq:
            return record.id, seq
    raise ValueError(f"No sequence found in plasmid FASTA {plasmid_fasta}")


def circular_slice(seq, start0, length):
    """``length`` bases of circular ``seq`` starting at 0-based ``start0``."""
    n = len(seq)
    if length <= 0 or n == 0:
        return ""
    start0 %= n
    reps = (start0 + length) // n + 1
    return (seq * reps)[start0:start0 + length]


def norm_pos(pos1, plen):
    """Map a 1-based coordinate on the doubled plasmid into 1..plen."""
    return (int(pos1) - 1) % plen + 1


def _to_int(value):
    try:
        if pd.isna(value):
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _run(cmd, what, cwd=None):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    except FileNotFoundError:
        raise RuntimeError(f"'{cmd[0]}' not found in PATH; it is required for {what}.") from None
    if res.returncode != 0:
        raise RuntimeError(f"{what} failed (exit {res.returncode}): {' '.join(cmd)}\n{res.stderr.strip()}")
    return res


# makeblastdb/blastn cannot handle spaces in database paths, so both run inside
# a temporary work_dir (no spaces) with relative file names.
PLASMID_DB = "plasmid_x2_db"


def build_plasmid_db(plasmid_seq, work_dir, name="plasmid"):
    """Write the plasmid concatenated with itself to work_dir and make a BLAST DB of it."""
    os.makedirs(work_dir, exist_ok=True)
    with open(os.path.join(work_dir, "plasmid_x2.fasta"), "w") as fh:
        fh.write(f">{name}_x2\n{plasmid_seq}{plasmid_seq}\n")
    _run(["makeblastdb", "-in", "plasmid_x2.fasta", "-dbtype", "nucl", "-out", PLASMID_DB],
         "makeblastdb", cwd=work_dir)
    return os.path.join(work_dir, PLASMID_DB)


def run_blast(query_name, out_name, work_dir, threads=1):
    """blastn work_dir/query_name against the doubled plasmid DB in work_dir."""
    _run(["blastn", "-query", query_name, "-db", PLASMID_DB, "-out", out_name,
          "-outfmt", "6 " + " ".join(BLAST_FIELDS),
          "-evalue", "1e-5", "-dust", "no", "-soft_masking", "false",
          "-word_size", "16", "-num_threads", str(max(1, int(threads)))],
         "blastn (plasmid flanks)", cwd=work_dir)


def load_hits(blast_out):
    hits = {}
    if not os.path.exists(blast_out) or os.path.getsize(blast_out) == 0:
        return hits
    df = pd.read_csv(blast_out, sep="\t", header=None, names=BLAST_FIELDS, dtype={"qseqid": str})
    for qid, grp in df.groupby("qseqid"):
        hits[qid] = grp.to_dict("records")
    return hits


# --------------------------------------------------------------------------
# core logic (pure, unit-testable)
# --------------------------------------------------------------------------
def flank_junction(h, side, limit=None):
    """Plasmid coordinate (doubled space) of the flank base adjacent to the insertion.

    side "L": last plasmid base of the HSP whose read position is <= limit.
    side "R": first plasmid base of the HSP whose read position is >= limit.
    limit=None means no clipping. The gapped alignment (qseq/sseq) is walked so
    clipping is exact in plasmid coordinates; without it an ungapped offset is used.
    Returns (coordinate, plus_strand).
    """
    plus = h["sstart"] <= h["send"]
    sign = 1 if plus else -1
    rq_s, rq_e = int(h["rq_s"]), int(h["rq_e"])
    if limit is None:
        return int(h["send"] if side == "L" else h["sstart"]), plus
    # HSP stops short of the IS boundary (sequencing error at the junction):
    # extend ungapped up to MAX_EXTEND bases.
    if side == "L" and rq_e <= limit:
        gap = limit - rq_e
        return int(h["send"]) + sign * (gap if gap <= MAX_EXTEND else 0), plus
    if side == "R" and rq_s >= limit:
        gap = rq_s - limit
        return int(h["sstart"]) - sign * (gap if gap <= MAX_EXTEND else 0), plus
    qseq, sseq = h.get("qseq"), h.get("sseq")
    if not isinstance(qseq, str) or not isinstance(sseq, str) or len(qseq) != len(sseq):
        if side == "L":
            k = min(rq_e - limit, rq_e - rq_s)
            return int(h["send"]) - sign * k, plus
        k = min(limit - rq_s, rq_e - rq_s)
        return int(h["sstart"]) + sign * k, plus
    if side == "L":
        r, p = rq_s - 1, int(h["sstart"]) - sign       # last consumed read / plasmid base
        for qc, sc in zip(qseq, sseq):
            if r >= limit:      # limit base consumed: stop, do not take trailing gaps
                break
            if qc != "-":
                r += 1
            if sc != "-":
                p += sign
        return p, plus
    r, p = rq_e + 1, int(h["send"]) + sign
    for qc, sc in zip(reversed(qseq), reversed(sseq)):
        if r <= limit:
            break
        if qc != "-":
            r -= 1
        if sc != "-":
            p -= sign
    return p, plus


def _pick(cands):
    """Highest bit score, then closest to the junction."""
    if not cands:
        return None
    return sorted(cands, key=lambda c: (-c[0]["bitscore"], c[1]))[0][0]


def analyze_insertion(hits, ins_start, ins_end, offset, plasmid_seq, threshold, repeat_th,
                      is_read_start=None, is_read_end=None):
    """Find the DR of one insertion.

    hits: BLAST rows (dicts with BLAST_FIELDS) of the read window vs the doubled
      plasmid; query coordinates are 1-based on the window.
    ins_start, ins_end: 0-based half-open insertion coordinates on the read.
    offset: 0-based read coordinate of the first window base.
    is_read_start, is_read_end: optional 1-based inclusive read coordinates of
      the IS element; flank HSPs are clipped so they do not run into it.
    Returns a dict with the output columns that depend on the alignment.
    """
    plen = len(plasmid_seq)
    res = {"Start Subject": None, "End Subject": None, "Insertion Point": None,
           "Repeat Length": None, "Overlap Sequence": "N/A", "Read Strand": None,
           "Status": "no_flank_hit"}
    left_junction_read = ins_start          # 1-based read position of last base before insertion
    right_junction_read = ins_end + 1       # 1-based read position of first base after insertion

    left_c, right_c = [], []
    for h in hits or []:
        rq_s, rq_e = offset + int(h["qstart"]), offset + int(h["qend"])
        if rq_s <= left_junction_read and rq_e >= right_junction_read:
            continue    # spans the whole insertion: not a flank
        h = dict(h, rq_s=rq_s, rq_e=rq_e)
        d_left = abs(rq_e - left_junction_read)
        if d_left <= threshold and rq_s <= left_junction_read:
            left_c.append((h, d_left))
        d_right = abs(rq_s - right_junction_read)
        if d_right <= threshold and rq_e >= right_junction_read:
            right_c.append((h, d_right))
    left, right = _pick(left_c), _pick(right_c)

    def junction(h, side):
        if side == "L":
            return flank_junction(h, "L", None if is_read_start is None else is_read_start - 1)
        return flank_junction(h, "R", None if is_read_end is None else is_read_end + 1)

    if left is None and right is None:
        return res

    if left is None or right is None:
        h, side = (left, "L") if left is not None else (right, "R")
        pos, plus = junction(h, side)
        res["Read Strand"] = "+" if plus else "-"
        res["Status"] = "left_flank_only" if side == "L" else "right_flank_only"
        # the flank gives either the plasmid-left (a) or plasmid-right (b) end
        gives_a = (side == "L") == plus
        if gives_a:
            res["End Subject"] = norm_pos(pos, plen)
            res["Insertion Point"] = norm_pos(pos + 1, plen)   # exact only for blunt insertions
        else:
            res["Start Subject"] = norm_pos(pos, plen)
            res["Insertion Point"] = norm_pos(pos, plen)
        return res

    lpos, lplus = junction(left, "L")
    rpos, rplus = junction(right, "R")
    if lplus != rplus:
        res["Status"] = "strand_conflict"
        return res
    res["Read Strand"] = "+" if lplus else "-"
    a, b = (lpos, rpos) if lplus else (rpos, lpos)   # a: plasmid-left end, b: plasmid-right start
    d = (a - b) % plen
    if d > plen // 2:
        d -= plen
    dr_len = d + 1
    b_n = norm_pos(b, plen)
    res.update({"Start Subject": b_n, "End Subject": norm_pos(a, plen),
                "Insertion Point": b_n, "Repeat Length": dr_len})
    if dr_len > repeat_th:
        res["Status"] = "repeat_too_long"
    elif dr_len > 0:
        res["Status"] = "DR"
        res["Overlap Sequence"] = circular_slice(plasmid_seq, b_n - 1, dr_len)
    elif dr_len == 0:
        res["Status"] = "no_duplication"
    else:
        res["Status"] = "target_site_deletion"
    return res


def read_is_lengths(is_fasta):
    """IS name -> length; keyed by record ID and by the last header word (older annotations)."""
    lengths = {}
    if not is_fasta:
        return lengths
    for rec in SeqIO.parse(is_fasta, "fasta"):
        lengths.setdefault(rec.id, len(rec.seq))
        lengths.setdefault(rec.description.split()[-1], len(rec.seq))
    return lengths


def _is_bounds(row, q_st, is_lengths=None):
    """1-based inclusive read coordinates of an IS_DB hit, else (None, None).

    If the hit does not reach the IS ends (bases assigned to the flank alignment
    by the mapper, or sequencing errors at the IS end) the bounds are extended by
    the missing IS bases: towards IS position 1 always, towards the IS end when the
    IS length is known (``is_lengths``).
    """
    if str(row.get("Note", "")).strip() != "IS_DB":
        return None, None
    qs, qe = _to_int(row.get("Query Start")), _to_int(row.get("Query End"))
    ss, se = _to_int(row.get("Subject Start")), _to_int(row.get("Subject End"))
    if qs is None or qe is None:
        return None, None
    qs, qe = min(qs, qe), max(qs, qe)
    start, end = q_st + qs, q_st + qe
    if ss is not None and se is not None:
        is_len = (is_lengths or {}).get(str(row.get("Subject ID")))
        missing_first = min(ss, se) - 1
        missing_last = max(0, is_len - max(ss, se)) if is_len else 0
        if ss <= se:            # IS on the read's + strand: IS position 1 is on the left
            start, end = start - missing_first, end + missing_last
        else:
            start, end = start - missing_last, end + missing_first
    return start, end


def find_drs_for_sample(best_alignment_path, reads_fasta, plasmid_seq,
                        open_gap=500, threshold=200, repeat_th=30, threads=1, is_lengths=None):
    """Return a DataFrame (OUTPUT_COLUMNS) with one row per insertion."""
    df = pd.read_csv(best_alignment_path, sep="\t", dtype={"Query ID": str, "Subject ID": str})
    if df.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    n_before = len(df)
    df = df.drop_duplicates(subset="Query ID", keep="first").reset_index(drop=True)
    if len(df) != n_before:
        logger.warning("%s: %d duplicate Query ID rows ignored (first row kept)",
                       best_alignment_path, n_before - len(df))

    tmp = tempfile.TemporaryDirectory(prefix="trapo_dr_")
    work_dir = tmp.name
    if " " in work_dir:
        raise RuntimeError(f"BLAST cannot use a temp directory with spaces ({work_dir}); set TMPDIR")
    reads = SeqIO.index(reads_fasta, "fasta")
    query_fa = os.path.join(work_dir, "flank_windows.fasta")
    jobs, rows = [], []
    try:
        with open(query_fa, "w") as fh:
            for idx, row in df.iterrows():
                base = {"Read ID": None, "Subject ID": row.get("Subject ID"),
                        "Query Start CSV": None, "Query End CSV": None,
                        "Overlap Sequence": "N/A"}
                parsed = parse_query_id(row["Query ID"])
                if parsed is None:
                    rows.append(dict(base, **{"Read ID": row["Query ID"], "Status": "bad_query_id"}))
                    continue
                read_id, q_st, q_end = parsed
                base.update({"Read ID": read_id, "Query Start CSV": q_st, "Query End CSV": q_end})
                if read_id not in reads:
                    rows.append(dict(base, Status="read_not_found"))
                    continue
                seq = str(reads[read_id].seq).upper()
                w0, w1 = max(0, q_st - open_gap), min(len(seq), q_end + open_gap)
                qid = f"ins{idx}"
                fh.write(f">{qid}\n{seq[w0:w1]}\n")
                is_s, is_e = _is_bounds(row, q_st, is_lengths)
                jobs.append((len(rows), qid, q_st, q_end, w0, is_s, is_e))
                rows.append(base)
    finally:
        reads.close()

    hits = {}
    with tmp:
        if jobs:
            build_plasmid_db(plasmid_seq, work_dir)
            run_blast("flank_windows.fasta", "flank_windows_vs_plasmid.tsv", work_dir, threads)
            hits = load_hits(os.path.join(work_dir, "flank_windows_vs_plasmid.tsv"))
    for n, (ri, qid, q_st, q_end, w0, is_s, is_e) in enumerate(jobs, 1):
        rows[ri].update(analyze_insertion(hits.get(qid, []), q_st, q_end, w0, plasmid_seq,
                                          threshold, repeat_th, is_s, is_e))
        if n % 1000 == 0:
            logger.info("  %d/%d insertions analysed", n, len(jobs))

    out = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    for col in ("Query Start CSV", "Query End CSV", "Start Subject", "End Subject",
                "Insertion Point", "Repeat Length"):
        out[col] = out[col].astype("Int64")
    return out


def _find_best_alignment(folder_path):
    files = sorted(f for f in os.listdir(folder_path)
                   if f.endswith("best_alignment.tab") or f.endswith("best_alignments.tab"))
    exact = [f for f in files if f in ("best_alignment.tab", "best_alignments.tab")]
    chosen = (exact or files or [None])[0]
    return os.path.join(folder_path, chosen) if chosen else None


def main_dr_finder(current_dir, plasmid_fasta, open_gap=500, threshold=200, repeat_th=30,
                   force_w=True, threads=1, is_fasta=None):
    """Run the DR finder on every sample folder in ``current_dir``.

    open_gap: read bases taken on each side of the insertion as flank (bp).
    threshold: max distance (bp, inclusive) between a flank HSP end and the
      insertion boundary on the read.
    repeat_th: max DR length (bp, inclusive) reported as a DR.
    threads: blastn -num_threads.
    is_fasta: optional curated IS FASTA (the one given to blast_annot); its
      lengths make IS-boundary clipping exact at both IS ends.
    Returns the list of written CSV paths.
    """
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    for tool in ("blastn", "makeblastdb"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"'{tool}' not found in PATH; it is required by dr_finder.")
    if open_gap < 30:
        logger.warning("open_gap=%d is very short; flanks may not be found", open_gap)

    _, plasmid_seq = read_plasmid(plasmid_fasta)
    is_lengths = read_is_lengths(is_fasta)
    ins_dir = os.path.join(current_dir, "insertions")
    os.makedirs(ins_dir, exist_ok=True)

    written = []
    for folder in sorted(os.listdir(current_dir)):
        folder_path = os.path.join(current_dir, folder)
        if folder in SKIP_FOLDERS or not os.path.isdir(folder_path):
            continue
        out_file = os.path.join(ins_dir, f"{folder}_insertions.csv")
        if os.path.exists(out_file) and not force_w:
            logger.info("%s exists, skipping", out_file)
            continue
        best = _find_best_alignment(folder_path)
        if best is None:
            logger.info("%s: no best_alignment.tab, skipping", folder_path)
            continue
        reads_fasta = os.path.join(folder_path, f"{folder}.fasta")
        if not os.path.exists(reads_fasta):
            logger.warning("%s: reads FASTA %s not found, skipping", folder_path, reads_fasta)
            continue
        logger.info("DR finder: %s", folder)
        result = find_drs_for_sample(best, reads_fasta, plasmid_seq, open_gap, threshold,
                                     repeat_th, threads, is_lengths)
        result.to_csv(out_file, index=False)
        counts = result["Status"].value_counts().to_dict() if not result.empty else {}
        logger.info("  %d insertions -> %s %s", len(result), out_file, counts)
        written.append(out_file)
    return written
