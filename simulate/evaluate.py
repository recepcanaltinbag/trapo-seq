#!/usr/bin/env python3
"""
Score trapo-seq output against the ground truth written by simulate_traposome.py.

Usage
-----
    python simulate/evaluate.py -s SIM_DIR -r RESULTS_DIR -o OUT_DIR [--tol 20]

SIM_DIR      contains truth.tsv, clones.tsv (and simulation.json, optional).
RESULTS_DIR  is the pipeline data directory: one folder per barcode with
             `insertion_from_bam.tab` (insert_finder_batch) and
             `best_alignment.tab` (blast_annot_batch), plus
             `insertions/<barcode>_insertions.csv` (dr_finder). Any missing
             file turns the corresponding metrics into N/A.
OUT_DIR      receives metrics.json, metrics.tsv (long format: barcode, section,
             metric, value), per_read.tsv (truth joined to the predictions)
             and sites.tsv (site-level matching).

Conventions scored against
--------------------------
* Coordinates. Truth `ref_pos` is 0-based; the element sits after plasmid
  base ref_pos + tsd_len - 1 and is followed by the duplicated TSD
  plasmid[ref_pos : ref_pos + tsd_len]. Because the TSD is a direct repeat, an
  aligner can legitimately place the insertion anywhere in the 0-based
  interval [ref_pos, ref_pos + tsd_len]. When the plasmid/IS/genome FASTA are
  available (paths from simulation.json, or -p/-i/-g) the window is widened
  by chance extensions of the repeat (plasmid base left of the TSD equal to
  the last insert base, etc.), which no method can resolve. The primary
  position error (`position.error`) is the circular distance from the
  predicted `ref_pos` (0-based, as written by insert_finder) to that window
  (0 if inside); `position.tsd_interval` uses [ref_pos, ref_pos + tsd_len] and
  `position.raw` |pred - ref_pos|. Site matching uses the same window +- --tol.
  DR-finder positions are 1-based BLAST coordinates and are converted (-1).
* TSD. `length_accuracy`/`sequence_accuracy` require the simulated TSD
  exactly; `consistent_accuracy` accepts a reported repeat that contains at
  least tsd_len bases and lies inside the extended observable repeat.
* Read identity. Predictions are matched to truth by read ID. Because
  v0.1.0's dr_finder cannot parse IDs containing '_' (AUDIT B34) the benchmark
  runner may rename `barcode01_r000001` to `barcode01-r000001`; both spellings
  are accepted.
* De-duplication. Pipelines may write several rows per read (several
  insertions, IN and SC duplicates, several BLAST HSPs; AUDIT B02/B03/B12).
  Each read is scored once, using a truth-agnostic "primary" call:
  insert finder -> the longest insertion; annotation -> priority
  IS_DB > Genome > Unannotated > Contamination, then highest BLAST score; DR finder -> the
  row of the primary insertion if present, otherwise the first row with a
  repeat. Row multiplicity is reported separately.
* Read-level detection. A truth read is positive if its category is
  is_insertion or novel_insertion. A read is "called" by the insert finder if
  it has any row in insertion_from_bam.tab; it is "called" by the annotation
  if it has a best_alignment.tab row whose Note is IS_DB, Genome or
  Unannotated (plasmid-flanked insert without IS/genome hit); Contamination
  rows are reported, not counted as biological calls.
* Classification. is_insertion is correct only as IS_DB with the true IS.
  novel_insertion is correct as Genome with a hit overlapping the true genomic
  source; with --novel_absent (the genome given to blast_annot lacks the novel
  inserts, e.g. --novel_mode element annotated with the original genome)
  Unannotated is the correct label instead. Contamination for a true
  insertion is always wrong (contamination_note_true_insertions).
  Fragmented reads may not contain the insertion at all, so recall is also
  reported on unfragmented reads that pass --min_read_len.
* Site level. Predicted sites are built without the truth: reads with the same
  annotation label (IS name for IS_DB, otherwise the Note) are chained by
  position (gap <= --site_link). Each site is matched one-to-one (greedy:
  closest first, ties broken by matching label, then support) to a true insertion event (clones.tsv, n_reads > 0) within
  --tol of its ambiguity window; unmatched sites are spurious.

Only the Python standard library is used.
"""

import argparse
import csv
import glob
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict

INSERTION_CATS = ("is_insertion", "novel_insertion")
NEGATIVE_CATS = ("empty_plasmid", "deletion", "contamination")
NOTE_PRIORITY = {"IS_DB": 0, "Genome": 1, "Unannotated": 2, "Contamination": 3}
CALL_NOTES = ("IS_DB", "Genome", "Unannotated")   # biological insertion calls
NOTE_NAMES = ("IS_DB", "Genome", "Unannotated", "Contamination")
NA = None  # serialised as "N/A" in TSV, null in JSON


# ---------------------------------------------------------------- helpers
def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def to_float(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def ratio(a, b):
    return a / b if b else NA


def f1(p, r):
    if p is None or r is None or (p + r) == 0:
        return NA
    return 2 * p * r / (p + r)


def circ_dist(a, b, plen=None):
    d = abs(a - b)
    if plen:
        d %= plen
        d = min(d, plen - d)
    return d


def interval_dist(pred, lo, hi, plen=None):
    """Distance from `pred` to the closed interval [lo, hi] (circular if plen)."""
    if lo <= pred <= hi:
        return 0
    return min(circ_dist(pred, lo, plen), circ_dist(pred, hi, plen))


def dist_summary(values, windows=(0, 5, 10, 20)):
    vals = [v for v in values if v is not None]
    out = {"n": len(vals)}
    if not vals:
        out.update({"median": NA, "mean": NA, "max": NA})
        for w in windows:
            out[f"pct_within_{w}bp"] = NA
        return out
    out["median"] = statistics.median(vals)
    out["mean"] = round(statistics.mean(vals), 3)
    out["max"] = max(vals)
    for w in windows:
        out[f"pct_within_{w}bp"] = round(100.0 * sum(v <= w for v in vals) / len(vals), 2)
    return out


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- truth
def read_fasta(path):
    recs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        recs[name] = "".join(chunks).upper()
    return recs


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def find_ref(path, sim_dir):
    if not path:
        return None
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for cand in (path, os.path.join(sim_dir, path), os.path.join(here, path)):
        if os.path.exists(cand):
            return cand
    return None


def repeat_extension(plasmid, pos, tsd_len, insert):
    """Bases by which the observable direct repeat extends beyond the TSD.

    Construct = P[:pos+t] + E + P[pos:pos+t] + P[pos+t:]. If P[pos-1] == E[-1]
    the repeat (and the equivalent insertion placements) extends one base to the
    left, if E[0] == P[pos+t] one base to the right, and so on. Such extensions
    are indistinguishable from the TSD in any read, so they are scored as an
    ambiguity window, not as errors.
    """
    n = len(plasmid)
    left = 0
    while left < len(insert) and plasmid[(pos - left - 1) % n] == insert[-left - 1]:
        left += 1
    right = 0
    while right < len(insert) and insert[right] == plasmid[(pos + tsd_len + right) % n]:
        right += 1
    ext = "".join(plasmid[(pos - left + k) % n] for k in range(left + tsd_len + right))
    return left, right, ext


def annotate_repeats(sim_dir, meta, truth, clones, refs=None):
    """Attach the ambiguity window [rep_lo, rep_hi] and extended repeat to every insertion."""
    refs = refs or {}
    paths = {k: find_ref(refs.get(k) or meta.get(k), sim_dir) for k in ("plasmid", "is_fasta", "genome")}
    rows = [t for t in truth.values() if t["category"] in INSERTION_CATS] + \
        [c for c in clones if c["category"] in INSERTION_CATS]
    for r in rows:
        r["rep_lo"], r["rep_hi"] = r["ref_pos"], (r["ref_pos"] or 0) + (r["tsd_len"] or 0)
        r["ext_seq"] = r.get("tsd_seq") or ""
    if not paths["plasmid"] or not paths["is_fasta"]:
        return False
    plasmid = next(iter(read_fasta(paths["plasmid"]).values()))
    is_db = read_fasta(paths["is_fasta"])
    genome = read_fasta(paths["genome"]) if paths["genome"] and any(
        r["category"] == "novel_insertion" for r in rows) else {}
    cache = {}
    for r in rows:
        key = (r["barcode"], r["clone_id"])
        if key not in cache:
            if r["category"] == "is_insertion":
                ins = is_db.get(r["element"])
            else:
                g = genome.get(r.get("genome_contig") or "")
                ins = g[r["genome_start"]:r["genome_end"]] if g and r["genome_start"] is not None else None
            if not ins or r["ref_pos"] is None:
                cache[key] = None
            else:
                if r["orientation"] == "-":
                    ins = revcomp(ins)
                cache[key] = repeat_extension(plasmid, r["ref_pos"], r["tsd_len"] or 0, ins)
        if cache[key]:
            left, right, ext = cache[key]
            r["rep_lo"], r["rep_hi"] = r["ref_pos"] - left, r["ref_pos"] + (r["tsd_len"] or 0) + right
            r["ext_seq"] = ext
    return True


def load_truth(sim_dir, refs=None):
    truth = {}
    for row in read_tsv(os.path.join(sim_dir, "truth.tsv")):
        for k in ("ref_pos", "tsd_len", "insert_len", "genome_start", "genome_end", "read_len"):
            row[k] = to_int(row.get(k))
        row["fragmented"] = str(row.get("fragmented")).strip() in ("True", "true", "1")
        truth[row["read_id"]] = row
    clones = []
    cpath = os.path.join(sim_dir, "clones.tsv")
    if os.path.exists(cpath):
        for row in read_tsv(cpath):
            for k in ("ref_pos", "tsd_len", "insert_len", "genome_start", "genome_end", "n_reads"):
                row[k] = to_int(row.get(k))
            clones.append(row)
    meta = {}
    jpath = os.path.join(sim_dir, "simulation.json")
    if os.path.exists(jpath):
        with open(jpath) as fh:
            meta = json.load(fh)
    plen = meta.get("plasmid_len")
    have_seq = annotate_repeats(sim_dir, meta, truth, clones, refs)
    return truth, clones, plen, have_seq


def make_id_resolver(truth):
    alias = {}
    for rid in truth:
        alias[rid] = rid
        alias[rid.replace("_", "-")] = rid
    return lambda x: alias.get(x)


# ---------------------------------------------------------------- predictions
FINDER_COLS = ["read_id", "query_start", "query_end", "length", "ref_pos", "quality", "type", "is_reverse"]
FINDER_ALIASES = {"insertion_length": "length", "insert_type": "type", "insert_len": "length"}


def parse_finder(path):
    """insertion_from_bam.tab -> list of dicts. Headerless (v0.1.0) or with a header."""
    rows = []
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None and len(f) > 1 and to_int(f[1]) is None:
                header = [FINDER_ALIASES.get(h.strip().lower(), h.strip().lower()) for h in f]
                continue
            cols = header or FINDER_COLS
            rec = dict(zip(cols, f))
            rows.append({
                "read_id": rec.get("read_id"),
                "query_start": to_int(rec.get("query_start")),
                "query_end": to_int(rec.get("query_end")),
                "length": to_int(rec.get("length")),
                "ref_pos": to_int(rec.get("ref_pos")),
                "type": rec.get("type"),
            })
    return rows


def split_query_id(qid):
    """'<read>_<qstart>_<qend>' -> (read, qstart, qend); robust to '_' in read IDs."""
    parts = qid.rsplit("_", 2)
    if len(parts) == 3 and to_int(parts[1]) is not None and to_int(parts[2]) is not None:
        return parts[0], int(parts[1]), int(parts[2])
    return qid, None, None


def parse_annot(path):
    if not path or not os.path.exists(path):
        return None
    rows = []
    for rec in read_tsv(path):
        qid = rec.get("Query ID") or ""
        rid, qs, qe = split_query_id(qid)
        if rec.get("Read ID"):
            rid = rec["Read ID"]
        rows.append({
            "read_id": rid, "qs": qs, "qe": qe,
            "subject": (rec.get("Subject ID") or "").strip(),
            "score": to_float(rec.get("Score")),
            "s_start": to_int(rec.get("Subject Start")),
            "s_end": to_int(rec.get("Subject End")),
            "note": (rec.get("Note") or "").strip(),
            "ref_pos": to_int(rec.get("ref_pos")),
        })
    return rows


def parse_dr(path):
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path, newline="") as fh:
        for rec in csv.DictReader(fh):
            rows.append({
                "read_id": rec.get("Read ID"),
                "qs": to_int(rec.get("Query Start CSV")),
                "qe": to_int(rec.get("Query End CSV")),
                "point": to_int(rec.get("Insertion Point")),
                "rep_len": to_int(rec.get("Repeat Length")),
                "rep_seq": (rec.get("Overlap Sequence") or "").strip(),
                "start": to_int(rec.get("Start Subject")),
            })
    return rows


def first_existing(*patterns):
    for p in patterns:
        hits = sorted(glob.glob(p))
        if hits:
            return hits[0]
    return None


def load_predictions(results_dir, barcode):
    bdir = os.path.join(results_dir, barcode)
    finder = first_existing(os.path.join(glob.escape(bdir), "insertion_from_bam.tab"),
                            os.path.join(glob.escape(bdir), "*insertion_from_bam*.tab"))
    annot = first_existing(os.path.join(glob.escape(bdir), "best_alignment.tab"),
                           os.path.join(glob.escape(bdir), "best_alignment*.tab"))
    dr = first_existing(os.path.join(glob.escape(results_dir), "insertions", f"{barcode}_insertions.csv"),
                        os.path.join(glob.escape(bdir), "*insertions*.csv"))
    return {"finder": parse_finder(finder), "annot": parse_annot(annot), "dr": parse_dr(dr),
            "files": {"finder": finder, "annot": annot, "dr": dr}}


# ---------------------------------------------------------------- per-read selection
def group_by_read(rows, resolve):
    out, unmatched = defaultdict(list), 0
    for r in rows or []:
        rid = resolve(r["read_id"]) if r["read_id"] is not None else None
        if rid is None:
            unmatched += 1
            continue
        out[rid].append(r)
    return out, unmatched


def primary_finder(rows):
    return max(rows, key=lambda r: (r["length"] or 0))


def primary_annot(rows):
    return min(rows, key=lambda r: (NOTE_PRIORITY.get(r["note"], 9), -(r["score"] or 0)))


def primary_dr(rows, finder_row=None):
    if finder_row is not None:
        same = [r for r in rows if r["qs"] == finder_row["query_start"] and r["qe"] == finder_row["query_end"]]
        if same:
            rows = same
    with_rep = [r for r in rows if (r["rep_len"] or 0) > 0]
    return (with_rep or rows)[0]


def genome_overlap(row, t, slack=50):
    if row["s_start"] is None or row["s_end"] is None or t["genome_start"] is None:
        return False
    lo, hi = sorted((row["s_start"], row["s_end"]))
    lo -= 1  # BLAST 1-based inclusive -> 0-based half-open
    return lo < t["genome_end"] + slack and hi > t["genome_start"] - slack


# ---------------------------------------------------------------- TSD
def site_label(r):
    if r.get("annot_note") == "IS_DB":
        return r.get("annot_subject")
    return r.get("annot_note") or None


def tsd_consistent(r):
    """Reported DR is compatible with the truth: at least the TSD, and inside the
    observable (extended) repeat. With tsd_len == 0, a 0-length report is correct."""
    if r["dr_len"] is None:
        return False
    if (r["tsd_len"] or 0) == 0 and r["dr_len"] == 0:
        return True
    seq = (r["dr_seq"] or "").upper()
    return r["dr_len"] >= (r["tsd_len"] or 0) and r["dr_len"] > 0 and seq in (r["ext_seq"] or "").upper()


def tsd_metrics(pos_reads, plen):
    rs = [r for r in pos_reads if r["n_dr_rows"]]
    det = [r for r in pos_reads if r["n_finder_rows"]]
    len_rs = [r for r in rs if r["dr_len"] is not None and r["tsd_len"] is not None]
    seq_rs = [r for r in rs if (r["tsd_len"] or 0) > 0]
    found = [r for r in seq_rs if (r["dr_len"] or 0) > 0]
    pt = [interval_dist(r["dr_point"] - 1, r["rep_lo"], r["rep_hi"], plen)
          for r in rs if r["dr_point"] is not None and r["dr_point"] > 0 and r["rep_lo"] is not None]
    return {
        "available": True, "reads_with_dr_rows": len(rs),
        "pct_detected_insertions_with_dr_row": ratio(100.0 * len(rs), len(det)),
        # strict: exactly the simulated TSD
        "length_accuracy": ratio(sum(r["dr_len"] == r["tsd_len"] for r in len_rs), len(len_rs)),
        "sequence_accuracy": ratio(sum((r["dr_seq"] or "").upper() == (r["tsd_seq"] or "").upper() for r in seq_rs),
                                   len(seq_rs)),
        # tolerant: compatible with the observable repeat (TSD plus chance extensions)
        "consistent_accuracy": ratio(sum(tsd_consistent(r) for r in len_rs), len(len_rs)),
        "pct_tsd_reported_nonzero": ratio(100.0 * len(found), len(seq_rs)),
        "length_accuracy_when_reported": ratio(sum(r["dr_len"] == r["tsd_len"] for r in found), len(found)),
        "consistent_accuracy_when_reported": ratio(sum(tsd_consistent(r) for r in found), len(found)),
        "length_error_when_reported": dist_summary([abs(r["dr_len"] - r["tsd_len"]) for r in found], windows=(0, 1, 2)),
        "insertion_point_error": dist_summary(pt),
    }


# ---------------------------------------------------------------- scoring
def novel_label_ok(label, novel_absent=False):
    """Correct site/read label for a novel insertion (see Classification)."""
    return label == ("Unannotated" if novel_absent else "Genome")


def note_pcts(notes, n):
    return {f"pct_note_{k}": ratio(100.0 * notes.get(k, 0), n) for k in NOTE_NAMES}


def score_barcode(truth_rows, preds, clones, resolve, plen, tol, min_read_len, min_support, site_link=15,
                  novel_absent=False):
    res = {"files": preds["files"]}
    per_read = []
    finder_by, finder_unmatched = group_by_read(preds["finder"], resolve)
    annot_by, annot_unmatched = group_by_read(preds["annot"], resolve)
    dr_by, dr_unmatched = group_by_read(preds["dr"], resolve)
    have_finder, have_annot, have_dr = (preds[k] is not None for k in ("finder", "annot", "dr"))

    cats = Counter(t["category"] for t in truth_rows)
    res["truth_counts"] = dict(cats)

    # ---- per-read join
    for t in truth_rows:
        rid = t["read_id"]
        frows = finder_by.get(rid, [])
        arows = annot_by.get(rid, [])
        drows = dr_by.get(rid, [])
        fp = primary_finder(frows) if frows else None
        ap = primary_annot(arows) if arows else None
        dp = primary_dr(drows, fp) if drows else None
        rec = {"read_id": rid, "barcode": t["barcode"], "category": t["category"],
               "element": t["element"], "ref_pos": t["ref_pos"], "tsd_len": t["tsd_len"],
               "tsd_seq": t["tsd_seq"], "fragmented": t["fragmented"], "read_len": t["read_len"],
               "n_finder_rows": len(frows), "n_annot_rows": len(arows), "n_dr_rows": len(drows),
               "pred_ref_pos": fp["ref_pos"] if fp else None, "pred_len": fp["length"] if fp else None,
               "pred_type": fp["type"] if fp else None,
               "annot_note": ap["note"] if ap else None, "annot_subject": ap["subject"] if ap else None,
               "annot_notes_all": ",".join(sorted({r["note"] for r in arows})),
               "dr_point": dp["point"] if dp else None, "dr_len": dp["rep_len"] if dp else None,
               "dr_seq": dp["rep_seq"] if dp else None,
               "rep_lo": t.get("rep_lo"), "rep_hi": t.get("rep_hi"), "ext_seq": t.get("ext_seq"),
               "pos_err": None, "pos_err_tsd": None, "pos_err_raw": None,
               "class_ok": None, "class_ok_any_row": None}
        is_pos = t["category"] in INSERTION_CATS
        if is_pos and fp and fp["ref_pos"] is not None and t["ref_pos"] is not None:
            rec["pos_err"] = interval_dist(fp["ref_pos"], t["rep_lo"], t["rep_hi"], plen)
            rec["pos_err_tsd"] = interval_dist(fp["ref_pos"], t["ref_pos"], t["ref_pos"] + (t["tsd_len"] or 0), plen)
            rec["pos_err_raw"] = circ_dist(fp["ref_pos"], t["ref_pos"], plen)
        if is_pos and ap:
            if t["category"] == "is_insertion":
                rec["class_ok"] = ap["note"] == "IS_DB" and ap["subject"] == t["element"]
                rec["class_ok_any_row"] = any(r["note"] == "IS_DB" and r["subject"] == t["element"] for r in arows)
            else:
                def novel_ok(r):
                    if novel_absent:
                        return r["note"] == "Unannotated"
                    return r["note"] == "Genome" and genome_overlap(r, t)
                rec["class_ok"] = novel_ok(ap)
                rec["class_ok_any_row"] = any(novel_ok(r) for r in arows)
        rec["_ap"], rec["_arows"] = ap, arows
        per_read.append(rec)

    pos_reads = [r for r in per_read if r["category"] in INSERTION_CATS]
    neg_reads = [r for r in per_read if r["category"] in NEGATIVE_CATS]
    eligible = [r for r in pos_reads if not r["fragmented"] and (r["read_len"] or 0) >= min_read_len]

    # ---- read-level detection, insert finder and annotation stages
    def detection(called):
        if not called[1]:
            return {"available": False}
        is_called = called[0]
        tp = sum(is_called(r) for r in pos_reads)
        fp_ = sum(is_called(r) for r in neg_reads)
        fn = len(pos_reads) - tp
        p, rc = ratio(tp, tp + fp_), ratio(tp, len(pos_reads))
        tp_el = sum(is_called(r) for r in eligible)
        tp_pos = sum(is_called(r) and r["pos_err"] is not None and r["pos_err"] <= tol for r in pos_reads)
        d = {"available": True, "tp": tp, "fp": fp_, "fn": fn, "precision": p, "recall": rc, "f1": f1(p, rc),
             "recall_unfragmented_eligible": ratio(tp_el, len(eligible)), "n_unfragmented_eligible": len(eligible),
             "tp_correct_position": tp_pos, "recall_correct_position": ratio(tp_pos, len(pos_reads))}
        for cat in NEGATIVE_CATS:
            rs = [r for r in neg_reads if r["category"] == cat]
            n = sum(is_called(r) for r in rs)
            d[f"false_calls_{cat}"] = n
            d[f"false_call_rate_{cat}"] = ratio(n, len(rs))
        for cat in INSERTION_CATS:
            rs = [r for r in pos_reads if r["category"] == cat]
            d[f"recall_{cat}"] = ratio(sum(is_called(r) for r in rs), len(rs))
        return d

    res["detection_finder"] = detection((lambda r: r["n_finder_rows"] > 0, have_finder))
    res["detection_annot"] = detection(
        (lambda r: any(x["note"] in CALL_NOTES for x in r["_arows"]), have_annot))

    # ---- multiplicity (B02/B03/B12)
    if have_finder:
        called = [r for r in per_read if r["n_finder_rows"]]
        res["multiplicity_finder"] = {
            "reads_with_rows": len(called), "rows": sum(r["n_finder_rows"] for r in called),
            "pct_reads_multiple_rows": ratio(100.0 * sum(r["n_finder_rows"] > 1 for r in called), len(called)),
            "unmatched_read_ids": finder_unmatched}
    if have_annot:
        called = [r for r in per_read if r["n_annot_rows"]]
        res["multiplicity_annot"] = {
            "reads_with_rows": len(called), "rows": sum(r["n_annot_rows"] for r in called),
            "pct_reads_multiple_rows": ratio(100.0 * sum(r["n_annot_rows"] > 1 for r in called), len(called)),
            "pct_reads_contradictory_notes": ratio(100.0 * sum("," in r["annot_notes_all"] for r in called), len(called)),
            "unmatched_read_ids": annot_unmatched}

    # ---- position
    res["position"] = {
        "error": dist_summary([r["pos_err"] for r in pos_reads if r["n_finder_rows"]]),
        "tsd_interval": dist_summary([r["pos_err_tsd"] for r in pos_reads if r["n_finder_rows"]]),
        "raw": dist_summary([r["pos_err_raw"] for r in pos_reads if r["n_finder_rows"]]),
    } if have_finder else {"available": False}

    # ---- classification
    if have_annot:
        cl = {}
        for cat in INSERTION_CATS:
            rs = [r for r in pos_reads if r["n_annot_rows"]]
            rs = [r for r in rs if r["category"] == cat]
            notes = Counter(r["annot_note"] for r in rs)
            cl[cat] = {"n_annotated": len(rs),
                       "accuracy": ratio(sum(bool(r["class_ok"]) for r in rs), len(rs)),
                       "accuracy_any_row": ratio(sum(bool(r.get("class_ok_any_row")) for r in rs), len(rs)),
                       **note_pcts(notes, len(rs))}
            if cat == "is_insertion":
                wrong = Counter((r["element"], r["annot_subject"]) for r in rs
                                if r["annot_note"] == "IS_DB" and r["annot_subject"] != r["element"])
                cl[cat]["top_confusions"] = [f"{a}->{b}:{n}" for (a, b), n in wrong.most_common(5)]
        # Contamination note broken down by truth category (AUDIT B10)
        contam = Counter(r["category"] for r in per_read if r["annot_note"] == "Contamination")
        cl["contamination_note_by_truth"] = dict(contam)
        cl["contamination_note_total_reads"] = sum(contam.values())
        cl["contamination_note_true_insertions"] = sum(contam[c] for c in INSERTION_CATS)
        rs = [r for r in per_read if r["category"] == "contamination"]
        cl["contamination_reads_not_called"] = ratio(sum(r["n_finder_rows"] == 0 and r["n_annot_rows"] == 0 for r in rs), len(rs))
        cl["contamination_reads_labelled_Contamination"] = sum(r["annot_note"] == "Contamination" for r in rs)
        unann = Counter(r["category"] for r in per_read if r["annot_note"] == "Unannotated")
        cl["unannotated_note_by_truth"] = dict(unann)
        res["classification"] = cl
    else:
        res["classification"] = {"available": False}

    # ---- TSD / DR
    res["tsd"] = dict(tsd_metrics(pos_reads, plen), unmatched_read_ids=dr_unmatched) if have_dr else {"available": False}

    # ---- site level
    sites = []
    if have_finder:
        # A predicted site = reads with the same annotation label (IS name, "Genome", ...) whose
        # positions chain with gaps <= site_link bp. Using the label keeps different elements
        # inserted at nearby positions apart, as a user of the pipeline would.
        groups = defaultdict(list)
        for r in per_read:
            if r["n_finder_rows"] and r["pred_ref_pos"] is not None:
                groups[site_label(r)].append((r["pred_ref_pos"], r))
        clusters = []
        for label in sorted(groups, key=str):
            pts = sorted(groups[label], key=lambda x: x[0])
            cs, cur = [], []
            for pos, r in pts:
                if cur and pos - cur[-1][0] > site_link:
                    cs.append(cur)
                    cur = []
                cur.append((pos, r))
            if cur:
                cs.append(cur)
            # merge first/last cluster across the plasmid origin
            if plen and len(cs) > 1 and cs[0][0][0] + plen - cs[-1][-1][0] <= site_link:
                cs[0] = cs.pop() + cs[0]
            clusters += cs
        clusters = [c for c in clusters if len(c) >= min_support]
        real = [c for c in clones if c["category"] in INSERTION_CATS and (c["n_reads"] or 0) > 0]
        cand = []
        for i, c in enumerate(clusters):
            centre = int(statistics.median(p for p, _ in c))
            for j, cl in enumerate(real):
                d = interval_dist(centre, cl["rep_lo"], cl["rep_hi"], plen)
                if d <= tol:
                    lab = site_label(c[0][1])
                    same = lab == cl["element"] if cl["category"] == "is_insertion" else novel_label_ok(lab, novel_absent)
                    # closest first; ties -> matching label, then larger support
                    cand.append((d, not same, -len(c), i, j))
        used_c, used_t, match = set(), set(), {}
        for d, _, _, i, j in sorted(cand):
            if i in used_c or j in used_t:
                continue
            used_c.add(i)
            used_t.add(j)
            match[i] = (j, d)
        for i, c in enumerate(clusters):
            centre = int(statistics.median(p for p, _ in c))
            label = site_label(c[0][1])
            truth_cats = Counter(r["category"] for _, r in c)
            row = {"barcode": truth_rows[0]["barcode"] if truth_rows else "", "cluster": i, "pred_pos": centre,
                   "support": len(c), "majority_label": label, "truth_categories": dict(truth_cats),
                   "clone_id": None, "clone_pos": None, "clone_element": None, "clone_reads": None,
                   "dist": None, "label_ok": None}
            if i in match:
                cl, d = real[match[i][0]], match[i][1]
                row.update(clone_id=cl["clone_id"], clone_pos=cl["ref_pos"], clone_element=cl["element"],
                           clone_reads=cl["n_reads"], dist=d)
                if label is not None:
                    row["label_ok"] = ((label == cl["element"]) if cl["category"] == "is_insertion"
                                       else novel_label_ok(label, novel_absent))
            sites.append(row)
        matched_clones = {match[i][0] for i in match}
        big = [j for j, cl in enumerate(real) if cl["n_reads"] >= 3]
        lab = [s["label_ok"] for s in sites if s["label_ok"] is not None]
        res["sites"] = {
            "true_events": len(real), "predicted_sites": len(clusters), "matched": len(match),
            "event_recall": ratio(len(matched_clones), len(real)),
            "event_recall_ge3_reads": ratio(len(matched_clones & set(big)), len(big)),
            "site_precision": ratio(len(match), len(clusters)),
            "spurious_sites": len(clusters) - len(match),
            "spurious_sites_support_ge2": sum(1 for i, c in enumerate(clusters) if i not in match and len(c) >= 2),
            "label_accuracy_matched": ratio(sum(lab), len(lab)),
            "missed_events": [real[j]["clone_id"] for j in range(len(real)) if j not in matched_clones]}
        res["sites"]["f1"] = f1(res["sites"]["site_precision"], res["sites"]["event_recall"])
    else:
        res["sites"] = {"available": False}

    for r in per_read:
        r.pop("_ap", None)
        r.pop("_arows", None)
    return res, per_read, sites


def aggregate(per_barcode, per_read_all, sites_all, clones, resolve, plen, tol, min_read_len, min_support):
    """Overall metrics: re-score the pooled reads (counts add up; ratios recomputed)."""
    # Pool by re-running the read-level scorer on the union; sites stay per-barcode, so sum them.
    overall = {}
    for key in ("detection_finder", "detection_annot"):
        ds = [b[key] for b in per_barcode.values() if b[key].get("available")]
        if not ds:
            overall[key] = {"available": False}
            continue
        tp, fp_, fn = (sum(d[k] for d in ds) for k in ("tp", "fp", "fn"))
        p, rc = ratio(tp, tp + fp_), ratio(tp, tp + fn)
        o = {"available": True, "tp": tp, "fp": fp_, "fn": fn, "precision": p, "recall": rc, "f1": f1(p, rc)}
        n_el = sum(d["n_unfragmented_eligible"] for d in ds)
        o["recall_unfragmented_eligible"] = ratio(sum(d["recall_unfragmented_eligible"] * d["n_unfragmented_eligible"]
                                                      for d in ds if d["recall_unfragmented_eligible"] is not None), n_el)
        o["tp_correct_position"] = sum(d["tp_correct_position"] for d in ds)
        o["recall_correct_position"] = ratio(o["tp_correct_position"], tp + fn)
        for cat in NEGATIVE_CATS:
            o[f"false_calls_{cat}"] = sum(d[f"false_calls_{cat}"] for d in ds)
            n = sum(1 for r in per_read_all if r["category"] == cat)
            o[f"false_call_rate_{cat}"] = ratio(o[f"false_calls_{cat}"], n)
        for cat in INSERTION_CATS:
            rs = [r for r in per_read_all if r["category"] == cat]
            if key == "detection_finder":
                o[f"recall_{cat}"] = ratio(sum(r["n_finder_rows"] > 0 for r in rs), len(rs))
            else:
                o[f"recall_{cat}"] = ratio(sum(r["annot_note"] in CALL_NOTES or
                                               any(n in CALL_NOTES for n in r["annot_notes_all"].split(","))
                                               for r in rs), len(rs))
        overall[key] = o
    pos = [r for r in per_read_all if r["category"] in INSERTION_CATS and r["n_finder_rows"]]
    if any(b["position"].get("available", True) for b in per_barcode.values()):
        overall["position"] = {"error": dist_summary([r["pos_err"] for r in pos]),
                               "tsd_interval": dist_summary([r["pos_err_tsd"] for r in pos]),
                               "raw": dist_summary([r["pos_err_raw"] for r in pos])}
    else:
        overall["position"] = {"available": False}
    cls = [b["classification"] for b in per_barcode.values() if b["classification"].get("available", True)]
    if cls:
        o = {}
        for cat in INSERTION_CATS:
            rs = [r for r in per_read_all if r["category"] == cat and r["n_annot_rows"]]
            notes = Counter(r["annot_note"] for r in rs)
            o[cat] = {"n_annotated": len(rs),
                      "accuracy": ratio(sum(bool(r["class_ok"]) for r in rs), len(rs)),
                      "accuracy_any_row": ratio(sum(bool(r.get("class_ok_any_row")) for r in rs), len(rs)),
                      **note_pcts(notes, len(rs))}
        contam, unann = Counter(), Counter()
        for c in cls:
            contam.update(c["contamination_note_by_truth"])
            unann.update(c.get("unannotated_note_by_truth", {}))
        o["unannotated_note_by_truth"] = dict(unann)
        o["contamination_note_by_truth"] = dict(contam)
        o["contamination_note_total_reads"] = sum(contam.values())
        o["contamination_note_true_insertions"] = sum(contam[c] for c in INSERTION_CATS)
        rs = [r for r in per_read_all if r["category"] == "contamination"]
        o["contamination_reads_not_called"] = ratio(sum(r["n_finder_rows"] == 0 and r["n_annot_rows"] == 0 for r in rs), len(rs))
        overall["classification"] = o
    else:
        overall["classification"] = {"available": False}
    if any(b["tsd"].get("available") for b in per_barcode.values()):
        overall["tsd"] = tsd_metrics([r for r in per_read_all if r["category"] in INSERTION_CATS], plen)
    else:
        overall["tsd"] = {"available": False}
    ss = [b["sites"] for b in per_barcode.values() if b["sites"].get("available", True)]
    if ss:
        te, ps, m = (sum(s[k] for s in ss) for k in ("true_events", "predicted_sites", "matched"))
        o = {"true_events": te, "predicted_sites": ps, "matched": m, "event_recall": ratio(m, te),
             "site_precision": ratio(m, ps), "spurious_sites": ps - m,
             "spurious_sites_support_ge2": sum(s["spurious_sites_support_ge2"] for s in ss)}
        o["f1"] = f1(o["site_precision"], o["event_recall"])
        lab = [s["label_ok"] for s in sites_all if s["label_ok"] is not None]
        o["label_accuracy_matched"] = ratio(sum(lab), len(lab))
        overall["sites"] = o
    else:
        overall["sites"] = {"available": False}
    overall["truth_counts"] = dict(Counter(r["category"] for r in per_read_all))
    return overall


def evaluate(sim_dir, results_dir, tol=20, min_read_len=0, min_support=1, barcodes=None, refs=None, site_link=15,
             novel_absent=False):
    truth, clones, plen, have_seq = load_truth(sim_dir, refs)
    resolve = make_id_resolver(truth)
    by_bc = defaultdict(list)
    for t in truth.values():
        by_bc[t["barcode"]].append(t)
    bcs = sorted(barcodes or by_bc)
    per_barcode, per_read_all, sites_all = {}, [], []
    for bc in bcs:
        preds = load_predictions(results_dir, bc)
        res, per_read, sites = score_barcode(by_bc[bc], preds, [c for c in clones if c["barcode"] == bc],
                                             resolve, plen, tol, min_read_len, min_support, site_link,
                                             novel_absent)
        per_barcode[bc] = res
        per_read_all += per_read
        sites_all += sites
    overall = aggregate(per_barcode, per_read_all, sites_all, clones, resolve, plen, tol, min_read_len, min_support)
    meta = {"sim_dir": os.path.abspath(sim_dir), "results_dir": os.path.abspath(results_dir), "tol_bp": tol,
            "min_read_len": min_read_len, "min_support": min_support, "site_link_bp": site_link, "plasmid_len": plen,
            "novel_absent": novel_absent,
            "sequences_loaded": have_seq,
            "position_convention": ("0-based. position.error = circular distance from predicted ref_pos to the "
                                    "ambiguity window [ref_pos - left_ext, ref_pos + tsd_len + right_ext] (TSD plus "
                                    "chance extensions of the direct repeat); position.tsd_interval uses "
                                    "[ref_pos, ref_pos + tsd_len]; position.raw = |pred - ref_pos|.")}
    return {"meta": meta, "overall": overall, "per_barcode": per_barcode}, per_read_all, sites_all


# ---------------------------------------------------------------- output
def flatten(d, prefix=""):
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            yield from flatten(v, key)
        elif isinstance(v, list):
            yield key, ";".join(map(str, v))
        else:
            yield key, v


def fmt(v):
    if v is None:
        return "N/A"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def write_outputs(report, per_read, sites, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "metrics.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    with open(os.path.join(out_dir, "metrics.tsv"), "w") as fh:
        fh.write("barcode\tsection\tmetric\tvalue\n")
        scopes = [("overall", report["overall"])] + list(report["per_barcode"].items())
        for bc, d in scopes:
            for key, v in flatten({k: x for k, x in d.items() if k != "files"}):
                section, _, metric = key.partition(".")
                fh.write(f"{bc}\t{section}\t{metric or section}\t{fmt(v)}\n")
    if per_read:
        cols = list(per_read[0].keys())
        with open(os.path.join(out_dir, "per_read.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in per_read:
                w.writerow({k: fmt(r.get(k)) for k in cols})
    with open(os.path.join(out_dir, "sites.tsv"), "w", newline="") as fh:
        cols = ["barcode", "cluster", "pred_pos", "support", "majority_label", "truth_categories", "clone_id",
                "clone_pos", "clone_element", "clone_reads", "dist", "label_ok"]
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        for s in sites:
            w.writerow({k: fmt(s.get(k)) for k in cols})


SUMMARY_KEYS = [
    ("detection_finder", "precision"), ("detection_finder", "recall"), ("detection_finder", "f1"),
    ("detection_finder", "recall_unfragmented_eligible"), ("detection_finder", "false_call_rate_empty_plasmid"),
    ("detection_annot", "precision"), ("detection_annot", "recall"), ("detection_annot", "f1"),
    ("position", "error.median"), ("position", "error.pct_within_0bp"), ("position", "error.pct_within_5bp"),
    ("position", "error.pct_within_10bp"), ("position", "error.pct_within_20bp"),
    ("classification", "is_insertion.accuracy"), ("classification", "is_insertion.pct_note_Contamination"),
    ("classification", "novel_insertion.accuracy"), ("classification", "novel_insertion.pct_note_Unannotated"),
    ("classification", "contamination_note_true_insertions"),
    ("sites", "event_recall"), ("sites", "site_precision"), ("sites", "spurious_sites"),
    ("tsd", "pct_tsd_reported_nonzero"), ("tsd", "length_accuracy"), ("tsd", "sequence_accuracy"),
    ("tsd", "consistent_accuracy"), ("tsd", "consistent_accuracy_when_reported"),
]


def summary_lines(report):
    flat = {}
    for scope, d in [("overall", report["overall"])] + list(report["per_barcode"].items()):
        flat[scope] = dict(flatten({k: x for k, x in d.items() if k != "files"}))
    scopes = list(flat)
    lines = ["metric\t" + "\t".join(scopes)]
    for sec, m in SUMMARY_KEYS:
        key = f"{sec}.{m}"
        lines.append(key + "\t" + "\t".join(fmt(flat[s].get(key)) for s in scopes))
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--sim_dir", required=True, help="Simulation directory (truth.tsv, clones.tsv)")
    ap.add_argument("-r", "--results_dir", required=True, help="Pipeline data directory (barcode folders)")
    ap.add_argument("-o", "--out_dir", required=True, help="Where to write metrics.json/tsv, per_read.tsv, sites.tsv")
    ap.add_argument("--tol", type=int, default=20, help="Position tolerance (bp) for correct position / site matching (default 20)")
    ap.add_argument("--min_read_len", type=int, default=0, help="Length filter used in the run; defines the 'eligible' recall set")
    ap.add_argument("--site_link", type=int, default=15,
                    help="Reads with the same label whose positions chain within this many bp form one predicted site (default 15; > the longest TSD so left- and right-aligned calls of one event stay together)")
    ap.add_argument("--min_support", type=int, default=1, help="Minimum reads for a predicted site (default 1)")
    ap.add_argument("--barcodes", nargs="*", help="Only score these barcodes")
    ap.add_argument("-p", "--plasmid", help="Plasmid FASTA (default: path in simulation.json)")
    ap.add_argument("-i", "--is_fasta", help="IS FASTA (default: path in simulation.json)")
    ap.add_argument("-g", "--genome", help="Genome FASTA (default: path in simulation.json)")
    ap.add_argument("--novel_absent", action="store_true",
                    help="The genome given to blast_annot lacks the novel inserts (e.g. --novel_mode element "
                         "annotated with the original genome): Unannotated, not Genome, is then the correct "
                         "label for novel insertions")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)
    for f in ("truth.tsv",):
        if not os.path.exists(os.path.join(args.sim_dir, f)):
            sys.exit(f"{f} not found in {args.sim_dir}")
    report, per_read, sites = evaluate(args.sim_dir, args.results_dir, args.tol, args.min_read_len,
                                       args.min_support, args.barcodes,
                                       {"plasmid": args.plasmid, "is_fasta": args.is_fasta, "genome": args.genome}, args.site_link,
                                       args.novel_absent)
    write_outputs(report, per_read, sites, args.out_dir)
    lines = summary_lines(report)
    with open(os.path.join(args.out_dir, "summary.tsv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    if not args.quiet:
        print("\n".join(lines))
        print(f"\nWritten to {args.out_dir}")
    return report


if __name__ == "__main__":
    main()
