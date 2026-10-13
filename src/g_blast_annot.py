#-------------------------------------
# Annotate insertions with BLAST.
#
# Every insert sequence is compared with (1) a manually curated IS database and
# (2) the host genome, and classified exactly once:
#   IS_DB         - best IS hit covers >= `threshold` % of the insert
#   Genome        - otherwise, the insert has a hit in the genome
#   Unannotated   - otherwise (no IS/genome hit), and the insert is flanked by
#                   plasmid on both sides, or resembles such an insert
#   Contamination - otherwise: no hit AND evidence of a non-transposition origin
#
# CONTAMINATION RULE (data available at annotation time: insert_type and the
# read length from the mapped-reads FASTA). An insert without IS/genome hit is
#   * Unannotated if insert_type is IN: a CIGAR insertion inside one plasmid
#     alignment, so plasmid flanks both sides by construction;
#   * Unannotated if insert_type is SC and the insert is interior to the read
#     (query_start > READ_END_TOLERANCE and query_end < read length -
#     READ_END_TOLERANCE): e_insert_finder reports maximal read stretches not
#     covered by any plasmid alignment, so an interior stretch has plasmid
#     alignments on both sides;
#   * for an SC insert that reaches a read end (plasmid on one side only):
#     Unannotated if it matches a both-flanked Unannotated insert of the same
#     sample over >= `threshold` % of the clip or of that insert (BLAST, same
#     E-value; the same unknown element seen with both flanks in another read;
#     the second case is a clip holding the element plus a plasmid stub too
#     short to align, e.g. next to the linearisation site), otherwise
#     Contamination. A one-sided unknown sequence carries no evidence of an
#     insertion into the plasmid; it is typically a chimeric read (ligation of
#     a plasmid molecule to foreign DNA), adapter/barcode sequence or foreign
#     DNA;
#   * Unannotated for any other insert_type (no information: conservative).
# Limitation: two plasmid alignments around an interior SC stretch are not
# checked for adjacency/strand (a plasmid-foreign-plasmid chimera looks like an
# insertion). insertion_from_bam.tab carries the reference position of one
# flank only; a column with both flank positions/strands would allow it.
#
# PLASMID-GENOME HOMOLOGY (GitHub issue #5): insert_finder writes
# <tab stem>_homology_flags.tsv next to the tab file when the plasmid and the
# genome share segments (src/plasmid_genome_homology.py). A flagged call is a
# genomic read that only looks plasmid-flanked: every plasmid flank of its
# insert ends inside a shared segment and the read's genome alignment runs
# contiguously across the flank/insert junction. By default
# (homology_action="relabel") a flagged call is labelled Contamination (it is
# genomic DNA, not a transposition), so stats, dr_finder and novel_elements do
# not count it; the column Homology keeps the classification it would have had
# ("shared_segment:<Note>"); "-" = not flagged / not checked.
# homology_action="keep" only fills the column. Flags are matched by Query ID
# and ignored if the flags file is older than the tab file.
#
# ARTEFACT FLAGS (docs/ARTEFACTS.md, src/artefacts.py): insert_finder writes
# <tab stem>_artefacts.tsv (flank adjacency/strands, cut-site distance,
# adapters, fold-back, plasmid self-match per call). The last column Artefact
# lists the flags of each call ("-" = none). Here three flags are added:
#   * one_sided_unsupported (artefact): a one-sided call labelled Genome whose
#     genome hit overlaps no two-sided Genome call of the sample (a genome
#     fragment joined to a plasmid read end - chimera - not a transposed
#     segment seen with both flanks elsewhere);
#   * single_molecule_support (caution): the call's site (same label, positions
#     within artefacts.SITE_TOLERANCE bp) is supported by reads that could all
#     be one molecule (one read; or a duplex/fold-back/split-read pair);
#   * possible_bleed_through (caution; blast_annot_batch): the site has
#     < bleed_pct % of the reads it has in another sample.
# one_sided_cut_site is an artefact flag unless the insert is IS_DB and
# explained >= artefacts.B1_EXPLAINED % (B1: an element next to the cut site
# leaves a plasmid stub too short to align); then it is written as
# one_sided_cut_site_is (caution).
# artefact_action="flag" (default) keeps the Note; "exclude" relabels calls
# with an artefact flag (artefacts.ARTEFACT_FLAGS + one_sided_unsupported)
# Contamination and writes "excluded:<Note>:<flags>" in Artefact. Caution
# flags never change the Note. Per-site support is written to
# <output stem>_sites.tsv (artefacts.SITE_COLUMNS).
#
# NEEDS: blastn and makeblastdb in PATH, Biopython.
#
# INPUT: TAB file written by e_insert_finder_from_bam:
#   "{read_id}\t{query_start}\t{query_end}\t{insertion_length}\t{ref_pos}\t{quality}\t{insert_type}\t{is_reverse}\n"
#
# OUTPUT: TAB file, exactly one row per insertion (Query ID = {read_id}_{start}_{end}):
#   Query ID, Subject ID, Identity (%), Score, E-value, Query Start, Query End,
#   Subject Start, Subject End, Note, Explained, ref_pos, is_Reverse, Homology,
#   Artefact
#   - Subject ID is the first token of the database FASTA header.
#   - Identity/Score/E-value/coordinates describe the best (highest bitscore)
#     HSP of the chosen subject; Score is the raw BLAST score.
#   - Explained is the % of the insert covered by the union of non-overlapping
#     HSPs of the chosen subject (0 for Unannotated and Contamination, whose
#     Subject ID is "no blast hit").
#
# Recep Can Altınbağ, 07 11 2024, v0.2; v0.3 rewritten (audit B10-B19, B52)
#-------------------------------------

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO

from . import artefacts as art
from . import plasmid_genome_homology as homology

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DUMMY_IS_FASTA = os.path.join(MODULE_DIR, "dummy_is_file.fasta")

DEFAULT_EVALUE = 1e-10

NOTES = ("IS_DB", "Genome", "Unannotated", "Contamination")
HOMOLOGY_ACTIONS = ("relabel", "keep")
HOMOLOGY_PREFIX = "shared_segment:"
ARTEFACT_ACTIONS = ("flag", "exclude")
ARTEFACT_EXCLUDED_PREFIX = "excluded:"
# flags added here (see ARTEFACT FLAGS)
EXTRA_ARTEFACT_FLAGS = ("one_sided_unsupported",)
SITES_SUFFIX = "_sites.tsv"

# An SC insert starting/ending within this many bp of a read end is a read-end
# clip (e_insert_finder reports such gaps starting at 0 / ending at the read length).
READ_END_TOLERANCE = 20

OUTPUT_COLUMNS = ["Query ID", "Subject ID", "Identity (%)", "Score", "E-value",
                  "Query Start", "Query End", "Subject Start", "Subject End",
                  "Note", "Explained", "ref_pos", "is_Reverse", "Homology", "Artefact"]

# Tabular BLAST columns; stitle must stay last (it may contain spaces).
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore", "score", "qlen",
                 "slen", "stitle"]


#----------------------------------------------------------------- input ----

def read_tab_file(file_path):
    tab_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            read_id = fields[0]
            query_start = int(fields[1])
            query_end = int(fields[2])
            insertion_length = int(fields[3])
            ref_pos = fields[4]
            quality = int(float(fields[5]))
            insert_type = fields[6]
            is_reverse = fields[7].strip() == "True"
            tab_data.append((read_id, query_start, query_end, insertion_length,
                             ref_pos, quality, insert_type, is_reverse))
    return tab_data


#----------------------------------------------------------- subprocess -----

def _run(cmd, what, cwd=None):
    """Run an external command (arg list, no shell); raise a clear error."""
    try:
        return subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise RuntimeError(f"{what}: executable '{cmd[0]}' not found in PATH. "
                           f"Install BLAST+ (e.g. conda install -c bioconda blast).") from None
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{what} failed (exit code {e.returncode}).\n"
                           f"Command: {' '.join(cmd)}\n"
                           f"stderr:\n{e.stderr.strip()}") from None


#------------------------------------------------------------ BLAST DBs -----

def _fasta_stamp(fasta):
    st = os.stat(fasta)
    return f"{os.path.abspath(fasta)}\t{st.st_mtime_ns}\t{st.st_size}\n"


def ensure_blast_db(fasta, db_dir):
    """Return the path prefix of a nucleotide BLAST DB for `fasta`, building it
    in `db_dir` when missing or when the FASTA changed since the last build.

    The DB is named after the FASTA basename plus a short hash of its absolute
    path, so two different FASTAs never share a DB. A stamp file records the
    FASTA path, mtime and size used for the build (B17).
    """
    if not os.path.isfile(fasta):
        raise FileNotFoundError(f"FASTA file not found: {fasta}")
    os.makedirs(db_dir, exist_ok=True)
    abs_fasta = os.path.abspath(fasta)
    base = os.path.splitext(os.path.basename(abs_fasta))[0]
    base = "".join(c if (c.isalnum() or c in "._-") else "_" for c in base)
    tag = hashlib.sha1(abs_fasta.encode()).hexdigest()[:8]
    db = os.path.join(db_dir, f"{base}_{tag}")
    stamp_path = db + ".stamp"

    current = _fasta_stamp(fasta)
    up_to_date = False
    if os.path.isfile(stamp_path):
        with open(stamp_path) as fh:
            old = fh.read()
        up_to_date = (old == current and
                      os.path.getmtime(fasta) <= os.path.getmtime(stamp_path))
    if not up_to_date:
        _build_db(abs_fasta, db_dir, os.path.basename(db))
        with open(stamp_path, "w") as fh:
            fh.write(current)
    return db


def _build_db(abs_fasta, db_dir, name):
    """makeblastdb that works when paths contain spaces.

    makeblastdb splits -in/-out on whitespace and fails a post-build check
    when the absolute output path contains a space. So the FASTA is streamed
    on stdin, the DB is built in a fresh system temp dir, and the files are
    then moved into db_dir (blastn is later run with cwd=db_dir, -db <name>).
    """
    build_dir = tempfile.mkdtemp(prefix="trapo_blastdb_")
    try:
        if any(c.isspace() for c in build_dir):
            raise RuntimeError(f"System temp dir '{build_dir}' contains whitespace; "
                               f"makeblastdb cannot write there. Set TMPDIR to a path without spaces.")
        with open(abs_fasta, "rb") as fin:
            try:
                subprocess.run(["makeblastdb", "-in", "-", "-dbtype", "nucl",
                                "-title", name, "-out", name],
                               cwd=build_dir, stdin=fin, check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except FileNotFoundError:
                raise RuntimeError("makeblastdb not found in PATH. Install BLAST+ "
                                   "(e.g. conda install -c bioconda blast).") from None
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"makeblastdb failed for {abs_fasta} (exit code {e.returncode}).\n"
                                   f"stderr:\n{(e.stderr or e.stdout).strip()}") from None
        for f in os.listdir(db_dir):
            if f.startswith(name + ".") and not f.endswith(".stamp"):
                os.remove(os.path.join(db_dir, f))
        for f in os.listdir(build_dir):
            shutil.move(os.path.join(build_dir, f), os.path.join(db_dir, f))
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


#---------------------------------------------------------------- BLAST -----

def run_blast(query_fasta, db, out_tab, evalue=DEFAULT_EVALUE, threads=1,
              task="megablast"):
    """One blastn call for a multi-FASTA query, tabular output (BLAST_COLUMNS).

    `db` is the prefix returned by ensure_blast_db. blastn runs inside the DB
    directory with a bare DB name because -db is split on whitespace.
    """
    _run(["blastn", "-task", task,
          "-query", os.path.abspath(query_fasta),
          "-db", os.path.basename(db),
          "-out", os.path.abspath(out_tab), "-evalue", str(evalue),
          "-num_threads", str(max(1, int(threads))),
          "-outfmt", "6 " + " ".join(BLAST_COLUMNS)],
         f"blastn against {os.path.basename(db)}", cwd=os.path.dirname(os.path.abspath(db)))
    return out_tab


def parse_blast_tab(path, evalue=DEFAULT_EVALUE):
    """Parse tabular BLAST output into {qseqid: [hsp dict, ...]}.

    HSPs with E-value > `evalue` are dropped (B10: no more E == 0 rule).
    Subject ID is the first token of the FASTA header (B13).
    """
    hits = {}
    if not os.path.isfile(path):
        return hits
    n = len(BLAST_COLUMNS)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t", n - 1)
            if len(f) < n - 1:
                continue
            f += [""] * (n - len(f))
            row = dict(zip(BLAST_COLUMNS, f))
            ev = float(row["evalue"])
            if ev > evalue:
                continue
            title = row["stitle"].strip()
            subject = title.split()[0] if title else row["sseqid"]
            hits.setdefault(row["qseqid"], []).append({
                "subject": subject,
                "pident": float(row["pident"]),
                "qstart": int(row["qstart"]),
                "qend": int(row["qend"]),
                "sstart": int(row["sstart"]),
                "send": int(row["send"]),
                "evalue": ev,
                "bitscore": float(row["bitscore"]),
                "score": int(float(row["score"])),
                "qlen": int(row["qlen"]),
                "slen": int(row["slen"]),
            })
    return hits


#------------------------------------------------------- classification -----

def hsp_query_len(h):
    """Bases of the query covered by an HSP (inclusive coordinates, B14)."""
    return abs(h["qend"] - h["qstart"]) + 1


def hsp_subject_len(h):
    """Bases of the subject covered by an HSP; works on both strands (B14)."""
    return abs(h["send"] - h["sstart"]) + 1


def _overlap(a, b):
    a1, a2 = sorted((a["qstart"], a["qend"]))
    b1, b2 = sorted((b["qstart"], b["qend"]))
    return a1 <= b2 and b1 <= a2


def best_subject(hsps):
    """Pick the best subject among HSPs of one query.

    For each subject, HSPs are taken greedily by bitscore, keeping only those
    that do not overlap (on the query) an already kept one. Subjects are ranked
    by covered query bases, then by summed bitscore.

    Returns (top_hsp, coverage_percent) or (None, 0.0).
    """
    by_subject = {}
    for h in hsps:
        by_subject.setdefault(h["subject"], []).append(h)
    best = None
    for subject, hs in by_subject.items():
        hs = sorted(hs, key=lambda h: (-h["bitscore"], h["evalue"]))
        kept = []
        for h in hs:
            if not any(_overlap(h, k) for k in kept):
                kept.append(h)
        covered = sum(hsp_query_len(k) for k in kept)
        bits = sum(k["bitscore"] for k in kept)
        key = (covered, bits)
        if best is None or key > best[0]:
            best = (key, kept[0])
    if best is None:
        return None, 0.0
    top = best[1]
    covered = best[0][0]
    return top, min(100.0, covered / top["qlen"] * 100)


def filter_partial_is(hsps, partial_threshold, partial_len):
    """Drop IS HSPs covering < partial_threshold % of a long (> partial_len) IS."""
    out = []
    for h in hsps:
        if (h["slen"] > partial_len and
                hsp_subject_len(h) / h["slen"] * 100 < partial_threshold):
            continue
        out.append(h)
    return out


def classify(is_hsps, genome_hsps, threshold=70, partial_threshold=80,
             partial_len=8000, flanked=True):
    """Classify one insertion. Returns (note, hsp_or_None, explained_percent).

    ``flanked``: whether the insert has plasmid on both sides (see
    insert_flanked); decides Unannotated vs Contamination when nothing hits.
    """
    is_hsps = filter_partial_is(is_hsps or [], partial_threshold, partial_len)
    top, cov = best_subject(is_hsps)
    if top is not None and cov >= threshold:
        return "IS_DB", top, cov
    top, cov = best_subject(genome_hsps or [])
    if top is not None:
        return "Genome", top, cov
    return ("Unannotated" if flanked else "Contamination"), None, 0.0


def insert_flanked(insert_type, query_start, query_end, read_len,
                   tolerance=READ_END_TOLERANCE):
    """True if the insert has plasmid alignments on both sides (or nothing
    says otherwise); False only for an SC insert reaching a read end."""
    if str(insert_type).strip().upper() != "SC" or read_len is None:
        return True
    return query_start > tolerance and query_end < read_len - tolerance


def _format_row(query_id, note, hsp, explained, ref_pos, is_reverse, homology_value="-",
                artefact_value="-"):
    if hsp is None:
        vals = [query_id, "no blast hit", "N/A", "N/A", "N/A", 0, 0, 0, 0]
    else:
        vals = [query_id, hsp["subject"], hsp["pident"], hsp["score"],
                hsp["evalue"], hsp["qstart"], hsp["qend"], hsp["sstart"],
                hsp["send"]]
    vals += [note, round(explained, 2), ref_pos, is_reverse, homology_value, artefact_value]
    return "\t".join(str(v) for v in vals) + "\n"


def rescue_end_clips(query_fa, clip_names, anchor_names, work_dir, threshold=70,
                     evalue=DEFAULT_EVALUE, threads=1, task="megablast"):
    """Names of read-end clips (``clip_names``) that match a both-flanked
    Unannotated insert (``anchor_names``) over >= threshold % of the clip or
    of the anchor (a clip holding the whole element plus an unaligned plasmid
    stub, e.g. an insertion next to the linearisation site).

    Both name lists refer to records of ``query_fa``. Returns a set of names.
    """
    os.makedirs(work_dir, exist_ok=True)
    clip_set, anchor_set = set(clip_names), set(anchor_names)
    clips_fa = os.path.join(work_dir, "end_clips.fasta")
    anchors_fa = os.path.join(work_dir, "flanked_unannotated.fasta")
    with open(clips_fa, "w") as fc, open(anchors_fa, "w") as fa:
        for rec in SeqIO.parse(query_fa, "fasta"):
            if rec.id in clip_set:
                fc.write(f">{rec.id}\n{rec.seq}\n")
            elif rec.id in anchor_set:
                fa.write(f">{rec.id}\n{rec.seq}\n")
    db = ensure_blast_db(anchors_fa, work_dir)
    out = os.path.join(work_dir, "end_clips_vs_unannotated.tsv")
    run_blast(clips_fa, db, out, evalue, threads, task)
    rescued = set()
    for qname, hsps in parse_blast_tab(out, evalue).items():
        top, cov = best_subject(hsps)
        if top is not None and (cov >= threshold or subject_coverage(hsps) >= threshold):
            rescued.add(qname)
    return rescued


def subject_coverage(hsps):
    """Best % of a subject covered by the union of its HSPs (on the subject)."""
    best = 0.0
    by_subject = {}
    for h in hsps:
        by_subject.setdefault(h["subject"], []).append(sorted((h["sstart"], h["send"])))
    for subject, ivs in by_subject.items():
        slen = next(h["slen"] for h in hsps if h["subject"] == subject)
        covered, end = 0, 0
        for a, b in sorted(ivs):
            a = max(a, end + 1)
            if b >= a:
                covered += b - a + 1
            end = max(end, b)
        best = max(best, 100.0 * covered / slen)
    return best


#------------------------------------------------------------- homology -----

def load_homology_flags(tab_file, homology_flags="auto"):
    """Query IDs flagged as plasmid-genome homology artefacts (empty set if
    there is no usable flags file)."""
    if homology_flags is None:
        return set()
    path = homology.flags_path_for(tab_file) if homology_flags == "auto" else homology_flags
    if not os.path.isfile(path):
        if homology_flags != "auto":
            raise FileNotFoundError(f"blast_annot: homology flags file not found: {path}")
        return set()
    if os.path.getmtime(path) < os.path.getmtime(tab_file):
        print(f"Warning: {path} is older than {tab_file}; homology flags ignored "
              f"(re-run insert_finder)", file=sys.stderr)
        return set()
    return homology.read_flagged_ids(path)


#------------------------------------------------------------- artefacts ----

def load_artefacts(tab_file, artefacts_tsv="auto"):
    """{Query ID: artefacts row} from <tab stem>_artefacts.tsv ({} if absent,
    off (None) or older than the tab file)."""
    if artefacts_tsv is None:
        return {}
    path = art.artefacts_path_for(tab_file) if artefacts_tsv == "auto" else artefacts_tsv
    if not os.path.isfile(path):
        if artefacts_tsv != "auto":
            raise FileNotFoundError(f"blast_annot: artefacts table not found: {path}")
        return {}
    if os.path.getmtime(path) < os.path.getmtime(tab_file):
        print(f"Warning: {path} is older than {tab_file}; artefact flags ignored "
              f"(re-run insert_finder)", file=sys.stderr)
        return {}
    return art.read_artefacts(path)


def _plasmid_len(tab_file):
    path = art.cut_site_path_for(tab_file)
    if os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                k, _, v = line.rstrip("\n").partition("\t")
                if k == "plasmid_len" and v.strip():
                    try:
                        return int(v)
                    except ValueError:
                        return None
    return None


def _genome_interval(hsp):
    lo, hi = sorted((hsp["sstart"], hsp["send"]))
    return hsp["subject"], lo, hi


def artefact_flags_for(calls, evidence, plasmid_len=None):
    """Final flags per call (see ARTEFACT FLAGS).

    calls: dicts with qname, query_id, read_id, note, hsp, explained, ref_pos.
    evidence: {Query ID: artefacts row}. Returns ({qname: [flags]}, site rows)."""
    flags = {}
    two_sided_genome = []
    for c in calls:
        ev = evidence.get(c["query_id"])
        a, ca = art.split_flags(ev["flags"]) if ev else ([], [])
        fl = a + ca
        if "one_sided_cut_site" in fl and c["note"] == "IS_DB" and c["explained"] >= art.B1_EXPLAINED:
            fl[fl.index("one_sided_cut_site")] = "one_sided_cut_site_is"
        flags[c["qname"]] = fl
        if (c["note"] == "Genome" and c["hsp"] is not None and ev and ev.get("flanks") == "both"
                and not a):
            two_sided_genome.append(_genome_interval(c["hsp"]))
    for c in calls:
        fl = flags[c["qname"]]
        if c["note"] == "Genome" and c["hsp"] is not None and "one_sided" in fl:
            subj, lo, hi = _genome_interval(c["hsp"])
            if not any(s == subj and lo <= h and l <= hi for s, l, h in two_sided_genome):
                fl.append("one_sided_unsupported")
    # sites: biological calls without artefact flags count as support
    site_calls = []
    for c in calls:
        if c["note"] not in ("IS_DB", "Genome", "Unannotated"):
            continue
        ev = evidence.get(c["query_id"]) or {}
        label = c["hsp"]["subject"] if c["note"] in ("IS_DB", "Genome") and c["hsp"] else c["note"]
        pos = art._to_int(ev.get("junction"))
        if pos is None:
            pos = art._to_int(c["ref_pos"])
        site_calls.append({"query_id": c["query_id"], "read_id": c["read_id"], "label": label,
                           "pos": pos, "strand": ev.get("read_strand") or None,
                           "start": ev.get("read_plasmid_start"), "end": ev.get("read_plasmid_end"),
                           "artefact": is_artefact(flags[c["qname"]])})
    sites, member = art.site_table(site_calls, plasmid_len)
    single = {r["site_id"] for r in sites if r["single_molecule_support"]}
    for c in calls:
        sid = member.get(c["query_id"])
        if sid in single and not is_artefact(flags[c["qname"]]):
            flags[c["qname"]].append("single_molecule_support")
    return flags, sites


def is_artefact(flags):
    return any(f in art.ARTEFACT_FLAGS or f in EXTRA_ARTEFACT_FLAGS for f in flags)


def write_sites(sites, path):
    art.write_rows(sites, path, art.SITE_COLUMNS)


def sites_path_for(output_csv):
    return os.path.splitext(output_csv)[0] + SITES_SUFFIX


#----------------------------------------------------------------- main -----

def _remove(paths):
    for p in paths:
        try:
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.remove(p)
        except OSError:
            pass


def main_annot(tab_file, mapped_fasta, genom_fasta, is_fasta=None,
               temp_dir="temp", output_csv="best_alignment.tab", threshold=70,
               threads=2, debug=False, keep_temp=False, partial_threshold=80,
               partial_len=8000, evalue=DEFAULT_EVALUE, db_dir=None,
               task="megablast", homology_flags="auto", homology_action="relabel",
               artefacts_tsv="auto", artefact_action="flag"):
    """Annotate all insertions of one sample. Writes `output_csv` and returns
    a dict {Query ID: Note}.

    homology_flags: flags table of plasmid-genome homology artefacts ("auto":
    <tab stem>_homology_flags.tsv next to tab_file if present; None: off).
    homology_action: "relabel" (flagged calls -> Contamination) or "keep"
    (Note unchanged); both fill the Homology column (see PLASMID-GENOME HOMOLOGY).
    artefacts_tsv: artefact evidence of insert_finder ("auto": <tab stem>_artefacts.tsv
    next to tab_file if present; None: off). artefact_action: "flag" (Note kept,
    flags in column Artefact) or "exclude" (calls with an artefact flag are
    relabelled Contamination); see ARTEFACT FLAGS. Writes <output stem>_sites.tsv.

    is_fasta=None uses src/dummy_is_file.fasta. BLAST DBs are built in
    `db_dir` (default: <temp_dir>/blastdb). Unless keep_temp, the files this
    function creates in temp_dir are removed (temp_dir itself only if empty).
    """
    if is_fasta is None:
        is_fasta = DUMMY_IS_FASTA
    for label, path in (("tab file", tab_file), ("mapped FASTA", mapped_fasta),
                        ("genome FASTA", genom_fasta), ("IS FASTA", is_fasta)):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"blast_annot: {label} not found: {path}")
    if homology_action not in HOMOLOGY_ACTIONS:
        raise ValueError(f"homology_action must be one of {HOMOLOGY_ACTIONS}, not {homology_action!r}")
    if artefact_action not in ARTEFACT_ACTIONS:
        raise ValueError(f"artefact_action must be one of {ARTEFACT_ACTIONS}, not {artefact_action!r}")
    flagged = load_homology_flags(tab_file, homology_flags)
    evidence = load_artefacts(tab_file, artefacts_tsv)
    print(f'Settings: Threshold: {threshold}, Threads: {threads}, E-value: {evalue}, '
          f'Debugging: {debug}, Partial Threshold: {partial_threshold}, Partial Len: {partial_len}')

    temp_existed = os.path.isdir(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)
    own_db_dir = db_dir is None
    if own_db_dir:
        db_dir = os.path.join(temp_dir, "blastdb")
    created = [db_dir] if own_db_dir else []

    # One multi-FASTA of all inserts (B19). Reads are indexed once.
    tab_data = read_tab_file(tab_file)
    print('Blast Process of ', len(tab_data))
    tag = os.path.splitext(os.path.basename(output_csv))[0]
    query_fa = os.path.join(temp_dir, f"{tag}_inserts.fasta")
    is_out = os.path.join(temp_dir, f"{tag}_is_blast.tsv")
    genome_out = os.path.join(temp_dir, f"{tag}_genome_blast.tsv")
    created += [query_fa, is_out, genome_out]

    insertions = []   # (qname, query_id, ref_pos, is_reverse, qual, flanked)
    seen = set()
    reads = SeqIO.index(mapped_fasta, "fasta")
    try:
        with open(query_fa, "w") as fh:
            for read_id, qs, qe, _qlen, ref_pos, qual, itype, is_reverse in tab_data:
                query_id = f"{read_id}_{qs}_{qe}"
                if query_id in seen:
                    continue
                seen.add(query_id)
                if read_id not in reads:
                    print(f"Warning: read {read_id} not found in {mapped_fasta}; skipping {query_id}",
                          file=sys.stderr)
                    continue
                read_seq = reads[read_id].seq
                seq = str(read_seq[qs:qe])
                if not seq:
                    print(f"Warning: empty insert sequence for {query_id}; skipping", file=sys.stderr)
                    continue
                qname = f"q{len(insertions)}"
                fh.write(f">{qname}\n{seq}\n")
                flanked = insert_flanked(itype, qs, qe, len(read_seq))
                insertions.append((qname, query_id, ref_pos, is_reverse, qual, flanked))
    finally:
        reads.close()

    is_hits, genome_hits = {}, {}
    if insertions:
        is_db = ensure_blast_db(is_fasta, db_dir)
        genome_db = ensure_blast_db(genom_fasta, db_dir)
        run_blast(query_fa, is_db, is_out, evalue, threads, task)
        run_blast(query_fa, genome_db, genome_out, evalue, threads, task)
        is_hits = parse_blast_tab(is_out, evalue)
        genome_hits = parse_blast_tab(genome_out, evalue)
    print('Blast Processes were Ended')

    classified = []
    for qname, query_id, ref_pos, is_reverse, qual, flanked in insertions:
        note, hsp, explained = classify(is_hits.get(qname), genome_hits.get(qname),
                                        threshold, partial_threshold, partial_len, flanked)
        classified.append([qname, query_id, note, hsp, explained, ref_pos, is_reverse, qual])

    # Read-end clips without hit: rescued as Unannotated when they match a
    # both-flanked Unannotated insert of this sample (see CONTAMINATION RULE).
    rescued = set()
    is_flagged = {c[0] for c in classified if c[1] in flagged}
    contam = [c[0] for c in classified if c[2] == "Contamination" and c[0] not in is_flagged]
    anchors = [c[0] for c in classified if c[2] == "Unannotated" and c[0] not in is_flagged]
    if artefact_action == "exclude":
        # chimeric inserts must not vouch for read-end clips
        bad = {c[0] for c in classified
               if art.split_flags((evidence.get(c[1]) or {}).get("flags"))[0]}
        anchors = [q for q in anchors if q not in bad]
    if contam and anchors:
        rescue_dir = os.path.join(temp_dir, f"{tag}_unannotated_db")
        created.append(rescue_dir)
        rescued = rescue_end_clips(query_fa, contam, anchors, rescue_dir, threshold,
                                   evalue, threads, task)
    notes = {}
    rows = []
    n_flagged = 0
    final = []
    for qname, query_id, note, hsp, explained, ref_pos, is_reverse, qual in classified:
        if qname in rescued:
            note = "Unannotated"
        homology_value = "-"
        if qname in is_flagged:
            n_flagged += 1
            homology_value = HOMOLOGY_PREFIX + note
            if homology_action == "relabel":
                note = "Contamination"
        final.append({"qname": qname, "query_id": query_id, "read_id": query_id.rsplit("_", 2)[0],
                      "note": note, "hsp": hsp, "explained": explained, "ref_pos": ref_pos,
                      "is_reverse": is_reverse, "qual": qual, "homology": homology_value})
    art_flags, sites = artefact_flags_for(final, evidence, _plasmid_len(tab_file))
    n_art, n_excluded = 0, 0
    for c in final:
        note, fl = c["note"], art_flags.get(c["qname"], [])
        artefact_value = ",".join(fl) or "-"
        if is_artefact(fl):
            n_art += 1
            if artefact_action == "exclude":
                artefact_value = f"{ARTEFACT_EXCLUDED_PREFIX}{note}:{artefact_value}"
                if note != "Contamination":
                    n_excluded += 1
                note = "Contamination"
        if debug:
            print(c["query_id"], note, f"explained={c['explained']:.2f}",
                  c["hsp"]["subject"] if c["hsp"] else "", f"quality={c['qual']}",
                  "(read-end clip matching a flanked Unannotated insert)" if c["qname"] in rescued else "",
                  artefact_value)
        notes[c["query_id"]] = note
        rows.append(_format_row(c["query_id"], note, c["hsp"], c["explained"], c["ref_pos"],
                                c["is_reverse"], c["homology"], artefact_value))

    out_dir = os.path.dirname(os.path.abspath(output_csv))
    os.makedirs(out_dir, exist_ok=True)
    with open(output_csv, 'w') as out_file:
        out_file.write("\t".join(OUTPUT_COLUMNS) + "\n")
        out_file.write(''.join(rows))
    write_sites(sites, sites_path_for(output_csv))

    if not keep_temp:
        _remove(created)
        if not temp_existed:
            try:
                os.rmdir(temp_dir)
            except OSError:
                pass

    counts = {n: list(notes.values()).count(n) for n in NOTES}
    print(f"Classified: {counts}" + (f" ({len(rescued)} read-end clip(s) Unannotated by similarity "
                                      f"to flanked Unannotated inserts)" if rescued else ""))
    if n_flagged:
        print(f"Plasmid-genome homology: {n_flagged} call(s) are genomic reads that only look "
              f"plasmid-flanked through a shared segment" +
              (" and were labelled Contamination" if homology_action == "relabel" else
               " (Note kept, see column Homology)"))
    if evidence:
        print(f"Artefact flags: {n_art} call(s) carry an artefact flag" +
              (f"; {n_excluded} relabelled Contamination (--artefact_action exclude)"
               if artefact_action == "exclude" else " (Note kept, see column Artefact)"))
    else:
        print("Artefact flags: no artefacts table from insert_finder (re-run it to get flank/chimera checks)")
    print(f"Outputs {output_csv} were recorded.")
    return notes
