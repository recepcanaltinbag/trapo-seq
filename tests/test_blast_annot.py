"""Tests for src/g_blast_annot.py and src/g_blast_annot_batch.py (audit B10-B19, B52)."""
import os
import random
import shutil
import time

import pandas as pd
import pytest

from src import g_blast_annot as ga
from src.g_blast_annot import (OUTPUT_COLUMNS, classify, ensure_blast_db,
                               main_annot, parse_blast_tab)
from src.g_blast_annot_batch import main_annot_batch

needs_blast = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("makeblastdb")),
    reason="blastn/makeblastdb not in PATH")

COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(rng, s, rate):
    s = list(s)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


# ------------------------------------------------------------ pure Python ---

def hsp(subject, qs, qe, ss, se, bits=1000.0, ev=0.0, qlen=1000, slen=1000):
    return dict(subject=subject, pident=99.0, qstart=qs, qend=qe, sstart=ss,
                send=se, evalue=ev, bitscore=bits, score=int(bits), qlen=qlen,
                slen=slen)


def test_parse_blast_tab_subject_first_token_and_evalue(tmp_path):
    p = tmp_path / "b.tsv"
    p.write_text(
        "q0\tchr1\t99.5\t100\t1\t100\t500\t401\t1e-50\t180\t200\t100\t9000\tchr1 E. coli complete genome\n"
        "q0\tIS1\t90.0\t100\t1\t100\t1\t100\t1e-5\t50\t60\t100\t100\tIS1 weak hit\n")
    hits = parse_blast_tab(str(p), evalue=1e-10)
    assert list(hits) == ["q0"]
    assert len(hits["q0"]) == 1                         # 1e-5 > 1e-10 dropped
    h = hits["q0"][0]
    assert h["subject"] == "chr1"                       # B13: not "genome"
    assert (h["sstart"], h["send"]) == (500, 401)
    assert len(parse_blast_tab(str(p), evalue=1e-3)["q0"]) == 2


def test_coverage_plus_one_and_minus_strand():
    # minus-strand HSP covering the whole 1000 bp query and subject
    h = hsp("IS1", 1, 1000, 1000, 1)
    assert ga.hsp_query_len(h) == 1000 and ga.hsp_subject_len(h) == 1000
    note, top, cov = classify([h], [], threshold=70)
    assert note == "IS_DB" and cov == 100.0


def test_classify_union_of_non_overlapping_hsps_one_result():
    # two non-overlapping HSPs of IS1 (40 % + 40 %) beat a single 60 % IS2 hit
    is_hsps = [hsp("IS1", 1, 400, 1, 400, bits=700),
               hsp("IS1", 601, 1000, 601, 1000, bits=700),
               hsp("IS1", 350, 450, 350, 450, bits=100),   # overlaps, ignored
               hsp("IS2", 1, 600, 1, 600, bits=1100)]
    note, top, cov = classify(is_hsps, [hsp("chr1", 1, 1000, 1, 1000)], threshold=70)
    assert note == "IS_DB" and top["subject"] == "IS1" and cov == 80.0


def test_classify_genome_best_not_last_and_contamination():
    genome = [hsp("chr1", 1, 900, 1, 900, bits=1600), hsp("chr1", 1, 100, 5000, 5099, bits=180)]
    note, top, cov = classify([hsp("IS1", 1, 300, 1, 300)], genome, threshold=70)
    assert note == "Genome" and top["bitscore"] == 1600 and cov == 90.0   # B11
    assert classify([], [], threshold=70) == ("Unannotated", None, 0.0)
    assert classify([], [], threshold=70, flanked=False) == ("Contamination", None, 0.0)
    assert classify([], genome, flanked=False)[0] == "Genome"        # a hit always wins


def test_insert_flanked():
    f = ga.insert_flanked
    assert f("IN", 0, 900, 900)                  # CIGAR insertion: both flanks by construction
    assert f("SC", 500, 1400, 3000)              # interior uncovered stretch
    assert not f("SC", 0, 1286, 7523)            # read-start clip
    assert not f("SC", 8080, 11688, 11688)       # read-end clip
    assert not f("SC", 10, 1286, 7523)           # within READ_END_TOLERANCE of the read start
    assert f("ins", 0, 900, 900)                 # unknown type: conservative, Unannotated
    assert f("SC", 0, 900, None)                 # unknown read length


def test_partial_is_filter():
    # 300 bp of a 10 kb IS: dropped when slen > partial_len and subject cov < partial_threshold
    h = hsp("Tn1", 1, 300, 1, 300, qlen=300, slen=10000)
    assert classify([h], [], partial_threshold=80, partial_len=8000)[0] == "Unannotated"
    assert classify([h], [], partial_threshold=80, partial_len=20000)[0] == "IS_DB"


def test_no_deprecated_biopython_import():
    for mod in ("g_blast_annot.py", "g_blast_annot_batch.py"):
        src = open(os.path.join(os.path.dirname(ga.__file__), mod)).read()
        assert "Bio.Blast.Applications" not in src and "os.system" not in src


# --------------------------------------------------------- with BLAST+ ---

@pytest.fixture
def dataset(tmp_path):
    """Synthetic IS DB, genome, reads and insertion table in a dir with spaces."""
    rng = random.Random(42)
    d = tmp_path / "dir with spaces"
    d.mkdir()
    is1 = rand_seq(rng, 1300)
    is2 = rand_seq(rng, 1000)
    is_div = rand_seq(rng, 600)
    is_short = rand_seq(rng, 400)
    genome = rand_seq(rng, 20000)
    genome_with_is1 = genome[:15000] + is1 + genome[15000:]
    write_fasta(d / "IS db.fasta", [("IS3_ISAba1 transposase family IS3", is1),
                                    ("IS5_IS5 description words", is2),
                                    ("IS1_ISdiv", is_div),
                                    ("IS1_ISshort", is_short)])
    write_fasta(d / "genome ref.fasta",
                [("chr1 Escherichia coli complete genome", genome_with_is1)])
    inserts = {
        "full": is1,
        "partial": is2[:300] + genome[5000:6000],   # IS2 is not in the genome
        "none": rand_seq(rng, 900),
        "minus": revcomp(is2),
        "divergent": mutate(rng, is_div, 0.12),     # 600 bp, 12 % substitutions
        "divshort": mutate(rng, is_short, 0.12),    # 400 bp, 12 %: E-value clearly > 0
    }
    reads, rows = [], []
    for i, (name, ins) in enumerate(inserts.items()):
        left, right = rand_seq(rng, 500), rand_seq(rng, 500)
        rid = f"read{name}"
        reads.append((f"{rid} runid=x", left + ins + right))
        rows.append(f"{rid}\t500\t{500 + len(ins)}\t{len(ins)}\t{1000 + i}\t20\tins\t{i % 2 == 0}\n")
    write_fasta(d / "mapped reads.fasta", reads)
    (d / "ins.tab").write_text("".join(rows))
    return dict(dir=d, is_fa=str(d / "IS db.fasta"), genome=str(d / "genome ref.fasta"),
                reads=str(d / "mapped reads.fasta"), tab=str(d / "ins.tab"),
                inserts=inserts, is_div=is_div)


def run(ds, **kw):
    out = str(ds["dir"] / "out dir" / "best_alignment.tab")
    kw.setdefault("is_fasta", ds["is_fa"])
    main_annot(ds["tab"], ds["reads"], ds["genome"], temp_dir=str(ds["dir"] / "tmp work"),
               output_csv=out, threads=1, **kw)
    return pd.read_csv(out, sep="\t")


def by_name(df, name):
    sub = df[df["Query ID"].str.startswith(f"read{name}_")]
    assert len(sub) == 1
    return sub.iloc[0]


@needs_blast
def test_main_annot_classification(dataset):
    df = run(dataset)
    assert list(df.columns) == OUTPUT_COLUMNS
    assert df["Query ID"].is_unique and len(df) == 6          # B12
    ins = dataset["inserts"]

    full = by_name(df, "full")                                 # also in genome: IS_DB wins, one row
    assert full["Note"] == "IS_DB" and full["Subject ID"] == "IS3_ISAba1"
    assert full["Explained"] == 100.0
    assert full["Query ID"] == f"readfull_500_{500 + len(ins['full'])}"

    partial = by_name(df, "partial")
    assert partial["Note"] == "Genome" and partial["Subject ID"] == "chr1"   # B13
    span = partial["Query End"] - partial["Query Start"] + 1                 # B14: +1
    assert partial["Explained"] == round(span / 1300 * 100, 2)
    assert 1000 <= span <= 1020                                 # genome part (+ chance extension)

    none = by_name(df, "none")                                  # flanked, no hit
    assert none["Note"] == "Unannotated" and none["Subject ID"] == "no blast hit"
    assert none["Explained"] == 0

    minus = by_name(df, "minus")
    assert minus["Note"] == "IS_DB" and minus["Subject ID"] == "IS5_IS5"
    assert minus["Explained"] == 100.0
    assert minus["Subject Start"] > minus["Subject End"]
    assert minus["Subject Start"] - minus["Subject End"] + 1 == len(ins["minus"])
    assert (minus["Query Start"], minus["Query End"]) == (1, len(ins["minus"]))

    div = by_name(df, "divergent")                              # B10
    assert div["Note"] == "IS_DB" and div["Subject ID"] == "IS1_ISdiv"
    assert 80 < float(div["Identity (%)"]) < 95

    # temp files cleaned up when keep_temp is False
    assert not (dataset["dir"] / "tmp work").exists()


@needs_blast
def test_divergent_insert_with_nonzero_evalue_is_is_db(dataset):
    # B10: the old rule kept only HSPs with E == 0 (< ~1e-180). This insert's
    # best HSP has E > 0, so the old code labelled it Contamination.
    df = run(dataset, task="dc-megablast")
    row = by_name(df, "divshort")
    assert row["Note"] == "IS_DB" and row["Subject ID"] == "IS1_ISshort"
    assert 0 < float(row["E-value"]) <= 1e-10
    assert row["Explained"] >= 70


@needs_blast
def test_evalue_parameter_is_applied(dataset):
    df = run(dataset, evalue=1e-300)
    assert by_name(df, "divergent")["Note"] != "IS_DB"
    assert by_name(df, "full")["Note"] == "IS_DB"


@needs_blast
def test_keep_temp_and_dummy_is_fasta(dataset):
    df = run(dataset, is_fasta=None, keep_temp=True)
    assert set(df["Note"]) <= {"Genome", "Unannotated", "Contamination"}
    assert by_name(df, "full")["Note"] == "Genome"
    work = dataset["dir"] / "tmp work"
    assert (work / "blastdb").is_dir()
    assert any(p.name.startswith("dummy_is_file_") for p in (work / "blastdb").iterdir())


@needs_blast
def test_db_rebuilt_after_fasta_change(dataset):
    d = dataset["dir"]
    db_dir = str(d / "db dir")
    is_fa = str(d / "small IS.fasta")
    write_fasta(is_fa, [("ISx_other", rand_seq(random.Random(1), 800))])
    db = ensure_blast_db(is_fa, db_dir)
    assert os.path.basename(db).startswith("small_IS_")         # named after FASTA
    stamp = db + ".stamp"
    t0 = os.path.getmtime(stamp)
    assert ensure_blast_db(is_fa, db_dir) == db and os.path.getmtime(stamp) == t0   # cached

    df = run(dataset, is_fasta=is_fa, db_dir=db_dir)
    assert by_name(df, "divergent")["Note"] == "Unannotated"

    write_fasta(is_fa, [("IS1_ISdiv", dataset["is_div"])])
    later = time.time() + 5
    os.utime(is_fa, (later, later))
    df = run(dataset, is_fasta=is_fa, db_dir=db_dir)
    assert by_name(df, "divergent")["Note"] == "IS_DB"
    assert os.path.getmtime(stamp) != t0


@needs_blast
def test_missing_read_is_skipped(dataset):
    with open(dataset["tab"], "a") as fh:
        fh.write("readMISSING\t10\t500\t490\t5\t20\tins\tFalse\n")
    df = run(dataset)
    assert len(df) == 6


@needs_blast
def test_batch_reuses_single_and_temp_handling(dataset, tmp_path):
    data = tmp_path / "batch data"
    for sample in ("bc 01", "bc02"):
        s = data / sample
        s.mkdir(parents=True)
        shutil.copy(dataset["reads"], s / f"{sample}.fasta")
        shutil.copy(dataset["tab"], s / "insertion_from_bam.tab")
    (data / "empty").mkdir()
    res = main_annot_batch(str(data), dataset["genome"], dataset["is_fa"], threads=1)
    assert set(res) == {"bc 01", "bc02"}
    for sample in ("bc 01", "bc02"):
        df = pd.read_csv(data / sample / "best_alignment.tab", sep="\t")
        assert len(df) == 6 and df["Query ID"].is_unique
    assert not (data / "temp").exists()                          # default temp under data_dir, removed

    # explicit temp dir inside data_dir is skipped (path comparison) and kept
    res = main_annot_batch(str(data), dataset["genome"], dataset["is_fa"],
                           temp_dir=str(data / "work"), threads=1, keep_temp=True)
    assert set(res) == {"bc 01", "bc02"}
    assert (data / "work" / "blastdb").is_dir()
    res = main_annot_batch(str(data), dataset["genome"], dataset["is_fa"],
                           temp_dir=str(data / "work"), threads=1)
    assert set(res) == {"bc 01", "bc02"}


@needs_blast
def test_contamination_rule_for_read_end_clips(tmp_path):
    """SC inserts at a read end without hit are Contamination unless they match a
    both-flanked Unannotated insert of the sample; IN and interior SC are Unannotated."""
    rng = random.Random(7)
    d = tmp_path / "rule dir"
    d.mkdir()
    novel = rand_seq(rng, 1200)            # unknown element, in neither database
    junk = rand_seq(rng, 1000)             # unrelated foreign sequence
    write_fasta(d / "is.fasta", [("IS1_ISx", rand_seq(rng, 1000))])
    write_fasta(d / "genome.fasta", [("chr1", rand_seq(rng, 20000))])
    p = lambda n: rand_seq(rng, n)         # plasmid-like flank (only its length matters here)
    reads = {
        "in_novel": p(500) + mutate(rng, novel, 0.05) + p(500),          # IN
        "sc_inner": p(500) + mutate(rng, novel, 0.05) + p(500),          # interior SC
        "sc_start_novel": mutate(rng, novel[300:], 0.05) + p(3000),      # read starts inside the element
        "sc_end_novel": p(3000) + mutate(rng, revcomp(novel), 0.05),     # same element, read end, minus
        "sc_end_junk": p(3000) + junk,                                   # chimera / foreign DNA
        "sc_end_stub": p(3000) + mutate(rng, novel, 0.05) + p(900),      # element + unaligned plasmid stub
    }
    write_fasta(d / "reads.fasta", list(reads.items()))
    rows = [("in_novel", 500, 1700, "IN"), ("sc_inner", 500, 1700, "SC"),
            ("sc_start_novel", 0, 900, "SC"), ("sc_end_novel", 3000, 4200, "SC"),
            ("sc_end_junk", 3000, 4000, "SC"), ("sc_end_stub", 3000, 5100, "SC")]
    (d / "ins.tab").write_text("".join(f"{r}\t{a}\t{b}\t{b - a}\t100\t20\t{t}\tFalse\n" for r, a, b, t in rows))
    out = d / "best.tab"
    notes = main_annot(str(d / "ins.tab"), str(d / "reads.fasta"), str(d / "genome.fasta"),
                       is_fasta=str(d / "is.fasta"), temp_dir=str(d / "tmp"), output_csv=str(out), threads=1)
    got = {q.rsplit("_", 2)[0]: n for q, n in notes.items()}
    assert got == {"in_novel": "Unannotated", "sc_inner": "Unannotated", "sc_start_novel": "Unannotated",
                   "sc_end_novel": "Unannotated", "sc_end_junk": "Contamination",
                   "sc_end_stub": "Unannotated"}       # anchor covered >= threshold by the clip
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == OUTPUT_COLUMNS and len(df) == 6
    assert set(df["Subject ID"]) == {"no blast hit"}
    assert not (d / "tmp").exists()

    # without a flanked anchor, the novel read-end clips are Contamination
    (d / "ins.tab").write_text("".join(f"{r}\t{a}\t{b}\t{b - a}\t100\t20\t{t}\tFalse\n"
                                       for r, a, b, t in rows[2:5]))
    notes = main_annot(str(d / "ins.tab"), str(d / "reads.fasta"), str(d / "genome.fasta"),
                       is_fasta=str(d / "is.fasta"), temp_dir=str(d / "tmp"), output_csv=str(out), threads=1)
    assert set(notes.values()) == {"Contamination"}


def test_subject_coverage():
    h = [hsp("a", 1, 100, 1, 100, slen=400), hsp("a", 101, 200, 300, 201, slen=400),
         hsp("b", 1, 50, 1, 50, slen=50)]
    assert ga.subject_coverage(h[:2]) == 50.0
    assert ga.subject_coverage(h) == 100.0
