"""Tests for scripts/build_element_db.py on tiny synthetic inputs (no IS_Db/ needed)."""
import importlib.util
import os
import random

import pytest

from conftest import REPO
from src.h_is_stats import split_subject_id

_spec = importlib.util.spec_from_file_location(
    "build_element_db", os.path.join(REPO, "scripts", "build_element_db.py"))
bdb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdb)


def rseq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def read_fasta(path):
    recs, h = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            h = line[1:]
            recs[h] = []
        elif h is not None:
            recs[h].append(line)
    return {k: "".join(v) for k, v in recs.items()}


def read_tsv(path):
    lines = open(path).read().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, l.split("\t"))) for l in lines[1:]]


S1, S2, S3, S4, S5 = (rseq(300 + 10 * i, i) for i in range(5))


@pytest.fixture
def inputs(tmp_path):
    """ISfinder FASTA + taxonomy + a TnCentral-style file."""
    isf = tmp_path / "IS.fna"
    isf.write_text(
        f">ISAb1_IS200/IS605_IS1341\n{wrap(S1)}\n"
        f">ISAb1_IS200/IS605_IS1341\n{wrap(S1)}\n"          # duplicated ID, identical (as in IS.fna)
        f">MITEAb1_IS4_IS4\n{wrap(S2[:150])}\n"
        f">ISBj2_B_IS5_unknown\n{wrap(S3)}\n"                # name with underscore, group unknown
        f">IS26_IS6_IS6\n{wrap(S4)}\n")
    tax = tmp_path / "isfinder_taxonomy.tsv"
    tax.write_text("IS_Name\tIS_Family\tIS_Group\nISAb1\tIS200/IS605\tIS1341\n"
                   "MITEAb1\tIS4\tIS4\nISBj2_B\tIS5\tIS903\nIS26\tIS6\tIS6\nIS10R\tIS4\tIS10\n")
    tn = tmp_path / "tncentral_test.fa"
    tn.write_text(
        f">ISFinder_ISAb1\n{bdb.revcomp(S1)}\n"              # RC duplicate of an ISfinder record
        f">IS26.var1-GU562437.2\n{S4[:-1] + ('A' if S4[-1] != 'A' else 'C')}\n"   # variant
        f">IS10_p-MH626558\n{rseq(200, 50)}\n"
        f">Tn21-AF071413\n{rseq(9000, 51)}\n"
        f">ISFinder_TnAb1\n{rseq(4000, 52)}\n"
        f">In2-AF071413\n{rseq(3000, 53)}\n"
        f">INTEGRALL-blaOXA-1-AF255921\n{rseq(500, 54).lower()}\n"
        f">IS1326_IS1353-AF071413\n{rseq(4000, 55)}\n"
        f">ISCR2-CP117330\n{rseq(1800, 56)}\n"
        f">IS26-X00011\n{rseq(820, 57)}\n")                  # same name as ISfinder IS26, other seq
    return isf, tax, tn


def build(tmp_path, isf, tax, extras, types="default", **kw):
    kept, report = bdb.build(str(isf), str(tax), [str(e) for e in extras], **kw)
    out = tmp_path / "out"
    bdb.write_outputs(kept, report, str(out), bdb.parse_types(types))
    return kept, report, out


#----------------------------------------------------------- malformed ------

def test_repair_wrapped_glued_headers(tmp_path):
    """Headers glued to the end of a sequence line, wrapped over two lines and
    glued to the next sequence (like _tncentral_nointegrall_isfinder-TNs.fasta)."""
    a, b, c = rseq(150, 1).lower(), rseq(130, 2).lower(), rseq(170, 3).lower()
    stream = a + ">In1223-KX784502" + b + ">Tn6023-GU562437.2" + c
    # first header is clean, the stream is wrapped at 60 columns
    text = ">ISX1-AB000001\n" + wrap(stream) + "\n>ISY1-AB000002\n" + wrap(rseq(90, 4)) + "\n"
    f = tmp_path / "bad.fasta"
    f.write_text(text)
    recs, problems = bdb.read_fasta_repair(str(f))
    assert problems == []
    got = {r.header: (r.seq, r.repaired) for r in recs}
    assert got["ISX1-AB000001"] == (a, False)
    assert got["In1223-KX784502"] == (b, True)
    assert got["Tn6023-GU562437.2"] == (c, True)
    assert got["ISY1-AB000002"][0] == rseq(90, 4)


def test_repair_header_split_across_lines(tmp_path):
    """'...acgt>In_Tn6234-HG9' / '34082tgtc...' (header wrapped mid-accession)."""
    a, b = rseq(44, 5).lower(), rseq(200, 6).lower()
    width = 58
    glued = a + ">In_Tn6234-HG9"                       # exactly one full-width line
    assert len(glued) == width
    rest = "34082" + b
    body = [rseq(width, 9).lower(), rseq(width, 10).lower(), glued]
    body += [rest[i:i + width] for i in range(0, len(rest), width)]
    f = tmp_path / "bad2.fasta"
    f.write_text(">Tn1-X00001\n" + "\n".join(body) + "\n")
    recs, problems = bdb.read_fasta_repair(str(f))
    assert problems == []
    assert [r.header for r in recs] == ["Tn1-X00001", "In_Tn6234-HG934082"]
    assert recs[0].seq == rseq(58, 9).lower() + rseq(58, 10).lower() + a
    assert recs[1].seq == b and recs[1].repaired


def test_repair_single_line_fasta(tmp_path):
    """Single-line FASTA: the glued header ends the line, sequence on the next line
    (like tncentral_integrall_isfinder.fa); an uppercase next line stays sequence."""
    a, b = rseq(100, 11).lower(), rseq(100, 12)
    f = tmp_path / "bad3.fa"
    f.write_text(f">ISA-AB000001\n{a}>In992-MF344578\n{b}\n>ISB-AB000003\n{rseq(80, 13)}\n")
    recs, problems = bdb.read_fasta_repair(str(f))
    assert problems == []
    assert [(r.header, r.seq) for r in recs][:2] == [("ISA-AB000001", a), ("In992-MF344578", b)]


def test_unrepairable_record_reported(tmp_path):
    f = tmp_path / "bad4.fa"
    f.write_text(f">ISA-AB000001\n{rseq(60, 14).lower()}>no header here!\n{rseq(60, 15)}\n"
                 f">ISB-AB000003\n{rseq(80, 16)}\n")
    recs, problems = bdb.read_fasta_repair(str(f))
    assert len(problems) == 1
    assert [r.header for r in recs] == ["ISA-AB000001", "ISB-AB000003"]


def test_clean_sequence():
    seq, notes = bdb.clean_sequence("acgu nNr\t+-ACGT")
    assert seq == "ACGTNNRACGT"
    assert any("removed" in n for n in notes) and any("ambiguity" in n for n in notes)


#------------------------------------------------------- build / dedupe -----

def test_dedupe_identical_and_reverse_complement(tmp_path, inputs):
    kept, report, out = build(tmp_path, *inputs[:2], [inputs[2]])
    recs = read_fasta(out / "elements_all.fasta")
    # ISAb1 twice in IS.fna (identical) + its reverse complement in the TnCentral file -> one record
    assert report["dup_identical"] == 1 and report["dup_rc"] == 1
    hits = [k for k, s in recs.items() if s in (S1, bdb.revcomp(S1))]
    assert hits == ["IS200-IS605_ISAb1"]
    assert recs["IS200-IS605_ISAb1"] == S1           # the ISfinder record (orientation) is kept
    row = {r["id"]: r for r in read_tsv(out / "elements.tsv")}["IS200-IS605_ISAb1"]
    assert row["origin"] == "ISfinder"
    assert row["sources"] == "IS.fna;tncentral_test.fa"
    assert "reverse-complement" in row["notes"]


def test_family_mapping_and_types(tmp_path, inputs):
    kept, report, out = build(tmp_path, *inputs[:2], [inputs[2]])
    rows = {r["name"]: r for r in read_tsv(out / "elements.tsv")}
    expect = {
        "ISAb1": ("IS200-IS605", "IS"),
        "MITEAb1": ("MITE", "MITE"),
        "ISBj2_B": ("IS5", "IS"),
        "IS26.var1": ("IS6", "IS"),
        "IS10_p": ("IS4", "IS"),            # root IS10 -> IS10R in the taxonomy
        "Tn21": ("Tn", "transposon"),
        "TnAb1": ("Tn3", "transposon"),     # Tn distributed by ISfinder
        "In2": ("Integron", "integron"),
        "blaOXA-1": ("Cassette", "cassette"),
        "IS1326_IS1353": ("Compound", "compound"),
        "ISCR2": ("ISCR", "ISCR"),
    }
    for name, (fam, etype) in expect.items():
        assert (rows[name]["family"], rows[name]["element_type"]) == (fam, etype), name
    assert rows["ISBj2_B"]["group"] == "IS903"          # taxonomy fills an unknown group
    assert rows["MITEAb1"]["isfinder_family"] == "IS4"
    assert rows["IS26.var1"]["accession"] == "GU562437.2"
    assert rows["blaOXA-1"]["accession"] == "AF255921"
    # default types: IS, MITE, ISCR only
    sel = read_fasta(out / "elements.fasta")
    assert {rows[n]["id"] for n in ("Tn21", "TnAb1", "In2", "blaOXA-1", "IS1326_IS1353")}.isdisjoint(sel)
    assert {"IS200-IS605_ISAb1", "MITE_MITEAb1", "ISCR_ISCR2"} <= set(sel)
    # ISfinder-only subset
    iso = read_fasta(out / "elements_is_only.fasta")
    assert set(iso) == {"IS200-IS605_ISAb1", "MITE_MITEAb1", "IS5_ISBj2_B", "IS6_IS26"}
    assert all(sel[k] == v for k, v in iso.items())


def test_mite_parent_family(tmp_path, inputs):
    kept, _ = bdb.build(str(inputs[0]), str(inputs[1]), [], mite_family="parent")
    assert {r.name: r.family for r in kept}["MITEAb1"] == "IS4"


def test_unique_ids_and_header_convention(tmp_path, inputs):
    kept, report, out = build(tmp_path, *inputs[:2], [inputs[2]], types="all")
    lines = [l[1:].rstrip("\n") for l in open(out / "elements.fasta") if l.startswith(">")]
    assert len(lines) == len(set(lines)) == len(kept)
    for h in lines:
        assert h == h.split()[0]                        # no description / whitespace
        family, name = split_subject_id(h)
        assert "_" not in family and family and name
    # IS26-X00011 has the ISfinder name IS26 but another sequence: accession appended
    assert "IS6_IS26" in lines and "IS6_IS26-X00011" in lines
    seqs = read_fasta(out / "elements.fasta")
    assert all(set(s) <= set("ACGTRYKMSWBDHVN") for s in seqs.values())
    assert seqs["Cassette_blaOXA-1"] == rseq(500, 54)   # uppercased


def test_name_clash_without_accession_gets_counter(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(f">PICOTA_ISQ1\n{rseq(300, 60)}\n>ISQ1-ND\n{rseq(300, 61)}\n")
    kept, report = bdb.build(None, None, [str(f)])
    out = tmp_path / "o"
    bdb.write_outputs(kept, report, str(out), set(bdb.TYPES))
    ids = sorted(read_fasta(out / "elements.fasta"))
    assert ids == ["unknown_ISQ1", "unknown_ISQ1.2"]


def test_superseded_isfinder_copy_dropped(tmp_path, inputs):
    isf, tax, _ = inputs
    old = tmp_path / "picota_old.fasta"
    old.write_text(f">ISFinder_IS26\n{rseq(820, 70)}\n>ISFinder_ISNew1\n{rseq(900, 71)}\n")
    kept, report = bdb.build(str(isf), str(tax), [str(old)])
    names = [r.name for r in kept]
    assert names.count("IS26") == 1
    assert [r.name for r in report["superseded"]] == ["IS26"]
    new = next(r for r in kept if r.name == "ISNew1")
    assert new.family == "unknown" and new.origin == "ISfinder-copy"


def test_family_map_override(tmp_path, inputs):
    fm = tmp_path / "fm.tsv"
    fm.write_text("Tn21\tTn3\n")
    kept, _ = bdb.build(str(inputs[0]), str(inputs[1]), [str(inputs[2])], family_map=str(fm))
    assert {r.name: r.family for r in kept}["Tn21"] == "Tn3"


def test_parse_types():
    assert bdb.parse_types("default") == {"IS", "MITE", "ISCR"}
    assert bdb.parse_types("all") == set(bdb.TYPES)
    assert bdb.parse_types("is,Transposon") == {"IS", "transposon"}
    with pytest.raises(SystemExit):
        bdb.parse_types("IS,plasmid")


def test_cli_autodiscovery(tmp_path, inputs, capsys):
    out = tmp_path / "db"
    assert bdb.main(["--src", str(tmp_path), "--out", str(out), "--types", "IS,transposon"]) == 0
    for f in ("elements.fasta", "elements_all.fasta", "elements_is_only.fasta",
              "elements.tsv", "summary.txt"):
        assert (out / f).is_file(), f
    sel = read_fasta(out / "elements.fasta")
    assert "Tn_Tn21" in sel and "MITE_MITEAb1" not in sel
    summary = (out / "summary.txt").read_text()
    assert "reverse-complement identical" in summary and "Per family" in summary
