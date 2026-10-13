"""Simulator --novel_mode element: synthetic mobile elements with TIRs, an ORF
and a genomic direct repeat, and the default mode left unchanged."""
import csv
import hashlib
import json
import os

import pytest

from simulate_traposome import STOP_CODONS, main, read_fasta, revcomp

from conftest import TEST_DATA

ARGS = ["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
        "-g", os.path.join(TEST_DATA, "genome.fasta"),
        "-i", os.path.join(TEST_DATA, "IS_curated.fasta")]


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _reads(fastq):
    with open(fastq) as fh:
        lines = fh.read().split("\n")
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 4)}


@pytest.fixture(scope="module")
def element_sim(tmp_path_factory):
    out = tmp_path_factory.mktemp("element sim")
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "150", "--novel_fraction", "0.5",
                 "--novel_mode", "element", "--n_novel_elements", "3",
                 "--sub", "0", "--ins", "0", "--del", "0", "--fragment_rate", "0", "--seed", "4"])
    return out


def test_element_structure_matches_truth(element_sim):
    elements = _tsv(element_sim / "elements.tsv")
    seqs = dict(read_fasta(str(element_sim / "elements.fasta")))
    assert len(elements) == 3
    for e in elements:
        seq = seqs[e["element_id"]]
        n = int(e["length"])
        t = int(e["tir_len"])
        assert len(seq) == n and 1200 <= n <= 2500 and 15 <= t <= 40
        # TIRs: right = revcomp(left) with the recorded number of mismatches,
        # outer 2 bp always intact
        assert seq.startswith(e["tir_left"]) and seq.endswith(e["tir_right"])
        assert len(e["tir_left"]) == len(e["tir_right"]) == t
        rc = revcomp(e["tir_right"])
        mm = sum(a != b for a, b in zip(e["tir_left"], rc))
        assert mm == int(e["tir_mismatches"]) == round(t * 0.1)
        assert e["tir_left"][:2] == rc[:2]
        # ORF: ATG ... stop, no internal in-frame stop, >= 900 bp
        a, b = int(e["orf_start"]), int(e["orf_end"])
        orf = seq[a:b] if e["orf_strand"] == "+" else revcomp(seq[a:b])
        assert len(orf) >= 900 and len(orf) % 3 == 0
        assert orf.startswith("ATG") and orf[-3:] in STOP_CODONS
        assert not any(orf[i:i + 3] in STOP_CODONS for i in range(0, len(orf) - 3, 3))
        assert t <= a and b <= n - t


def test_element_in_modified_genome_with_dr(element_sim):
    genome = dict(read_fasta(str(element_sim / "genome_with_elements.fasta")))
    original = dict(read_fasta(os.path.join(TEST_DATA, "genome.fasta")))
    seqs = dict(read_fasta(str(element_sim / "elements.fasta")))
    elements = _tsv(element_sim / "elements.tsv")
    total = sum(int(e["length"]) + int(e["genome_dr_len"]) for e in elements)
    for name, seq in original.items():
        assert len(genome[name]) == len(seq) + total
    for e in elements:
        g = genome[e["genome_contig"]]
        s, en, k = int(e["genome_insert_pos"]), int(e["genome_end"]), int(e["genome_dr_len"])
        assert g[s:en] == seqs[e["element_id"]]
        assert 2 <= k <= 9
        assert g[s - k:s] == g[en:en + k] == e["genome_dr_seq"]
        # removing the element and one DR copy restores the original genome locally
        assert g[s - 200:s] + g[en + k:en + 200] in original[e["genome_contig"]]


def test_clones_and_reads_carry_elements(element_sim):
    plasmid = read_fasta(os.path.join(TEST_DATA, "plasmid.fasta"))[0][1]
    seqs = dict(read_fasta(str(element_sim / "elements.fasta")))
    elements = {e["element_id"]: e for e in _tsv(element_sim / "elements.tsv")}
    reads = _reads(element_sim / "barcode01" / "barcode01_raw.fastq")
    novel = [r for r in _tsv(element_sim / "truth.tsv") if r["category"] == "novel_insertion"]
    assert novel
    for r in novel:
        e = elements[r["element"]]
        assert r["element_family"] == "novel_element"
        assert int(r["tsd_len"]) == int(e["genome_dr_len"])
        assert (r["genome_contig"], r["genome_start"], r["genome_end"]) == \
            (e["genome_contig"], e["genome_insert_pos"], e["genome_end"])
        p, t = int(r["ref_pos"]), int(r["tsd_len"])
        el = seqs[r["element"]] if r["orientation"] == "+" else revcomp(seqs[r["element"]])
        seq = reads[r["read_id"]]
        if r["read_strand"] == "-":
            seq = revcomp(seq)
        tsd = plasmid[p:p + t]
        assert tsd == r["tsd_seq"] and tsd + el + tsd in seq + seq
    meta = json.loads((element_sim / "simulation.json").read_text())
    assert meta["novel_mode"] == "element" and len(meta["elements"]) == 3


# sha256 of the outputs of the simulator BEFORE --novel_mode was added
# (same arguments, seed 21). The default mode must not change.
REFERENCE = {
    "barcode01/barcode01_raw.fastq": "1d33546a680a53e860c8b4af3aa943bfd885dabd11ea60f31d61ed8562c29d66",
    "barcode02/barcode02_raw.fastq": "6a0df2ee4d71e6abc54875383aa148ca9b33c54e560ad32f44f0cca0d6f4c62d",
    "truth.tsv": "a833e6aa25f53bb5809e0af42af1afc44e10001d41841231721ecc0333506bd6",
    "clones.tsv": "938565ac5c863ea0ba6e850526b9c0298ebac0f88935fd65f9159a24687931d2",
}


@pytest.mark.parametrize("explicit", [False, True])
def test_default_mode_is_byte_identical(tmp_path, explicit):
    extra = ["--novel_mode", "segment"] if explicit else []
    out = tmp_path / "seg"
    main(ARGS + ["-o", str(out), "--barcodes", "2", "--reads", "60", "--novel_fraction", "0.4",
                 "--seed", "21"] + extra)
    for rel, digest in REFERENCE.items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
    assert not (out / "elements.tsv").exists()
    assert not (out / "genome_with_elements.fasta").exists()
    assert "novel_mode" not in json.loads((out / "simulation.json").read_text())
