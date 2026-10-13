"""The simulator is the ground truth for every benchmark, so test it first."""
import csv
import os

from simulate_traposome import build_construct, linearise, main, read_fasta, revcomp

from conftest import TEST_DATA

PLASMID = "ACGTTGCAAGCTTGACCTAGGTACCGATCGATTACG"


def test_build_construct_duplicates_target_site():
    construct, tsd = build_construct(PLASMID, 10, "NNNNN", "+", 4)
    assert tsd == PLASMID[10:14]
    assert construct == PLASMID[:14] + "NNNNN" + PLASMID[10:14] + PLASMID[14:]
    assert len(construct) == len(PLASMID) + 5 + 4


def test_build_construct_minus_orientation():
    construct, _ = build_construct(PLASMID, 10, "AAACC", "-", 3)
    assert "GGTTT" in construct and "AAACC" not in construct


def test_linearise_shifts_cut_after_insertion():
    construct, _ = build_construct(PLASMID, 5, "GGGGG", "+", 2)
    mol = linearise(construct, 20, 5 + 2 - 1, 5 + 2)
    # the molecule starts at plasmid position 20 regardless of the insertion
    assert mol.startswith(PLASMID[20:26])
    assert sorted(mol) == sorted(construct)


def _reads(fastq):
    with open(fastq) as fh:
        lines = fh.read().split("\n")
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 4)}


def test_error_free_reads_match_truth(tmp_path):
    out = tmp_path / "sim"
    main(["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
          "-g", os.path.join(TEST_DATA, "genome.fasta"),
          "-i", os.path.join(TEST_DATA, "IS_curated.fasta"),
          "-o", str(out), "--barcodes", "1", "--reads", "120",
          "--sub", "0", "--ins", "0", "--del", "0", "--fragment_rate", "0", "--seed", "3"])
    plasmid = read_fasta(os.path.join(TEST_DATA, "plasmid.fasta"))[0][1]
    is_db = dict(read_fasta(os.path.join(TEST_DATA, "IS_curated.fasta")))
    reads = _reads(out / "barcode01" / "barcode01_raw.fastq")
    n_ins = 0
    with open(out / "truth.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq = reads[row["read_id"]]
            if row["read_strand"] == "-":
                seq = revcomp(seq)
            if row["category"] == "is_insertion":
                n_ins += 1
                p, t = int(row["ref_pos"]), int(row["tsd_len"])
                el = is_db[row["element"]]
                el = el if row["orientation"] == "+" else revcomp(el)
                tsd = plasmid[p:p + t]
                assert tsd == row["tsd_seq"]
                assert tsd + el + tsd in seq + seq
            elif row["category"] == "empty_plasmid":
                assert len(seq) == len(plasmid) and seq in plasmid + plasmid
    assert n_ins > 0
    assert (out / "clones.tsv").exists() and (out / "simulation.json").exists()


def test_seed_is_reproducible(tmp_path):
    args = ["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
            "-g", os.path.join(TEST_DATA, "genome.fasta"),
            "-i", os.path.join(TEST_DATA, "IS_curated.fasta"),
            "--barcodes", "1", "--reads", "30", "--seed", "5"]
    main(args + ["-o", str(tmp_path / "a")])
    main(args + ["-o", str(tmp_path / "b")])
    a = (tmp_path / "a" / "barcode01" / "barcode01_raw.fastq").read_text()
    b = (tmp_path / "b" / "barcode01" / "barcode01_raw.fastq").read_text()
    assert a == b
