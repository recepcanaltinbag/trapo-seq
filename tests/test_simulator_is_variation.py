"""Simulator --is_divergence / --is_truncate: per-clone variation of curated IS
elements, recorded in truth, and outputs unchanged when the options are off."""
import csv
import hashlib
import json
import os

import pytest

from simulate_traposome import (IS_DIV_FIELDS, IS_TRUNC_FIELDS, diverge_element, main, read_fasta,
                                revcomp, truncate_element)

from conftest import TEST_DATA

ARGS = ["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
        "-g", os.path.join(TEST_DATA, "genome.fasta"),
        "-i", os.path.join(TEST_DATA, "IS_curated.fasta")]
# perfect ONT reads: the ONT model draws errors from its own stream, so the read
# structure (main stream) does not depend on the inserted element's length
ERROR_FREE = ["--error_model", "ont", "--read_accuracy_mean", "1", "--read_accuracy_sd", "0",
              "--fragment_rate", "0"]
BASE_TRUTH = ["barcode", "read_id", "category", "clone_id", "element", "element_family", "ref_pos",
              "orientation", "tsd_len", "tsd_seq", "genome_contig", "genome_start", "genome_end",
              "read_strand", "fragmented"]


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _reads(fastq):
    with open(fastq) as fh:
        lines = fh.read().split("\n")
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 4)}


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# defaults: byte-identical to the simulator before these options existed
# --------------------------------------------------------------------------
# sha256 of outputs of the simulator as committed in c785d1b (before --is_divergence /
# --is_truncate), 1 barcode x 80 reads. The uniform model is covered by the
# reference tests in test_simulator_elements.py / test_simulator_errors.py.
ONT_REFERENCE = {
    "--error_model ont --fragment_rate 0.3 --seed 13": {
        "barcode01/barcode01_raw.fastq": "8c05841f81b81e68b20d2cded724e97941f3ab4a71ee1463a7016d551a151328",
        "truth.tsv": "abf263778376b18c08fc6677a7760f31d77167761337e6d32d1d53c09d1a2d0f",
        "clones.tsv": "5e90e6c613d99b94f418119a28573e9ea4ae913ddc82ccbc47f84483e253a260",
    },
    "--error_model ont --ont_profile r10.4.1 --novel_fraction 0.5 --novel_mode element --seed 14": {
        "barcode01/barcode01_raw.fastq": "baa660f61cfd349301857f5bafad4d748a55045442b109aa290644b1fcb579f7",
        "truth.tsv": "78828d99c35d84379931c0cd64432f9b503d44824f6b65eba3449cebebd2d2f7",
        "clones.tsv": "a0788a14e8b71336535d66df5d511b8ca3f399fc4a06b16cc1de85165d0a3777",
    },
    "--novel_fraction 0.5 --novel_mode element --seed 15": {
        "barcode01/barcode01_raw.fastq": "084fd49ffa12d33a7aa0e13e4bbb102aa25ea8b760c4b71190ce4966b5c8c0ad",
        "truth.tsv": "248652ce47298e5829f4a4da18a8e6dabc5de1a72de55386b32d04efaf662b28",
        "clones.tsv": "0e0539cd9204796ed56041a870817c94eacd7efe0c8d595a21710b4ca45a9ca0",
    },
}


@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("cfg", sorted(ONT_REFERENCE))
def test_defaults_are_byte_identical(tmp_path, cfg, explicit):
    extra = ["--is_divergence", "0", "--is_truncate_frac", "0.9"] if explicit else []
    out = tmp_path / "o"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "80"] + cfg.split() + extra)
    for rel, digest in ONT_REFERENCE[cfg].items():
        assert _sha(out / rel) == digest, rel
    meta = json.loads((out / "simulation.json").read_text())
    assert not {"is_divergence", "is_truncate", "is_truncate_frac"} & set(meta)
    assert not set(IS_TRUNC_FIELDS + IS_DIV_FIELDS) & set(_tsv(out / "truth.tsv")[0])


# --------------------------------------------------------------------------
# unit level
# --------------------------------------------------------------------------
def test_diverge_element_rate_and_composition():
    import random
    rng = random.Random(1)
    src = "".join(rng.choice("ACGT") for _ in range(20000))
    new, st = diverge_element(random.Random(2), src, 0.1)
    events = st["is_div_sub"] + st["is_div_ins"] + st["is_div_del"]
    assert abs(st["is_divergence"] - 0.1) < 0.01 and st["is_divergence"] == round(events / len(src), 5)
    assert abs(st["is_div_sub"] / events - 0.8) < 0.03
    assert abs(st["is_div_ins"] - st["is_div_del"]) < 0.25 * st["is_div_ins"]
    assert new != src and abs(len(new) - len(src)) < 0.02 * len(src)
    same, st0 = diverge_element(random.Random(2), src, 0.0)
    assert same == src and st0["is_divergence"] == 0


def test_truncate_element_bounds():
    import random
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(1500))
    for _ in range(200):
        piece, a, b = truncate_element(rng, seq, 300, 600)
        assert 300 <= len(piece) <= 600 and piece == seq[a:b] and 0 <= a < b <= len(seq)
    piece, a, b = truncate_element(rng, seq[:200], 300, 600)  # shorter than MIN: whole element
    assert (a, b) == (0, 200)


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------
def _inserts(out):
    """{read_id: forward-oriented inserted element} for error-free IS reads."""
    plasmid = read_fasta(os.path.join(TEST_DATA, "plasmid.fasta"))[0][1]
    plen = len(plasmid)
    pp = plasmid * 3
    reads = _reads(out / "barcode01" / "barcode01_raw.fastq")
    got = {}
    for r in _tsv(out / "truth.tsv"):
        if r["category"] != "is_insertion":
            continue
        seq = reads[r["read_id"]]
        if r["read_strand"] == "-":
            seq = revcomp(seq)
        p, t, n = int(r["ref_pos"]), int(r["tsd_len"]), int(r["insert_len"])
        flank = pp[plen + p + t - 60:plen + p + t]
        ss = seq + seq
        i = ss.find(flank)
        assert i >= 0
        ins = ss[i + 60:i + 60 + n]
        assert ss[i + 60 + n:i + 60 + n + t] == r["tsd_seq"]
        got[r["read_id"]] = (r, ins if r["orientation"] == "+" else revcomp(ins))
    return got


@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("isvar") / "default"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "200", "--seed", "8", "--clones", "60"] + ERROR_FREE)
    return out


def test_divergence_per_clone(tmp_path, default_run):
    out = tmp_path / "div"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "200", "--seed", "8",
                 "--is_divergence", "0.15", "--clones", "60"] + ERROR_FREE)
    is_db = dict(read_fasta(os.path.join(TEST_DATA, "IS_curated.fasta")))
    clones = {c["clone_id"]: c for c in _tsv(out / "clones.tsv")}
    for c in clones.values():
        assert not set(IS_TRUNC_FIELDS) & set(c)
        if c["category"] == "is_insertion":
            assert 0.08 < float(c["is_divergence"]) < 0.25
            n_ev = int(c["is_div_sub"]) + int(c["is_div_ins"]) + int(c["is_div_del"])
            assert float(c["is_divergence"]) == pytest.approx(n_ev / len(is_db[c["element"]]), abs=1e-5)
        else:
            assert all(c[f] == "" for f in IS_DIV_FIELDS)
    by_clone = {}
    for rid, (r, ins) in _inserts(out).items():
        assert len(ins) == int(r["insert_len"]) and ins != is_db[r["element"]]
        assert r["is_divergence"] == clones[r["clone_id"]]["is_divergence"]
        by_clone.setdefault(r["clone_id"], set()).add(ins)
    assert by_clone and all(len(v) == 1 for v in by_clone.values())  # one divergent copy per clone
    # separate random stream: read structure identical to the default run
    a, b = _tsv(default_run / "truth.tsv"), _tsv(out / "truth.tsv")
    assert [[r[k] for k in BASE_TRUTH] for r in a] == [[r[k] for k in BASE_TRUTH] for r in b]
    meta = json.loads((out / "simulation.json").read_text())
    assert meta["is_divergence"] == 0.15 and "is_truncate" not in meta


def test_truncation_records_kept_interval(tmp_path, default_run):
    out = tmp_path / "trunc"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "200", "--seed", "8", "--clones", "60",
                 "--is_truncate", "300", "600", "--is_truncate_frac", "0.5"] + ERROR_FREE)
    is_db = dict(read_fasta(os.path.join(TEST_DATA, "IS_curated.fasta")))
    clones = [c for c in _tsv(out / "clones.tsv") if c["category"] == "is_insertion"]
    n_trunc = sum(c["is_truncated"] == "True" for c in clones)
    assert 0.25 * len(clones) < n_trunc < 0.75 * len(clones)
    for c in clones:
        full = is_db[c["element"]]
        a, b = int(c["is_kept_start"]), int(c["is_kept_end"])
        assert int(c["is_full_len"]) == len(full) and int(c["insert_len"]) == b - a
        if c["is_truncated"] == "True":
            assert 300 <= b - a <= 600 or (a, b) == (0, len(full))
        else:
            assert (a, b) == (0, len(full))
        assert not set(IS_DIV_FIELDS) & set(c)
    n = 0
    for rid, (r, ins) in _inserts(out).items():
        assert ins == is_db[r["element"]][int(r["is_kept_start"]):int(r["is_kept_end"])]
        n += r["is_truncated"] == "True"
    assert n > 0
    a, b = _tsv(default_run / "truth.tsv"), _tsv(out / "truth.tsv")
    assert [[r[k] for k in BASE_TRUTH] for r in a] == [[r[k] for k in BASE_TRUTH] for r in b]


def test_truncation_then_divergence(tmp_path):
    out = tmp_path / "both"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "80", "--seed", "9", "--error_model", "ont",
                 "--is_truncate", "600", "1200", "--is_truncate_frac", "1", "--is_divergence", "0.05"])
    rows = _tsv(out / "truth.tsv")
    assert set(IS_TRUNC_FIELDS + IS_DIV_FIELDS) <= set(rows[0])
    is_rows = [r for r in rows if r["category"] == "is_insertion"]
    assert is_rows and all(r["is_truncated"] == "True" and r["is_divergence"] != "" for r in is_rows)
    assert all(r["is_truncated"] == "" for r in rows if r["category"] != "is_insertion")
    meta = json.loads((out / "simulation.json").read_text())
    assert meta["is_truncate"] == [600, 1200] and meta["is_truncate_frac"] == 1


@pytest.mark.parametrize("bad", [["--is_divergence", "1.2"], ["--is_divergence", "-0.1"],
                                 ["--is_truncate", "600", "300"], ["--is_truncate", "0", "10"],
                                 ["--is_truncate_frac", "1.5"]])
def test_invalid_options_exit(tmp_path, bad):
    with pytest.raises(SystemExit):
        main(ARGS + ["-o", str(tmp_path / "x")] + bad)
