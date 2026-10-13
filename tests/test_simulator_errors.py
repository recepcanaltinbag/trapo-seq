"""Simulator --error_model: the ONT-like model behaves as configured and the
default uniform model is unchanged."""
import csv
import hashlib
import json
import os
import random

import pytest

from simulate_traposome import (ONT_PROFILES, ONT_TRUTH_FIELDS, OntErrorModel, homopolymer_runs,
                                main, read_fasta, revcomp)

from conftest import TEST_DATA

ARGS = ["-p", os.path.join(TEST_DATA, "plasmid.fasta"),
        "-g", os.path.join(TEST_DATA, "genome.fasta"),
        "-i", os.path.join(TEST_DATA, "IS_curated.fasta")]
INDELS = ("I", "D", "hI", "hD")


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fastq(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    return {lines[i][1:]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 1, 4)}


def _model(seed=1, profile="r9.4.1", **over):
    p = dict(ONT_PROFILES[profile], **over)
    return OntErrorModel(random.Random(seed), p["read_accuracy_mean"], p["read_accuracy_sd"],
                         p["ont_sub_frac"], p["ont_ins_frac"], p["ont_del_frac"],
                         p["homopolymer_factor"])


def _random_seq(seed, n):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(scope="module")
def plasmid():
    return read_fasta(os.path.join(TEST_DATA, "plasmid.fasta"))[0][1]


# --------------------------------------------------------------------------
# error rate and composition
# --------------------------------------------------------------------------
@pytest.mark.parametrize("profile", sorted(ONT_PROFILES))
def test_error_rate_and_composition_match_profile(plasmid, profile):
    prof = ONT_PROFILES[profile]
    model = _model(profile=profile)
    accs, targets, counts = [], [], {"n_sub": 0, "n_ins": 0, "n_del": 0}
    for _ in range(400):
        _, _, st = model.read(plasmid)
        accs.append(st["read_accuracy"])
        targets.append(st["read_accuracy_target"])
        for k in counts:
            counts[k] += st[k]
    mean_acc = sum(accs) / len(accs)
    # realised accuracy tracks the per-read target and the configured mean
    assert abs(mean_acc - prof["read_accuracy_mean"]) < 0.25 * (1 - prof["read_accuracy_mean"])
    assert max(abs(a - t) for a, t in zip(accs, targets)) < 0.03
    # per-read spread, with a tail of poor reads below the mean
    sd = (sum((a - mean_acc) ** 2 for a in accs) / len(accs)) ** 0.5
    assert 0.5 * prof["read_accuracy_sd"] < sd < 1.5 * prof["read_accuracy_sd"]
    assert mean_acc - min(accs) > 2 * sd
    total = sum(counts.values())
    for key, frac in (("n_sub", "ont_sub_frac"), ("n_ins", "ont_ins_frac"), ("n_del", "ont_del_frac")):
        assert abs(counts[key] / total - prof[frac]) < 0.03, key
    assert (counts["n_ins"] + counts["n_del"]) / total >= 0.5


def test_custom_composition_and_accuracy_are_honoured(plasmid):
    model = _model(read_accuracy_mean=0.9, read_accuracy_sd=0.0, ont_sub_frac=0.1,
                   ont_ins_frac=0.1, ont_del_frac=0.8, homopolymer_factor=2)
    n = {"n_sub": 0, "n_ins": 0, "n_del": 0}
    accs = []
    for _ in range(200):
        _, _, st = model.read(plasmid)
        accs.append(st["read_accuracy"])
        for k in n:
            n[k] += st[k]
    total = sum(n.values())
    assert abs(sum(accs) / len(accs) - 0.9) < 0.005
    assert abs(n["n_del"] / total - 0.8) < 0.03 and abs(n["n_sub"] / total - 0.1) < 0.03


def test_error_log_agrees_with_alignment():
    """The logged errors explain the read: edit distance to the source equals
    (up to a few cancelling neighbours) the number of logged errors."""
    from Bio import Align
    aligner = Align.PairwiseAligner(mode="global", match_score=0, mismatch_score=-1,
                                    open_gap_score=-1, extend_gap_score=-1)
    model = _model(seed=5, read_accuracy_mean=0.9, read_accuracy_sd=0.02)
    total_dist = total_err = 0
    for k in range(15):
        src = _random_seq(100 + k, 600)
        log = []
        read, qual, st = model.read(src, log=log)
        assert len(read) == len(qual) == len(src) + st["n_ins"] - st["n_del"]
        assert len(log) == st["n_errors"]
        dist = -aligner.score(src, read)
        assert 0.8 * st["n_errors"] <= dist <= st["n_errors"]
        total_dist += dist
        total_err += st["n_errors"]
    # realised error by alignment ~ logged error (adjacent edits can merge)
    assert total_dist >= 0.9 * total_err
    assert abs(total_dist / (15 * 600) - 0.1) < 0.02


# --------------------------------------------------------------------------
# homopolymers
# --------------------------------------------------------------------------
def test_homopolymer_runs_get_more_indels():
    rng = random.Random(3)
    parts = []
    for _ in range(400):
        parts.append(_random_seq(rng.random(), 40))
        parts.append(rng.choice("ACGT") * rng.randint(5, 8))
    src = "".join(parts)
    runs = homopolymer_runs(src, min_run=2)
    in_long, in_short = set(), set()          # positions in runs >= 5 / in runs >= 2
    for L, starts in runs.items():
        for s in starts:
            in_short.update(range(s, s + L))
            if L >= 5:
                in_long.update(range(s, s + L))
    non_hp = len(src) - len(in_short)
    for profile in sorted(ONT_PROFILES):
        model = _model(seed=8, profile=profile, read_accuracy_sd=0.0)
        hp_indel = bg_indel = 0
        for _ in range(60):
            log = []
            model.read(src, log=log)
            for pos, kind, _ in log:
                if kind in INDELS:
                    if pos in in_long:
                        hp_indel += 1
                    elif pos not in in_short:
                        bg_indel += 1
        rate_hp = hp_indel / (60 * len(in_long))
        rate_bg = bg_indel / (60 * non_hp)
        assert rate_hp > 3 * rate_bg, (profile, rate_hp, rate_bg)


def test_homopolymer_errors_are_length_changes():
    src = "ACGT" * 20 + "AAAAAAA" + "TGCA" * 20
    model = _model(seed=2, read_accuracy_mean=0.9, read_accuracy_sd=0.0, ont_sub_frac=1,
                   ont_ins_frac=1, ont_del_frac=1, homopolymer_factor=50)
    seen = set()
    for _ in range(200):
        log = []
        read, _, _ = model.read(src, log=log)
        for pos, kind, _ in log:
            if kind in ("hI", "hD"):
                assert pos == 80
                seen.add(kind)
                if len(log) == 1:
                    assert read in (src[:80] + "A" * 6 + src[87:], src[:80] + "A" * 8 + src[87:])
    assert seen == {"hI", "hD"}


def test_homopolymer_factor_zero_disables_hp_errors(plasmid):
    model = _model(homopolymer_factor=0)
    assert all(model.read(plasmid)[2]["n_hp_indel"] == 0 for _ in range(20))


# --------------------------------------------------------------------------
# qualities
# --------------------------------------------------------------------------
def test_qualities_lower_at_errors_and_in_poor_reads(plasmid):
    model = _model(seed=4)
    at_err, overall = [], []
    by_read = []
    for _ in range(100):
        log = []
        read, qual, st = model.read(plasmid, log=log)
        q = [ord(c) - 33 for c in qual]
        assert all(2 <= x <= 40 for x in q)
        overall.append(sum(q) / len(q))
        at_err.extend(q[rp] for _, _, rp in log if rp < len(q))
        by_read.append((st["read_accuracy"], sum(q) / len(q)))
    assert sum(at_err) / len(at_err) < sum(overall) / len(overall) - 3
    # quality follows read accuracy: best quarter of reads have higher qualities
    by_read.sort()
    k = len(by_read) // 4
    assert sum(q for _, q in by_read[-k:]) / k > sum(q for _, q in by_read[:k]) / k + 1
    # not a constant string
    assert len(set(qual)) > 5


def test_minus_strand_is_reverse_complement(plasmid):
    a = _model(seed=9).read(plasmid, 100, 3000, "+", accuracy=0.95)
    b = _model(seed=9).read(plasmid, 100, 3000, "-", accuracy=0.95)
    assert b[0] == revcomp(a[0]) and b[1] == a[1][::-1] and a[2] == b[2]


# --------------------------------------------------------------------------
# end to end (main)
# --------------------------------------------------------------------------
def test_ont_run_is_reproducible_and_records_errors(tmp_path):
    common = ARGS + ["--barcodes", "1", "--reads", "60", "--error_model", "ont", "--seed", "5"]
    main(common + ["-o", str(tmp_path / "a")])
    main(common + ["-o", str(tmp_path / "b")])
    main(ARGS + ["--barcodes", "1", "--reads", "60", "--error_model", "ont", "--seed", "6",
                 "-o", str(tmp_path / "c")])
    fq = "barcode01/barcode01_raw.fastq"
    a = (tmp_path / "a" / fq).read_bytes()
    assert a == (tmp_path / "b" / fq).read_bytes()
    assert a != (tmp_path / "c" / fq).read_bytes()
    assert (tmp_path / "a" / "truth.tsv").read_bytes() == (tmp_path / "b" / "truth.tsv").read_bytes()
    rows = _tsv(tmp_path / "a" / "truth.tsv")
    reads = _fastq(tmp_path / "a" / fq)
    for r in rows:
        assert all(r[k] != "" for k in ONT_TRUTH_FIELDS)
        seq, qual = reads[r["read_id"]]
        assert len(seq) == len(qual) == int(r["read_len"])
        assert int(r["n_errors"]) == int(r["n_sub"]) + int(r["n_ins"]) + int(r["n_del"])
    meta = json.loads((tmp_path / "a" / "simulation.json").read_text())
    assert meta["error_model"] == "ont" and meta["ont_profile"] == "r9.4.1"
    assert meta["read_accuracy_mean"] == 0.92 and "approximation" in meta["ont_model_note"]
    assert "sub" not in meta and "mean_q" not in meta


def test_ont_overrides_and_profile(tmp_path):
    main(ARGS + ["--barcodes", "1", "--reads", "5", "--error_model", "ont",
                 "--ont_profile", "r10.4.1", "--read_accuracy_mean", "0.95",
                 "--ont_del_frac", "0.5", "-o", str(tmp_path / "o")])
    meta = json.loads((tmp_path / "o" / "simulation.json").read_text())
    assert meta["read_accuracy_mean"] == 0.95 and meta["ont_del_frac"] == 0.5
    assert meta["read_accuracy_sd"] == ONT_PROFILES["r10.4.1"]["read_accuracy_sd"]


def test_ont_structure_with_perfect_accuracy(tmp_path, plasmid):
    out = tmp_path / "perfect"
    main(ARGS + ["-o", str(out), "--barcodes", "1", "--reads", "150", "--error_model", "ont",
                 "--read_accuracy_mean", "1", "--read_accuracy_sd", "0", "--fragment_rate", "0.2",
                 "--novel_fraction", "0.3", "--seed", "3"])
    is_db = dict(read_fasta(os.path.join(TEST_DATA, "IS_curated.fasta")))
    reads = _fastq(out / "barcode01" / "barcode01_raw.fastq")
    n_ins = n_empty = 0
    for row in _tsv(out / "truth.tsv"):
        seq = reads[row["read_id"]][0]
        assert row["n_errors"] == "0" and float(row["read_accuracy"]) == 1.0
        if row["read_strand"] == "-":
            seq = revcomp(seq)
        if row["fragmented"] == "True":
            assert seq in plasmid + plasmid or row["category"] != "empty_plasmid"
            continue
        if row["category"] == "is_insertion":
            n_ins += 1
            p, t = int(row["ref_pos"]), int(row["tsd_len"])
            el = is_db[row["element"]]
            el = el if row["orientation"] == "+" else revcomp(el)
            tsd = plasmid[p:p + t]
            assert tsd == row["tsd_seq"] and tsd + el + tsd in seq + seq
        elif row["category"] == "empty_plasmid":
            n_empty += 1
            assert len(seq) == len(plasmid) and seq in plasmid + plasmid
    assert n_ins > 0 and n_empty > 0


def test_invalid_ont_parameters_exit(tmp_path):
    with pytest.raises(SystemExit):
        main(ARGS + ["-o", str(tmp_path / "x"), "--error_model", "ont",
                     "--read_accuracy_mean", "0.9", "--read_accuracy_sd", "0.5"])


# --------------------------------------------------------------------------
# uniform mode unchanged
# --------------------------------------------------------------------------
# sha256 of outputs of the simulator BEFORE --error_model was added.
UNIFORM_REFERENCE = {
    "--barcodes 1 --reads 80 --fragment_rate 0.4 --seed 7": {
        "barcode01/barcode01_raw.fastq": "92ad65f1688daa2f314eff24047b10a5936056d0ff1c2089462c3bc8801660f1",
        "truth.tsv": "2e2b4b5e7a1810ea0cb62fc638e03fb4130751c0edb15f64e9cb3f7b7b9c1b41",
        "clones.tsv": "d109106f63e14d85d1650ebc023220908e4cbe2b3ed3e0e141562731536b3fdd",
    },
    "--barcodes 1 --reads 60 --novel_fraction 0.5 --novel_mode element --sub 0.05 --ins 0.02 "
    "--del 0.03 --mean_q 12 --seed 9": {
        "barcode01/barcode01_raw.fastq": "b673396a306a705aa9b2856271a3d7a2c8145ea85c945338a97f496b34ec86fa",
        "truth.tsv": "56f0173ed51aaa01f1667f2355a38117578d0feaa095a0c2511d7f3b76fd332a",
        "clones.tsv": "3a5360c67a8c2492b6276cb047124816699cbd597d90aba7b59c1550735b7596",
    },
}


@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("cfg", sorted(UNIFORM_REFERENCE))
def test_uniform_mode_is_byte_identical(tmp_path, cfg, explicit):
    extra = ["--error_model", "uniform", "--ont_profile", "r10.4.1"] if explicit else []
    out = tmp_path / "u"
    main(ARGS + ["-o", str(out)] + cfg.split() + extra)
    for rel, digest in UNIFORM_REFERENCE[cfg].items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
    meta = json.loads((out / "simulation.json").read_text())
    assert not {"error_model", "ont_profile", "read_accuracy_mean"} & set(meta)
    assert "read_accuracy" not in _tsv(out / "truth.tsv")[0]


def test_existing_reference_hashes_with_explicit_uniform(tmp_path):
    from test_simulator_elements import REFERENCE
    out = tmp_path / "seg"
    main(ARGS + ["-o", str(out), "--barcodes", "2", "--reads", "60", "--novel_fraction", "0.4",
                 "--seed", "21", "--error_model", "uniform"])
    for rel, digest in REFERENCE.items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
