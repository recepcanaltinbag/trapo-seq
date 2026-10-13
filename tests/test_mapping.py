"""Tests for src/a_data_mapping.py + scripts/01_data_prep_from_py.sh, and an
end-to-end mapping + insertion-finding run on simulated data with ground truth."""
import csv
import os
import random
import shutil
import stat

import pytest

from conftest import REPO, TEST_DATA, requires_tools
from src import a_data_mapping
from src.a_data_mapping import MappingError, main_mapping
from src.e_insert_finder_from_bam import OUTPUT_NAME
from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch

MAP_OUTPUTS = ["02_mapped_to_plasmid.bam",
               "03_sorted_mapped_to_plasmid.bam", "03_sorted_mapped_to_plasmid.bam.bai",
               "04_mapped_to_plasmid.fastq", "07_sorted_aligned_to_other_genome.bam",
               "07_sorted_aligned_to_other_genome.bam.bai", "08_mapping_stats.txt"]
PLASMID = os.path.join(TEST_DATA, "plasmid.fasta")
GENOME = os.path.join(TEST_DATA, "genome.fasta")


def fastq_count(path):
    with open(path) as fh:
        return sum(1 for _ in fh) // 4


def fasta_records(path):
    with open(path) as fh:
        return [line[1:].strip() for line in fh if line.startswith(">")]


def stats(folder):
    out = {}
    with open(os.path.join(folder, "08_mapping_stats.txt")) as fh:
        for line in fh:
            k, v = line.rsplit(":", 1)
            out[k.strip()] = v.strip()
    return out


def leftover_temp_files(root):
    return [os.path.join(d, f) for d, _, files in os.walk(root) for f in files if f.startswith(".tmp_")]


# ------------------------------------------------------------ no tools needed

def test_script_path_independent_of_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = a_data_mapping.script_path()
    assert os.path.isabs(path)
    assert os.path.isfile(path)
    assert path == os.path.join(REPO, "scripts", "01_data_prep_from_py.sh")


def test_missing_inputs_raise(tmp_path):
    with pytest.raises(FileNotFoundError):
        main_mapping(str(tmp_path / "nope"), PLASMID, GENOME)
    with pytest.raises(FileNotFoundError):
        main_mapping(str(tmp_path), str(tmp_path / "missing.fasta"), GENOME)


def write_fake_tool(bindir, name, body):
    path = bindir / name
    path.write_text("#!/usr/bin/env bash\n" + body + "\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def test_failure_propagates_and_leaves_no_partial_output(tmp_path, monkeypatch):
    # B06: a failing minimap2 that has already written part of its output.
    bindir = tmp_path / "fake bin"
    bindir.mkdir()
    write_fake_tool(bindir, "minimap2",
                    'if [[ "$1" == --version ]]; then echo 0.fake; exit 0; fi\n'
                    'echo "@SQ partial output"\necho "fake minimap2 failure" >&2\nexit 3')
    write_fake_tool(bindir, "samtools", "exit 0")
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")

    data = tmp_path / "data dir"
    (data / "s 1").mkdir(parents=True)
    (data / "s 1" / "s 1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    with pytest.raises(MappingError) as exc:
        main_mapping(str(data), PLASMID, GENOME)
    assert exc.value.returncode != 0
    assert "fake minimap2 failure" in str(exc.value)
    assert sorted(os.listdir(data / "s 1")) == ["s 1.fastq"]


# ------------------------------------------------------------- with tools

def random_reads(path, n, length, seed=0):
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for i in range(n):
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * length}\n")


@requires_tools
def test_no_plasmid_reads_no_division_by_zero(tmp_path):
    # B07: nothing maps -> rate NA, exit status 0.
    data = tmp_path / "data dir"
    (data / "s1").mkdir(parents=True)
    random_reads(data / "s1" / "s1.fastq", 5, 300)
    main_mapping(str(data), PLASMID, GENOME)
    st = stats(data / "s1")
    assert st["Total reads"] == "5"
    assert st["Plasmid"] == "0"
    assert st["Transposition Rate"] == "NA"
    assert os.path.getsize(data / "s1" / "s1.fasta") == 0
    assert leftover_temp_files(data) == []


def plasmid_reads(path, n, seed=0):
    """Error-free plasmid fragments, half reverse-complemented."""
    with open(PLASMID) as fh:
        seq = "".join(line.strip() for line in fh if not line.startswith(">")).upper()
    comp = str.maketrans("ACGT", "TGCA")
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for i in range(n):
            start = rng.randrange(0, len(seq) - 3000)
            s = seq[start:start + 3000]
            if i % 2:
                s = s.translate(comp)[::-1]
            fh.write(f"@p{i} comment\n{s}\n+\n{'5' * len(s)}\n")


@requires_tools
def test_fasta_is_written_and_stale_outputs_rebuilt(tmp_path):
    data = tmp_path / "data dir"
    folder = data / "sample A"
    folder.mkdir(parents=True)
    plasmid_reads(folder / "sample A.fastq", 6)
    main_mapping(str(data), PLASMID, GENOME)

    for name in MAP_OUTPUTS + ["sample A.fasta"]:
        assert (folder / name).is_file(), name
    # B01: FASTA (portable conversion) holds the same reads as the mapped FASTQ.
    assert fasta_records(folder / "sample A.fasta") == [f"p{i}" for i in range(6)]
    assert fastq_count(folder / "04_mapped_to_plasmid.fastq") == 6
    assert stats(folder)["Plasmid"] == "6"
    assert leftover_temp_files(data) == []

    # Re-run without --force reuses everything ...
    mtimes = {n: os.stat(folder / n).st_mtime_ns for n in MAP_OUTPUTS}
    main_mapping(str(data), PLASMID, GENOME)
    assert all(os.stat(folder / n).st_mtime_ns == mtimes[n] for n in MAP_OUTPUTS[:-1])
    # ... but a removed upstream file triggers a rebuild of everything after it.
    os.remove(folder / "03_sorted_mapped_to_plasmid.bam")
    main_mapping(str(data), PLASMID, GENOME)
    assert os.stat(folder / "02_mapped_to_plasmid.bam").st_mtime_ns == mtimes["02_mapped_to_plasmid.bam"]
    # no intermediate SAM files are written any more
    assert not list(folder.glob("*.sam"))
    assert os.stat(folder / "03_sorted_mapped_to_plasmid.bam").st_mtime_ns != mtimes["03_sorted_mapped_to_plasmid.bam"]
    assert os.stat(folder / "07_sorted_aligned_to_other_genome.bam").st_mtime_ns != \
        mtimes["07_sorted_aligned_to_other_genome.bam"]


# ----------------------------------------------------------- end to end

POSITION_TOLERANCE = 20


@requires_tools
def test_end_to_end_recall_on_simulated_data(sim_dir, tmp_path, capsys):
    data = tmp_path / "sim copy"
    shutil.copytree(sim_dir, data)
    barcodes = sorted(d for d in os.listdir(data) if os.path.isdir(data / d))
    for bc in barcodes:
        shutil.copy(data / bc / f"{bc}_raw.fastq", data / bc / f"{bc}.fastq")

    main_mapping(str(data), PLASMID, GENOME)
    main_insert_finder_from_bam_batch(str(data), 500)

    calls = {}
    for bc in barcodes:
        with open(data / bc / OUTPUT_NAME) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                assert len(row) == 8
                calls.setdefault(row[0], []).append(row)
    with open(data / "truth.tsv") as fh:
        truth = list(csv.DictReader(fh, delimiter="\t"))

    def placed(t):
        rp, tsd = int(t["ref_pos"]), int(t["tsd_len"] or 0)
        return any(min(abs(int(c[4]) - rp), abs(int(c[4]) - rp - tsd)) <= POSITION_TOLERANCE
                   for c in calls.get(t["read_id"], []))

    ins = [t for t in truth if t["category"] == "is_insertion"]
    hits = [t for t in ins if placed(t)]
    recall = len(hits) / len(ins)
    # Strand of each call matches the read strand (B04).
    for t in hits:
        for c in calls[t["read_id"]]:
            assert c[7] == ("True" if t["read_strand"] == "-" else "False")
    empty = [t for t in truth if t["category"] == "empty_plasmid"]
    false_calls = sum(1 for t in empty if t["read_id"] in calls)
    with capsys.disabled():
        print(f"\n[end-to-end] is_insertion recall {len(hits)}/{len(ins)} = {recall:.2f} "
              f"(ref_pos within {POSITION_TOLERANCE} bp); calls on empty plasmids: {false_calls}/{len(empty)}")
    assert recall >= 0.85
    assert false_calls == 0
    assert leftover_temp_files(data) == []


# ------------------------------------------------- circular plasmid (B09)

COMP = str.maketrans("ACGT", "TGCA")


def first_fasta_seq(path):
    with open(path) as fh:
        lines = fh.read().split(">")[1].splitlines()
    return "".join(lines[1:]).upper()


def near_origin_sites(plasmid, element):
    """Insertion sites 0-800 bp from the origin on both sides. Each site is
    nudged so that the element's end bases differ from the plasmid bases they
    replace, so an error-free alignment cannot run a base into the element and
    ref_pos is exactly defined."""
    n = len(plasmid)
    sites = []
    for want in (15, 60, 150, 300, 500, 800, n - 15, n - 60, n - 150, n - 300, n - 500, n - 800):
        pos = want
        while plasmid[pos] == element[0] or plasmid[pos - 1] == element[-1] \
                or plasmid[pos + 1] == element[1] or plasmid[pos - 2] == element[-2]:
            pos += 1
        sites.append(pos)
    return sites


@requires_tools
def test_insertions_next_to_origin_exclude_plasmid(tmp_path):
    """B09: inserts within 0-800 bp of the plasmid origin (both sides), on reads
    linearised far from the origin (so every read crosses it), on both strands.
    The reported insert must be the element only (<= 10 bp overhang) and
    ref_pos must be the true site. (With the plain linear mapping,
    main_mapping(..., circular=False), about half of these reads fail.)"""
    plasmid = first_fasta_seq(PLASMID)
    element = first_fasta_seq(os.path.join(TEST_DATA, "IS_curated.fasta"))
    n, e = len(plasmid), len(element)
    cut = 3000                                    # restriction site, far from the origin
    truth = {}
    data = tmp_path / "data dir"
    (data / "s1").mkdir(parents=True)
    with open(data / "s1" / "s1.fastq", "w") as fh:
        for pos in near_origin_sites(plasmid, element):
            construct = plasmid[:pos] + element + plasmid[pos:]
            c = cut if cut <= pos else cut + e    # the cut site moves when the element precedes it
            read = construct[c:] + construct[:c]
            start = (pos - c) % len(construct)    # element interval on the forward read
            for strand in "+-":
                rid = f"site{pos}{'f' if strand == '+' else 'r'}"
                seq = read if strand == "+" else read.translate(COMP)[::-1]
                lo, hi = (start, start + e) if strand == "+" else (len(read) - start - e, len(read) - start)
                truth[rid] = (pos, lo, hi)
                fh.write(f"@{rid}\n{seq}\n+\n{'5' * len(seq)}\n")

    main_mapping(str(data), PLASMID, GENOME)
    main_insert_finder_from_bam_batch(str(data), 500)
    calls = {}
    with open(data / "s1" / OUTPUT_NAME) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            calls.setdefault(row[0], []).append((int(row[1]), int(row[2]), int(row[4])))

    problems = []
    for rid, (pos, lo, hi) in sorted(truth.items()):
        got = calls.get(rid, [])
        if len(got) != 1:
            problems.append(f"{rid}: {len(got)} calls {got}")
            continue
        qs, qe, ref_pos = got[0]
        if abs(qs - lo) > 10 or abs(qe - hi) > 10 or ref_pos != pos:
            problems.append(f"{rid}: insert {qs}-{qe} (true {lo}-{hi}), ref_pos {ref_pos} (true {pos})")
    assert problems == []
