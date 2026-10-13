"""Read-based characterisation of Unannotated inserts (src/l_novel_elements.py)."""
import csv
import os
import random
import shutil

import pytest

from src.l_novel_elements import (NA, TSV_COLUMNS, characterise_consensus, cluster_by_similarity,
                                  main_novel_elements, majority_consensus, medoid,
                                  plasmid_dr_summary, revcomp)

from conftest import requires_tools

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not in PATH")

ANNOT_HEADER = ("Query ID\tSubject ID\tIdentity (%)\tScore\tE-value\tQuery Start\tQuery End\t"
                "Subject Start\tSubject End\tNote\tExplained\tref_pos\tis_Reverse\n")


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def noisy(rng, seq, rate):
    """Uniform substitutions / insertions / deletions at `rate` each third."""
    out = []
    for b in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(b + rng.choice("ACGT"))
        else:
            out.append(b)
    return "".join(out)


def edits(a, b):
    """Edit operations of the best global alignment of a with b or revcomp(b)."""
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode, al.match_score, al.mismatch_score = "global", 1, -1
    al.open_gap_score, al.extend_gap_score = -2, -1
    best = None
    for q in (b, revcomp(b)):
        c = al.align(a, q)[0].counts()
        e = c.mismatches + c.gaps
        best = e if best is None else min(best, e)
    return best


# ------------------------------------------------------------- pure Python

def test_majority_consensus():
    assert majority_consensus(["ACGT"]) == "ACGT"
    # column 3: one read has an insertion (gap majority) -> dropped
    assert majority_consensus(["ACG-T", "ACG-T", "ACGAT"]) == "ACGT"
    # one deletion out of three -> base kept
    assert majority_consensus(["ACGT", "A-GT", "ACGT"]) == "ACGT"
    # terminal gaps of a fragmented read are no vote, but a column needs >= half of the reads
    assert majority_consensus(["ACGTAA", "--GTAA", "--GTAA", "ACGTAA"]) == "ACGTAA"
    assert majority_consensus(["ACGTAA", "--GTAA", "--GTAA", "--GTAA"]) == "GTAA"
    # tie base vs gap -> base
    assert majority_consensus(["ACGT", "A-GT"]) == "ACGT"


def _pair(i, j, ci, cj, matches=100, strand="+"):
    return {(min(i, j), max(i, j)): {"strand": strand, "matches": matches, "alen": matches,
                                     "cov": {i: ci, j: cj}}}


def test_cluster_by_similarity_components_chimera_partial():
    # element A: 0,1,2 (1000 bp; 0-2 never aligned directly: single linkage)
    # element B: 3,4,5 (1200 bp); 6 = chimera A+B (2200 bp): joins B (its longer part,
    # coverage 0.545) but must not absorb A's cluster; 7 = 300 bp fragment of A
    lengths = [1000, 1000, 1000, 1200, 1200, 1200, 2200, 300]
    pairs = {}
    for p in (_pair(0, 1, .95, .95), _pair(1, 2, .9, .9), _pair(3, 4, .97, .97), _pair(4, 5, .97, .97),
              _pair(3, 5, .97, .97), _pair(6, 0, 1000 / 2200, .99), _pair(6, 3, 1200 / 2200, .99),
              _pair(7, 1, .99, .3)):
        pairs.update(p)
    cl = cluster_by_similarity(lengths, pairs)
    groups = sorted((sorted(c["members"]), sorted(c["partial"])) for c in cl)
    assert groups == [([0, 1, 2], [7]), ([3, 4, 5, 6], [])]
    a = next(c for c in cl if 0 in c["members"])
    assert a["seed"] == 0
    assert medoid([0, 1, 2], pairs) == 1          # aligned to both others


def test_plasmid_dr_summary():
    p = [{"sample": "s", "read_id": f"r{i}", "qstart": 1, "qend": 9} for i in range(5)]
    drs = {("s", "r0", 1, 9): (4, "CCTA"), ("s", "r1", 1, 9): (4, "CCTA"), ("s", "r2", 1, 9): (3, "CTA"),
           ("s", "r3", 1, 9): (4, "")}
    d = plasmid_dr_summary(p, drs)
    assert d == {"plasmid_dr_len": 4, "plasmid_dr_seq": "CCTA", "plasmid_dr_support": 0.75,
                 "plasmid_dr_counts": "4:3;3:1"}
    assert plasmid_dr_summary(p, {})["plasmid_dr_len"] == NA


def test_characterise_consensus_trims_to_tir():
    rng = random.Random(3)
    tir = rand_seq(rng, 25)
    orf = "ATG" + "".join(rng.choice(["GCT", "GAA", "CTG", "AAA"]) for _ in range(250)) + "TAA"
    element = tir + rand_seq(rng, 60) + orf + rand_seq(rng, 60) + revcomp(tir)
    cons = "CCTA" + element + "G"                  # TSD copy + one plasmid base kept by the reads
    row, el = characterise_consensus(cons)
    assert el == element and (row["start"], row["end"]) == (4, 4 + len(element))
    assert row["boundary_source"] == "consensus+TIR" and row["tir_left_len"] >= 25
    assert row["longest_orf_aa"] >= 251 and row["dr_len"] == NA
    assert row["verdict"] == "mobile_element_candidate"          # TIR + ORF
    row, _ = characterise_consensus(rand_seq(rng, 900), min_orf_aa=10_000,
                                    plasmid_dr={"plasmid_dr_len": 4, "plasmid_dr_support": 0.9})
    assert row["boundary_source"] == "consensus" and row["verdict"] == "no_element_evidence"


# ------------------------------------------------------------ with tools

@requires_tools
@needs_mafft
def test_unannotated_clusters_end_to_end(tmp_path):
    """Two unknown elements absent from the genome, inserted at several plasmid
    sites with ~6 % read errors, plus a singleton: two read clusters whose
    consensus equals the element, TIR found, plasmid DR from the DR table."""
    rng = random.Random(5)
    genome = rand_seq(rng, 20000)
    plasmid = rand_seq(rng, 4000)
    tir = rand_seq(rng, 30)
    el1 = tir + rand_seq(rng, 1100) + revcomp(tir)
    el2 = rand_seq(rng, 1500)
    d = tmp_path / "my data"
    (d / "s 1").mkdir(parents=True)
    (d / "insertions").mkdir()
    gpath = tmp_path / "genome.fasta"
    gpath.write_text(f">chr\n{genome}\n")
    reads, rows, drs = [], [], []
    k = 0
    for el, n, sites in ((el1, 12, (800, 1500, 2600)), (el2, 6, (3000,)), (rand_seq(rng, 900), 1, (100,))):
        for i in range(n):
            pos = sites[i % len(sites)]
            ins = noisy(rng, el if i % 2 else revcomp(el), 0.06)
            left, right = noisy(rng, plasmid[:pos + 4], 0.06), noisy(rng, plasmid[pos:], 0.06)
            rid = f"read_{k}"
            k += 1
            reads.append(f">{rid}\n{left + ins + right}\n")
            qs, qe = len(left), len(left) + len(ins)
            rows.append(f"{rid}_{qs}_{qe}\tno blast hit\tN/A\tN/A\tN/A\t0\t0\t0\t0\tUnannotated\t0\t{pos}\tFalse\n")
            drs.append(f"{rid},Unannotated,{qs},{qe},{pos + 1},{pos + 4},{pos + 1},4,{plasmid[pos:pos + 4]},+,DR\n")
    (d / "s 1" / "s 1.fasta").write_text("".join(reads))
    (d / "s 1" / "best_alignment.tab").write_text(ANNOT_HEADER + "".join(rows))
    (d / "insertions" / "s 1_insertions.csv").write_text(
        "Read ID,Subject ID,Query Start CSV,Query End CSV,Start Subject,End Subject,Insertion Point,"
        "Repeat Length,Overlap Sequence,Read Strand,Status\n" + "".join(drs))
    res = main_novel_elements(str(d), str(gpath), threads=2)
    assert [r["source"] for r in res] == ["reads", "reads"]
    by_len = sorted(res, key=lambda r: r["n_reads"], reverse=True)
    r1, r2 = by_len
    assert r1["cluster_id"].startswith("unannotated_") and r1["n_reads"] == 12 and r2["n_reads"] == 6
    assert r1["n_plasmid_sites"] == 3 and r2["n_plasmid_sites"] == 1
    assert r1["tir_found"] and r1["boundary_source"] == "consensus+TIR" and abs(r1["length"] - len(el1)) <= 2
    assert r1["plasmid_dr_len"] == 4 and r1["plasmid_dr_support"] == 1.0
    assert r1["contig"] == NA and r1["dr_len"] == NA and r1["genome_copies"] == 0
    assert r1["verdict"] in ("mobile_element_candidate", "partial_evidence")
    assert r1["read_identity"] > 0.85
    seqs = {}
    with open(d / "novel_elements" / "novel_elements.fasta") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    # 12 reads at 6 % error: the consensus is within a few edits of the element
    assert edits(seqs[r1["cluster_id"]], el1) <= 3
    c2 = seqs[r2["cluster_id"]]
    assert abs(len(c2) - len(el2)) <= 12 and edits(c2, el2) <= 14     # no TIR: TSD copy kept
    with open(d / "novel_elements" / "novel_elements.tsv") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == TSV_COLUMNS
        assert [r["source"] for r in rdr] == ["reads", "reads"]
    # unannotated=False: no read clusters
    assert main_novel_elements(str(d), str(gpath), unannotated=False) == []
