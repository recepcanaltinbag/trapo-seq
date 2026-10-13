"""Artefact evidence and flags (src/artefacts.py, docs/ARTEFACTS.md), their use in
blast_annot / is_stat / dr_finder / novel_elements / evaluate, and the simulator's
artefact reads (--chimera_rate, --ligation_chimera_rate, --dimer_rate,
--foldback_rate, --overhang)."""
import csv
import hashlib
import os
import random
import shutil

import pytest

from conftest import TEST_DATA, requires_tools
from src import artefacts as art
from src.artefacts import Piece, evaluate_call, revcomp

L = 3000
RNG = random.Random(5)
PLASMID = "".join(RNG.choice("ACGT") for _ in range(L))


def rand(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def plus(q0, q1, r0, record=None, mapq=60):
    return Piece(q0, q1, r0, r0 + (q1 - q0), False, mapq, "p", record)


def minus(q0, q1, r0, mapq=60):
    return Piece(q0, q1, r0, r0 + (q1 - q0), True, mapq, "p")


def info(pieces, length, seq=None):
    return {"pieces": sorted(pieces, key=lambda p: p.q_start), "length": length, "seq": seq}


CUT = {"cut_site": 100, "cut_end": 104, "end_offset": 4, "site_len": 6}


def flags_of(row):
    return row["flags"].split(",")


# --------------------------------------------------------------- enzymes / cut site

def test_enzyme_table_overhangs():
    assert art.enzyme_info("ecori")["overhang"] == "AATT"
    assert art.enzyme_info("EcoRI")["offset"] == 4
    assert art.enzyme_info("PstI")["overhang"] == "-TGCA" and art.enzyme_info("PstI")["offset"] == -4
    assert art.enzyme_info("EcoRV")["offset"] == 0 and art.enzyme_info("EcoRV")["overhang"] == ""
    assert art.enzyme_info("NotI")["overhang"] == "GGCC"
    assert art.parse_overhang("AATT") == 4 and art.parse_overhang("-TGCA") == -4 and art.parse_overhang("") == 0
    with pytest.raises(ValueError):
        art.enzyme_info("NoSuchI")


def test_recognition_sites_are_circular():
    seq = "TTCAAAAAAAAAAAAAAAGAA"           # GAATTC across the origin starts at 18
    assert art.recognition_sites(seq, "GAATTC") == [18]
    assert art.recognition_sites("AGAATTCA", "GAATTC") == [1]


def test_infer_cut_site_from_read_ends():
    reads = {}
    for i in range(20):
        n = L + 4
        first = L - 100
        if i % 2:
            reads[f"r{i}"] = info([plus(0, first, 100), plus(first, n, 0)], n)       # ends at 104 (0..104)
        else:   # minus strand: read starts at the molecule end (104) and runs back to 100
            reads[f"r{i}"] = info([minus(0, 104, 0), minus(104, n, 100)], n)
    reads["frag"] = info([plus(0, 500, 1500)], 500)                                 # fragment: not eligible
    cut = art.infer_cut_site(reads, {"p": L})
    assert cut["cut_site"] == 100 and cut["cut_end"] == 104 and cut["end_offset"] == 4
    assert cut["eligible_reads"] == 20 and cut["start_support"] == 20
    # random read ends (transposase kit): no dominant site
    reads = {f"r{i}": info([plus(0, L - s, s), plus(L - s, L, 0)], L)
             for i, s in enumerate(range(10, 3000, 150))}
    assert art.infer_cut_site(reads, {"p": L})["cut_site"] is None


def test_resolve_cut_site_with_enzyme_and_guess():
    seq = PLASMID[:1000] + "GAATTC" + PLASMID[1006:]
    seq = seq.replace("GAATTC", "GAATTA")[:1000] + "GAATTC" + seq[1006:]
    inferred = {"plasmid_len": L, "cut_site": 1001, "cut_end": 1005, "end_offset": 4}
    res = art.resolve_cut_site(dict(inferred), seq)
    assert "EcoRI" in res["enzyme_guess"] and res["overhang"] == "AATT" and res["site_len"] == 6
    res = art.resolve_cut_site({"plasmid_len": L, "cut_site": None}, seq, enzyme="EcoRI")
    assert (res["cut_site"], res["cut_end"], res["source"]) == (1001, 1005, "enzyme")
    res = art.resolve_cut_site(dict(inferred), seq, cut_site=50, overhang="-TGCA")
    assert (res["cut_site"], res["cut_end"], res["end_offset"], res["source"]) == (50, 46, -4, "user")


# --------------------------------------------------------------- per-call flags

def test_genuine_insertion_has_no_flags():
    # + read: plasmid ..1000+5 | insert | plasmid 1000.. (5 bp TSD)
    row = evaluate_call(("r", 1005, 2005, "SC"), info([plus(0, 1005, 0), plus(2005, 4005, 1000)], 4005), CUT, L)
    assert row["flanks"] == "both" and row["flank_gap"] == 5 and row["strand_consistent"] is True
    assert row["plasmid_left_end"] == 1005 and row["plasmid_right_start"] == 1000
    assert row["flags"] == "-"
    # the same on the minus strand: read-right flank is plasmid-left
    row = evaluate_call(("r", 2000, 3000, "SC"), info([minus(0, 2000, 1000), minus(3000, 4005, 0)], 4005), CUT, L)
    assert row["flank_gap"] == 5 and row["read_strand"] == "-" and row["flags"] == "-"


def test_nonadjacent_and_strand_flip_and_low_mapq():
    row = evaluate_call(("r", 1000, 2000, "SC"), info([plus(0, 1000, 0), plus(2000, 3000, 1500)], 3000), CUT, L)
    assert row["flank_gap"] == -500 and "nonadjacent_flanks" in flags_of(row)
    row = evaluate_call(("r", 1000, 2000, "SC"), info([plus(0, 1000, 500), plus(2000, 3000, 1450)], 3000), CUT, L)
    assert row["flank_gap"] == 50 and "nonadjacent_flanks" in flags_of(row)
    row = evaluate_call(("r", 1000, 2000, "SC"),
                        info([plus(0, 1000, 0), minus(2000, 3000, 1500, mapq=3)], 3000), CUT, L)
    assert row["strand_consistent"] is False and {"strand_flip", "low_mapq"} <= set(flags_of(row))


def test_cut_site_junction_and_overhang_dr():
    # two full molecules ligated with foreign DNA in between: left flank ends at the
    # molecule end (104, filled-in overhang), right flank starts at the cut (100)
    pieces = [plus(0, 2900, 104 - 2900 + L), plus(2900, 3000, 4), plus(4000, 7000, 100)]
    row = evaluate_call(("r", 3000, 4000, "SC"), info(pieces, 7000), CUT, L)
    assert row["flank_gap"] == 4 and row["dist_to_cut"] == 0
    assert {"cut_site_junction", "overhang_dr"} <= set(flags_of(row))
    # one-sided: foreign DNA ligated to the molecule end, insert reaches the read end
    row = evaluate_call(("r", 2900, 4000, "SC"), info([plus(0, 2900, 104 - 2900 + L)], 4000), CUT, L)
    assert row["flanks"] == "left_only" and "one_sided_cut_site" in flags_of(row)
    assert art.split_flags(row["flags"])[0] == ["one_sided_cut_site"]
    # one-sided away from the cut site: caution only
    row = evaluate_call(("r", 1000, 2000, "SC"), info([plus(0, 1000, 1500)], 2000), CUT, L)
    assert row["flags"] == "one_sided" and art.split_flags(row["flags"]) == ([], ["one_sided"])


def test_cut_in_insert_b3():
    site = "GAATTC"
    read = rand(1000, 1) + site + rand(994, 2)
    read = site[1:] + read[5:]        # the read starts inside the recognition site (after G^)
    ins_end = 1000
    row = evaluate_call(("r", 0, ins_end, "SC"), info([plus(ins_end, len(read), 1500)], len(read), read),
                        CUT, L, recog=[site[1:]])
    assert {"one_sided", "cut_in_insert"} <= set(flags_of(row))


def test_internal_adapter_foldback_and_plasmid_self():
    scanner = art.AdapterScanner(art.ADAPTERS)
    left, right = PLASMID[0:800], PLASMID[800:1600]
    ins = rand(400, 3) + art.ADAPTERS["ONT_ligation_Y_top"].replace("GTT", "GT") + rand(400, 4)
    read = left + ins + right
    qs, qe = len(left), len(left) + len(ins)
    row = evaluate_call(("r", qs, qe, "SC"), info([plus(0, qs, 0), plus(qe, len(read), 800)], len(read), read),
                        None, L, scanner)
    assert "internal_adapter" in flags_of(row) and row["adapter_hits"].startswith("ONT_ligation_Y_top")
    # no adapter in a random insert
    read2 = left + rand(len(ins), 9) + right
    row = evaluate_call(("r", qs, qe, "SC"), info([plus(0, qs, 0), plus(qe, len(read2), 800)], len(read2), read2),
                        None, L, scanner)
    assert row["adapter_hits"] == "-" and row["flags"] == "-"
    # fold-back: X + junk + revcomp(X)
    x = PLASMID[0:1500]
    junk = rand(700, 5)
    read3 = x + junk + revcomp(x)
    qs, qe = len(x), len(x) + len(junk)
    row = evaluate_call(("r", qs, qe, "SC"),
                        info([plus(0, qs, 0), minus(qe, len(read3), 0)], len(read3), read3), None, L)
    assert {"foldback", "strand_flip"} <= set(flags_of(row)) and int(row["foldback_kmers"]) > 100
    # plasmid self: the insert is plasmid sequence
    pk = art.PlasmidKmers(PLASMID)
    read4 = PLASMID[0:800] + PLASMID[1500:2500] + PLASMID[800:1600]
    row = evaluate_call(("r", 800, 1800, "SC"),
                        info([plus(0, 800, 0), plus(1800, 2600, 800)], 2600, read4), None, L, pk=pk)
    assert "plasmid_self" in flags_of(row) and float(row["plasmid_cover"]) > 0.9


def test_in_call_overlap_recovers_tsd():
    tsd = PLASMID[1000:1005]
    element = "".join(random.Random(8).choice("ACGT") for _ in range(800))
    # the aligner put element + TSD copy inside a CIGAR insertion after plasmid 1005
    read = PLASMID[0:1005] + element + tsd + PLASMID[1005:2000]
    qs, qe = 1005, 1005 + len(element) + len(tsd)
    row = evaluate_call(("r", qs, qe, "IN"), info([plus(0, qs, 0), plus(qe, len(read), 1005)], len(read), read),
                        None, L, plasmid_seq=PLASMID)
    assert row["flank_gap"] >= 5 and row["flags"] == "-"
    assert art.in_call_overlap(element + tsd, PLASMID, 1005)[1] >= 5


def test_unaligned_flank():
    # a 30 bp "flank" between two large insertions (bridged chimera junction)
    pieces = [plus(0, 1000, 0), plus(1600, 1630, 1000), plus(2200, 3000, 1030)]
    row = evaluate_call(("r", 1000, 1600, "IN"), info(pieces, 3000), None, L)
    assert "unaligned_flank" in flags_of(row)
    # a short stub at a read end (B1) or a piece split at the origin is fine
    row = evaluate_call(("r", 30, 1000, "SC"), info([plus(0, 30, 70), plus(1000, 3000, 100)], 3000), None, L)
    assert "unaligned_flank" not in flags_of(row)
    pieces = [plus(0, 1000, L - 1000), plus(1000, 1040, 0), plus(1700, 3000, 40)]
    row = evaluate_call(("r", 1040, 1700, "SC"), info(pieces, 3000), None, L)
    assert "unaligned_flank" not in flags_of(row)
    # low alignment identity of the flank next to the insert
    read = PLASMID[0:1000] + rand(600, 11) + PLASMID[1000:2000]
    bad = read[:900] + rand(100, 12) + read[1000:]
    rec = (0, ((0, 1000), (4, 2600)), 2600)
    p_left = plus(0, 1000, 0, record=rec)
    p_right = plus(1600, 2600, 1000, record=(1000, ((4, 1600), (0, 1000)), 2600))
    row = evaluate_call(("r", 1000, 1600, "SC"), info([p_left, p_right], 2600, bad), None, L, plasmid_seq=PLASMID)
    assert float(row["left_flank_identity"]) < 0.65 and "unaligned_flank" in flags_of(row)
    row = evaluate_call(("r", 1000, 1600, "SC"), info([p_left, p_right], 2600, read), None, L, plasmid_seq=PLASMID)
    assert float(row["left_flank_identity"]) == 1.0 and row["flags"] == "-"


# --------------------------------------------------------------- sites

def _call(qid, read, pos, strand, start=100, end=104, label="IS5_X", artefact=False):
    return {"query_id": qid, "read_id": read, "label": label, "pos": pos, "strand": strand,
            "start": start, "end": end, "artefact": artefact}


def test_site_table_single_molecule_and_independence():
    calls = [_call("a", "r1", 500, "+"),                                    # one read
             _call("b", "r2", 900, "+"), _call("c", "r3", 905, "-"),      # duplex-like pair
             _call("d", "r4", 1500, "+"), _call("e", "r5", 1502, "+"),    # same strand, same ends
             _call("f", "r6", 2000, "+"), _call("g", "r7", 2001, "-"), _call("h", "r8", 2003, "+"),
             _call("i", "r9", 2500, "+", artefact=True)]
    rows, member = art.site_table(calls, L)
    by = {r["position"]: r for r in rows}
    assert by[500]["single_molecule_support"] is True
    assert by[905]["single_molecule_support"] is True and by[905]["n_plus"] == 1
    assert by[1502]["single_molecule_support"] is False
    assert by[2001]["single_molecule_support"] is False and by[2001]["min_molecules"] == 2
    assert by[2500]["n_reads"] == 0 and by[2500]["n_artefact_flagged"] == 1
    assert by[2500]["single_molecule_support"] is False
    assert member["b"] == member["c"]


def test_cross_barcode_bleed_through():
    per = {"bc1": [{"site_id": "s1", "label": "IS5_X", "position": 1000, "n_reads": 100},
                   {"site_id": "s2", "label": "IS5_X", "position": 2000, "n_reads": 10}],
           "bc2": [{"site_id": "s9", "label": "IS5_X", "position": 1003, "n_reads": 2},
                   {"site_id": "s8", "label": "IS5_X", "position": 2000, "n_reads": 8}]}
    rows = art.cross_barcode_sites(per, L, bleed_pct=5)
    by = {r["position"]: r for r in rows}
    hit = [r for r in rows if r["possible_bleed_through"] != "-"]
    assert len(hit) == 1 and hit[0]["possible_bleed_through"] == "bc2"
    assert hit[0]["reads_bc1"] == 100 and hit[0]["reads_bc2"] == 2
    assert by[2000]["possible_bleed_through"] == "-"


# --------------------------------------------------------------- downstream modules

def test_blast_annot_flags_unsupported_and_b1_exemption():
    from src.g_blast_annot import artefact_flags_for, is_artefact
    hsp = lambda s, a, b: {"subject": s, "sstart": a, "send": b}     # noqa: E731
    calls = [
        {"qname": "q0", "query_id": "r0_0_900", "read_id": "r0", "note": "Genome",
         "hsp": hsp("chr", 1000, 3000), "explained": 99, "ref_pos": "500"},
        {"qname": "q1", "query_id": "r1_0_900", "read_id": "r1", "note": "Genome",
         "hsp": hsp("chr", 2500, 2900), "explained": 99, "ref_pos": "800"},          # supported by q0
        {"qname": "q2", "query_id": "r2_0_900", "read_id": "r2", "note": "Genome",
         "hsp": hsp("chr", 90000, 91000), "explained": 99, "ref_pos": "900"},        # unsupported
        {"qname": "q3", "query_id": "r3_0_900", "read_id": "r3", "note": "IS_DB",
         "hsp": hsp("IS5_X", 1, 900), "explained": 98, "ref_pos": "100"},            # B1 element
        {"qname": "q4", "query_id": "r4_0_900", "read_id": "r4", "note": "Genome",
         "hsp": hsp("chr", 5000, 6000), "explained": 99, "ref_pos": "100"},          # ligation
    ]
    ev = {"r0_0_900": {"flags": "-", "flanks": "both", "junction": "500", "read_strand": "+"},
          "r1_0_900": {"flags": "one_sided", "flanks": "left_only", "junction": "800"},
          "r2_0_900": {"flags": "one_sided", "flanks": "left_only", "junction": "900"},
          "r3_0_900": {"flags": "one_sided_cut_site", "flanks": "left_only", "junction": "100"},
          "r4_0_900": {"flags": "one_sided_cut_site", "flanks": "left_only", "junction": "100"}}
    flags, sites = artefact_flags_for(calls, ev, L)
    assert "one_sided_unsupported" not in flags["q1"] and "one_sided_unsupported" in flags["q2"]
    assert flags["q3"][0] == "one_sided_cut_site_is" and not is_artefact(flags["q3"])
    assert is_artefact(flags["q4"]) and is_artefact(flags["q2"]) and not is_artefact(flags["q1"])
    assert "single_molecule_support" in flags["q0"]


def test_is_stat_counts_artefact_flags(tmp_path):
    import pandas as pd
    from src import h_is_stats
    from src.g_blast_annot import OUTPUT_COLUMNS
    rows = [["q1", "IS5_IS5", 99, 900, 0, 1, 900, 1, 900, "IS_DB", 99, 10, False, "-", "-"],
            ["q2", "IS5_IS5", 99, 900, 0, 1, 900, 1, 900, "IS_DB", 99, 10, False, "-",
             "nonadjacent_flanks,single_molecule_support"],
            ["q3", "chr", 99, 900, 0, 1, 900, 1, 900, "Contamination", 99, 10, False, "-",
             "excluded:Genome:strand_flip,internal_adapter"],
            ["q4", "IS5_IS5", 99, 900, 0, 1, 900, 1, 900, "IS_DB", 99, 10, False, "-", "one_sided"]]
    path = tmp_path / "best_alignment.tab"
    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(path, sep="\t", index=False)
    s = h_is_stats.summarize_condition(str(path))
    assert s["artefact_flagged"] == 2 and s["artefact_excluded"] == 1
    assert s["artefact_flags"] == {"nonadjacent_flanks": 1, "single_molecule_support": 1, "strand_flip": 1,
                                   "internal_adapter": 1, "one_sided": 1}
    text = h_is_stats.format_condition("bc", s)
    assert "#Artefact_flagged=2" in text and "#Artefact_excluded=1" in text and "#Artefact_strand_flip=1" in text
    assert s["is_count"] == 3


def test_dr_finder_and_novel_elements_handle_flags(tmp_path):
    from src.j_dr_finder import _dr_subject
    from src.l_novel_elements import read_inserts
    assert _dr_subject({"Note": "IS_DB", "Subject ID": "IS5_X", "Artefact": "strand_flip"}) == "Artefact"
    assert _dr_subject({"Note": "IS_DB", "Subject ID": "IS5_X", "Artefact": "one_sided"}) == "IS5_X"
    assert _dr_subject({"Note": "Unannotated", "Subject ID": "no blast hit", "Artefact": "-"}) == "Unannotated"
    assert _dr_subject({"Note": "IS_DB", "Subject ID": "IS5_X"}) == "IS5_X"     # older tables
    annot = tmp_path / "best_alignment.tab"
    annot.write_text("Query ID\tNote\tArtefact\nr1_0_10\tGenome\t-\nr2_0_10\tGenome\tcut_site_junction\n"
                     "r3_0_10\tGenome\tsingle_molecule_support\n")
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">r1\nACGTACGTACGT\n>r2\nACGTACGTACGT\n>r3\nACGTACGTACGT\n")
    got = read_inserts([("bc", str(annot), str(fasta))], ("Genome",))
    assert sorted(i["read_id"] for i in got) == ["r1", "r3"]


def test_evaluate_flag_list_in_sync():
    import evaluate as ev
    from src.g_blast_annot import EXTRA_ARTEFACT_FLAGS
    assert set(ev.ARTEFACT_FLAGS) == set(art.ARTEFACT_FLAGS) | set(EXTRA_ARTEFACT_FLAGS)


def test_evaluate_artefact_metrics():
    import evaluate as ev

    def r(cat, before, after, flagged, flags="", rows=1, unflagged=True):
        return {"category": cat, "n_finder_rows": rows, "finder_called_unflagged": unflagged,
                "art_called_before": before, "art_called_after": after, "art_call_flagged": flagged,
                "art_flags": flags, "finder_art_flags": flags}
    rows = [r("ligation_chimera", True, False, True, "cut_site_junction", unflagged=False),
            r("ligation_chimera", False, False, False, "one_sided_cut_site", unflagged=False),
            r("chimera_plasmid_genome", True, True, False, "one_sided"),
            r("is_insertion", True, True, False, "single_molecule_support"),
            r("is_insertion", True, True, True, "strand_flip,one_sided")]
    d = ev.artefact_metrics(rows, True, True)
    assert d["ligation_chimera"]["called_before"] == 1 and d["ligation_chimera"]["called_after"] == 0
    assert d["ligation_chimera"]["flag_sensitivity"] == 1.0
    assert d["ligation_chimera"]["finder_called_before"] == 2 and d["ligation_chimera"]["finder_called_after"] == 0
    assert d["false_calls_before_total"] == 2 and d["false_calls_after_total"] == 1
    assert d["flagged_true_insertions"] == 1 and d["lost_true_insertions"] == 0
    assert d["caution_flags_true_insertions"] == "one_sided:1;single_molecule_support:1"
    assert ev.parse_artefact_value("excluded:Genome:strand_flip,one_sided", "Contamination") == \
        ("Genome", ["strand_flip", "one_sided"])
    assert ev.parse_artefact_value("-", "IS_DB") == ("IS_DB", [])


# --------------------------------------------------------------- simulator

SIM = ["-p", os.path.join(TEST_DATA, "plasmid.fasta"), "-g", os.path.join(TEST_DATA, "genome.fasta"),
       "-i", os.path.join(TEST_DATA, "IS_curated.fasta"), "--barcodes", "1", "--reads", "100", "--seed", "3"]


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fastq(path):
    lines = open(path).read().split("\n")
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 4)}


@pytest.mark.parametrize("model", [[], ["--error_model", "ont"]])
def test_simulator_artefact_reads_are_extra_and_leave_normal_reads(tmp_path, model):
    from simulate_traposome import main
    base, extra = tmp_path / "base", tmp_path / "extra"
    main(SIM + model + ["-o", str(base)])
    main(SIM + model + ["-o", str(extra), "--chimera_rate", "0.2", "--ligation_chimera_rate", "0.1",
                        "--dimer_rate", "0.1", "--foldback_rate", "0.1"])
    r0, r1 = _fastq(base / "barcode01/barcode01_raw.fastq"), _fastq(extra / "barcode01/barcode01_raw.fastq")
    assert all(r1[k] == v for k, v in r0.items()) and len(r1) == 150
    assert (base / "clones.tsv").read_bytes() == (extra / "clones.tsv").read_bytes()
    truth = _tsv(extra / "truth.tsv")
    cats = {t["category"] for t in truth[100:]}
    assert {"ligation_chimera", "dimer", "foldback"} <= cats and any(c.startswith("chimera_") for c in cats)
    assert all(t["artefact_detail"] for t in truth[100:]) and not any(t["artefact_detail"] for t in truth[:100])
    assert all(len(r1[t["read_id"]]) == int(t["read_len"]) for t in truth)


def test_simulator_structures_error_free(tmp_path):
    from simulate_traposome import main, read_fasta
    out = tmp_path / "o"
    perfect = ["--error_model", "ont", "--read_accuracy_mean", "1", "--read_accuracy_sd", "0"]
    main(SIM + perfect + ["-o", str(out), "--cut_pos", "500", "--overhang", "4",
                          "--ligation_chimera_rate", "0.3", "--dimer_rate", "0.2", "--foldback_rate", "0.2"])
    plasmid = read_fasta(os.path.join(TEST_DATA, "plasmid.fasta"))[0][1]
    lin = plasmid[500:] + plasmid[:500] + plasmid[500:504]
    reads = _fastq(out / "barcode01/barcode01_raw.fastq")
    for t in _tsv(out / "truth.tsv"):
        seq = reads[t["read_id"]]
        if t["category"] == "empty_plasmid" and t["fragmented"] == "False":
            assert seq in (lin, revcomp(lin))
        elif t["category"] == "dimer" and "ligation:head_to_tail" in t["artefact_detail"]:
            assert seq in (lin + lin, revcomp(lin + lin))
        elif t["category"] == "ligation_chimera":
            assert lin in seq or revcomp(lin) in seq
        elif t["category"] == "foldback":
            b = int(t["artefact_detail"].split("|")[0].split("-")[1])
            assert len(seq) > b
    import json
    meta = json.loads((out / "simulation.json").read_text())
    assert meta["overhang_seq"] == plasmid[500:504] and meta["ligation_chimera_rate"] == 0.3


def test_simulator_artefact_validation(tmp_path):
    from simulate_traposome import main
    for bad in (["--chimera_rate", "1.5"], ["--overhang", "XYZ"], ["--chimera_adapter_frac", "-1"]):
        with pytest.raises(SystemExit):
            main(SIM + ["-o", str(tmp_path / "x")] + bad)


# --------------------------------------------------------------- end to end

@requires_tools
def test_artefacts_end_to_end(tmp_path):
    """Chimeras and ligation products are flagged, the cut site and overhang are
    inferred, true insertions stay unflagged; exclude relabels."""
    import pandas as pd
    from simulate_traposome import main as simulate
    from src.a_data_mapping import main_mapping
    from src.e_insert_finder_from_bam_batch import main_insert_finder_from_bam_batch
    from src.g_blast_annot_batch import CROSS_SITES_NAME, main_annot_batch
    import evaluate as ev

    sim = tmp_path / "sim"
    simulate(SIM[:-4] + ["--barcodes", "2", "--reads", "250", "--seed", "7", "-o", str(sim),
                         "--chimera_rate", "0.12", "--ligation_chimera_rate", "0.08", "--overhang", "4",
                         "--error_model", "ont", "--ont_profile", "r10.4.1"])
    data = tmp_path / "data dir"
    for bc in ("barcode01", "barcode02"):
        (data / bc).mkdir(parents=True)
        shutil.copy(sim / bc / f"{bc}_raw.fastq", data / bc / f"{bc}.fastq")
    shutil.copy(os.path.join(TEST_DATA, "plasmid.fasta"), data / "plasmid.fasta")
    genome = os.path.join(TEST_DATA, "genome.fasta")
    main_mapping(str(data), str(data / "plasmid.fasta"), genome)
    main_insert_finder_from_bam_batch(str(data), 500)
    import json
    used = json.loads((sim / "simulation.json").read_text())["cut_pos_used"]
    cut = dict(line.rstrip("\n").split("\t", 1) for line in open(data / "barcode01" / "insertion_from_bam_cut_site.tsv"))
    assert int(cut["cut_site"]) == used and int(cut["end_offset"]) == 4

    main_annot_batch(str(data), genome, is_fasta=os.path.join(TEST_DATA, "IS_curated.fasta"), artefact_action="exclude")
    assert (data / CROSS_SITES_NAME).is_file() and (data / "barcode01" / "best_alignment_sites.tsv").is_file()
    report, per_read, _ = ev.evaluate(str(sim), str(data))
    a = report["overall"]["artefacts"]
    assert a["ligation_chimera"]["called_before"] >= 1 and a["ligation_chimera"]["called_after"] == 0
    assert a["false_calls_before_total"] >= 10
    assert a["false_calls_after_total"] <= 0.1 * a["false_calls_before_total"]
    assert a["flagged_true_insertions"] <= 0.01 * a["true_insertion_reads_called"]
    df = pd.read_csv(data / "barcode01" / "best_alignment.tab", sep="\t")
    excluded = df[df["Artefact"].astype(str).str.startswith("excluded:")]
    assert len(excluded) and set(excluded["Note"]) == {"Contamination"}


def test_default_simulation_unchanged_by_new_code(tmp_path):
    """The sha256 reference tests cover this; here: no artefact columns or keys."""
    from simulate_traposome import main
    import json
    out = tmp_path / "o"
    main(SIM + ["-o", str(out)])
    assert "artefact_detail" not in _tsv(out / "truth.tsv")[0]
    meta = json.loads((out / "simulation.json").read_text())
    assert not {"chimera_rate", "overhang", "dimer_rate", "foldback_rate"} & set(meta)
    assert hashlib.sha256((out / "truth.tsv").read_bytes()).hexdigest()
