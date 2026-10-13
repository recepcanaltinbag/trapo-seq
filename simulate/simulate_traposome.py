#!/usr/bin/env python3
"""
Simulate long-read traposome sequencing data with a known ground truth.

The model follows the wet-lab protocol trapo-seq was built for:

  1. A trap plasmid is carried by a host. Some plasmid molecules acquire an
     insertion: a curated IS element, or a "novel" segment of the host genome
     that is not in the curated IS database. Each insertion duplicates a short
     target site (TSD / direct repeat) of the plasmid.
  2. Insertion events form clones whose abundances are skewed (a few clones
     dominate after selection, many are rare).
  3. Plasmids are linearised at a single restriction site and sequenced as
     whole molecules; a fraction of molecules is fragmented during library
     prep, and reads come from either strand.
  4. Background: empty plasmids, plasmids carrying deletions, and genomic DNA
     contamination.

Every read is written with a record in truth.tsv so pipeline output can be
scored (see simulate/evaluate.py).

Only the Python standard library is used, so the simulator runs anywhere.

Coordinates: `ref_pos` is the 0-based plasmid position at which the element
is inserted, i.e. the element is placed between plasmid[ref_pos + tsd_len - 1]
and the duplicated copy of plasmid[ref_pos : ref_pos + tsd_len].

Novel inserts (--novel_mode)
----------------------------
  segment (default)  A novel clone carries a random host-genome segment
                     (--novel_len). This is the original behaviour; with the
                     default arguments the output is byte-identical to earlier
                     versions for the same --seed.
  element            --n_novel_elements synthetic mobile elements are built and
                     inserted into a copy of the host genome first; novel clones
                     then carry one of these elements. Each element has
                       * length drawn from --element_len (default 1200-2500 bp),
                       * terminal inverted repeats (TIRs) of 15-40 bp
                         (--tir_len); the right TIR is the reverse complement of
                         the left one with --tir_mismatch (default 0.1) of its
                         positions mutated (never the outer 2 bp),
                       * one ORF (ATG ... stop, no internal in-frame stop,
                         >= --min_orf bp, default 900) on a random strand, with
                         an in-frame stop codon right before the ATG so it is
                         the longest ORF of its frame,
                       * a genomic direct repeat (target-site duplication) of
                         2-9 bp (--element_dr) at its genomic copy.
                     In the plasmid, an element duplicates a target site of its
                     own DR length. Extra outputs in --out_dir:
                       genome_with_elements.fasta  modified host genome; the
                                                   pipeline must use it as -g
                       elements.fasta              element sequences
                       elements.tsv                element truth: element_id,
                         length, tir_len, tir_left, tir_right, tir_mismatches,
                         orf_start, orf_end, orf_strand (0-based half-open on
                         the element, stop codon included; strand relative to
                         the element), genome_contig, genome_insert_pos,
                         genome_end (element = modified_genome[pos:end]),
                         genome_dr_len, genome_dr_seq
                     truth.tsv / clones.tsv: element = element_id,
                     element_family = novel_element, genome_start/genome_end =
                     element coordinates in genome_with_elements.fasta.
                     Contamination reads are drawn from the modified genome.
                     All element randomness comes from a separate random stream
                     (seeded from --seed), so --seed still fixes the whole run.
"""

import argparse
import csv
import json
import os
import random
import sys

# Typical target-site duplication lengths per IS family. These are
# approximate family-level values (ISfinder family summaries); individual
# elements vary. Override with --tsd FAMILY=LEN.
DEFAULT_TSD = {
    "IS1": 9, "IS3": 3, "IS4": 11, "IS5": 4, "IS6": 8, "IS21": 5,
    "IS30": 3, "IS66": 8, "IS110": 0, "IS200": 0, "IS256": 9, "IS481": 7,
    "IS630": 2, "IS982": 8, "ISL3": 8, "Tn3": 5, "ISNCY": 4,
}
FALLBACK_TSD = 4
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    if not records:
        sys.exit(f"No sequences found in {path}")
    return records


def is_family(name):
    """IS names in the curated FASTA follow the >Family_Name convention."""
    return name.split("_", 1)[0] if "_" in name else name


def add_errors(seq, rng, sub, ins, dele):
    """Apply independent per-base substitution / insertion / deletion errors."""
    out = []
    bases = "ACGT"
    for b in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + sub:
            out.append(rng.choice([x for x in bases if x != b] or bases))
        else:
            out.append(b)
        if rng.random() < ins:
            out.append(rng.choice(bases))
    return "".join(out)


def quality_string(length, rng, mean_q):
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(mean_q, 4))))) for _ in range(length))


def skewed_weights(n, rng, alpha):
    """Pareto-like clone abundances: a few large clones, a long tail."""
    return [rng.paretovariate(alpha) for _ in range(n)]


def build_construct(plasmid, pos, element, orientation, tsd_len):
    """Return the plasmid with `element` inserted at `pos` plus its TSD."""
    insert = element if orientation == "+" else revcomp(element)
    tsd = plasmid[pos:pos + tsd_len]
    return plasmid[:pos + tsd_len] + insert + tsd + plasmid[pos + tsd_len:], tsd


def linearise(construct, cut, event_pos, event_len):
    """Rotate a circular construct so it starts at plasmid coordinate `cut`.

    `event_pos`/`event_len` describe the extra sequence (insert + TSD) added
    after plasmid coordinate event_pos; cut sites after it are shifted.
    """
    c = cut if event_pos is None or cut <= event_pos else cut + event_len
    return construct[c:] + construct[:c]


STOP_CODONS = ("TAA", "TAG", "TGA")
SENSE_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
                if a + b + c not in STOP_CODONS]


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_element(rng, length, tir_len, tir_mismatch, min_orf):
    """Build one synthetic mobile element.

    Returns (sequence, info) where info holds tir_left, tir_right,
    tir_mismatches, orf_start, orf_end, orf_strand (coordinates on the element,
    0-based half-open, stop codon included).
    """
    tir_left = random_seq(rng, tir_len)
    right = list(revcomp(tir_left))
    # Mutate positions of the right TIR, never its outer 2 bp (the element end).
    n_mm = int(round(tir_len * tir_mismatch))
    candidates = list(range(0, tir_len - 2))
    mm_pos = sorted(rng.sample(candidates, min(n_mm, len(candidates))))
    for i in mm_pos:
        right[i] = rng.choice([b for b in "ACGT" if b != right[i]])
    tir_right = "".join(right)

    interior_len = length - 2 * tir_len
    margin = 30
    max_orf = max(min_orf, interior_len - 2 * margin - 3)
    orf_len = rng.randint(min_orf // 3, max_orf // 3) * 3
    if orf_len + 3 + 2 * margin > interior_len:
        raise ValueError(f"element length {length} too short for an ORF of {orf_len} bp")
    codons = [rng.choice(SENSE_CODONS) for _ in range(orf_len // 3 - 2)]
    orf = "ATG" + "".join(codons) + rng.choice(STOP_CODONS)
    unit = "TAA" + orf            # in-frame stop right before the ATG
    free = interior_len - len(unit)
    offset = rng.randint(margin, free - margin)
    strand = rng.choice("+-")
    placed = unit if strand == "+" else revcomp(unit)
    interior = random_seq(rng, offset) + placed + random_seq(rng, free - offset)
    seq = tir_left + interior + tir_right
    if strand == "+":
        orf_start = tir_len + offset + 3
    else:
        orf_start = tir_len + offset
    info = {"tir_left": tir_left, "tir_right": tir_right, "tir_mismatches": len(mm_pos),
            "orf_start": orf_start, "orf_end": orf_start + len(orf), "orf_strand": strand}
    return seq, info


def insert_elements(genome, elements, rng, dr_range, margin=5000):
    """Insert each element once into a copy of `genome` ([(name, seq)]).

    Each insertion duplicates a target site of length drawn from `dr_range`.
    Returns (new_genome, placements) where placements[i] has genome_contig,
    genome_insert_pos, genome_end (in the modified genome), genome_dr_len and
    genome_dr_seq for elements[i].
    """
    lengths = [len(s) for _, s in genome]
    sites = []      # (contig index, position in the original contig, element index, dr_len)
    for idx in range(len(elements)):
        for _ in range(1000):
            ci = rng.choices(range(len(genome)), weights=lengths, k=1)[0]
            if lengths[ci] <= 2 * margin:
                continue
            pos = rng.randrange(margin, lengths[ci] - margin)
            if all(c != ci or abs(p - pos) > margin for c, p, _, _ in sites):
                break
        else:
            raise ValueError("genome too small to place the novel elements")
        sites.append((ci, pos, idx, rng.randint(*dr_range)))

    placements = [None] * len(elements)
    new_genome = []
    for ci, (name, seq) in enumerate(genome):
        parts, prev, shift = [], 0, 0
        for _, pos, idx, dr_len in sorted(s for s in sites if s[0] == ci):
            el = elements[idx]
            dr = seq[pos:pos + dr_len]
            parts.append(seq[prev:pos + dr_len])
            parts.append(el)
            start = pos + dr_len + shift
            placements[idx] = {"genome_contig": name, "genome_insert_pos": start,
                               "genome_end": start + len(el), "genome_dr_len": dr_len,
                               "genome_dr_seq": dr}
            shift += len(el) + dr_len
            prev = pos       # the target site is copied again after the element
        parts.append(seq[prev:])
        new_genome.append((name, "".join(parts)))
    return new_genome, placements


def write_fasta(path, records, width=70):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


ELEMENT_FIELDS = ["element_id", "length", "tir_len", "tir_left", "tir_right", "tir_mismatches",
                  "orf_start", "orf_end", "orf_strand", "genome_contig", "genome_insert_pos",
                  "genome_end", "genome_dr_len", "genome_dr_seq"]


def build_novel_elements(args, genome):
    """Element mode: make the elements, insert them into the genome, write
    genome_with_elements.fasta, elements.fasta and elements.tsv.

    Uses its own random stream so the main stream is untouched.
    Returns (modified genome, list of element dicts incl. 'seq').
    """
    erng = random.Random(f"novel-elements-{args.seed}")
    elements = []
    for k in range(args.n_novel_elements):
        length = erng.randint(*args.element_len)
        tir_len = erng.randint(*args.tir_len)
        seq, info = make_element(erng, length, tir_len, args.tir_mismatch, args.min_orf)
        elements.append({"element_id": f"novel_element_{k + 1}", "length": len(seq),
                         "tir_len": tir_len, **info, "seq": seq})
    new_genome, placements = insert_elements(genome, [e["seq"] for e in elements], erng,
                                             args.element_dr)
    for e, p in zip(elements, placements):
        e.update(p)
    write_fasta(os.path.join(args.out_dir, "genome_with_elements.fasta"), new_genome)
    write_fasta(os.path.join(args.out_dir, "elements.fasta"),
                [(e["element_id"], e["seq"]) for e in elements])
    with open(os.path.join(args.out_dir, "elements.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, ELEMENT_FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for e in elements:
            w.writerow(e)
    return new_genome, elements


def pick_position(rng, plen, cut, margin):
    """Choose an insertion site that does not overlap the restriction site."""
    while True:
        p = rng.randrange(0, plen)
        if abs(p - cut) > margin and abs(p - cut) < plen - margin:
            return p


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-p", "--plasmid", required=True, help="Trap plasmid FASTA (first record used)")
    ap.add_argument("-g", "--genome", required=True, help="Host genome FASTA (source of novel inserts and contamination)")
    ap.add_argument("-i", "--is_fasta", required=True, help="Curated IS FASTA (>Family_Name headers)")
    ap.add_argument("-o", "--out_dir", required=True, help="Output directory (barcode folders are created inside)")
    ap.add_argument("--barcodes", type=int, default=2, help="Number of barcodes/conditions (default 2)")
    ap.add_argument("--reads", type=int, default=2000, help="Reads per barcode (default 2000)")
    ap.add_argument("--insertion_rate", type=float, nargs="+", default=[0.3],
                    help="Fraction of plasmid molecules carrying an insertion; one value, or one per barcode")
    ap.add_argument("--clones", type=int, default=40, help="Distinct insertion events per barcode (default 40)")
    ap.add_argument("--novel_fraction", type=float, default=0.15, help="Fraction of clones whose insert is a genome segment absent from the IS DB")
    ap.add_argument("--novel_len", type=int, nargs=2, default=[1200, 3000], metavar=("MIN", "MAX"))
    ap.add_argument("--novel_mode", choices=["segment", "element"], default="segment",
                    help="Novel inserts: random genome segments (default, original behaviour) or "
                         "synthetic mobile elements with TIRs/ORF/DR inserted into the genome first")
    ap.add_argument("--n_novel_elements", type=int, default=2, help="element mode: number of elements (default 2)")
    ap.add_argument("--element_len", type=int, nargs=2, default=[1200, 2500], metavar=("MIN", "MAX"),
                    help="element mode: element length range (default 1200 2500)")
    ap.add_argument("--tir_len", type=int, nargs=2, default=[15, 40], metavar=("MIN", "MAX"),
                    help="element mode: TIR length range (default 15 40)")
    ap.add_argument("--tir_mismatch", type=float, default=0.1,
                    help="element mode: fraction of mismatched TIR positions (default 0.1)")
    ap.add_argument("--min_orf", type=int, default=900, help="element mode: minimum ORF length in bp (default 900)")
    ap.add_argument("--element_dr", type=int, nargs=2, default=[2, 9], metavar=("MIN", "MAX"),
                    help="element mode: genomic direct-repeat length range (default 2 9)")
    ap.add_argument("--deletion_rate", type=float, default=0.02, help="Fraction of plasmid molecules with a deletion")
    ap.add_argument("--contamination_rate", type=float, default=0.03, help="Fraction of reads from host genomic DNA")
    ap.add_argument("--fragment_rate", type=float, default=0.1, help="Fraction of molecules fragmented in library prep")
    ap.add_argument("--cut_pos", type=int, default=None, help="0-based restriction site on the plasmid (default: random)")
    ap.add_argument("--sub", type=float, default=0.02, help="Per-base substitution rate (default 0.02)")
    ap.add_argument("--ins", type=float, default=0.01, help="Per-base insertion error rate (default 0.01)")
    ap.add_argument("--del", dest="dele", type=float, default=0.01, help="Per-base deletion error rate (default 0.01)")
    ap.add_argument("--mean_q", type=float, default=18, help="Mean Phred quality written to FASTQ")
    ap.add_argument("--tsd", nargs="*", default=[], metavar="FAMILY=LEN", help="Override TSD lengths, e.g. IS4=10 Tn3=5")
    ap.add_argument("--seed", type=int, default=1, help="Random seed (default 1)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rng = random.Random(args.seed)

    plasmid_name, plasmid = read_fasta(args.plasmid)[0]
    genome = read_fasta(args.genome)
    is_db = read_fasta(args.is_fasta)
    tsd_table = dict(DEFAULT_TSD)
    for item in args.tsd:
        fam, _, length = item.partition("=")
        tsd_table[fam] = int(length)

    plen = len(plasmid)
    cut = args.cut_pos if args.cut_pos is not None else rng.randrange(plen)
    rates = args.insertion_rate
    if len(rates) not in (1, args.barcodes):
        sys.exit("--insertion_rate takes one value or one per barcode")

    os.makedirs(args.out_dir, exist_ok=True)
    elements = []
    if args.novel_mode == "element":
        genome, elements = build_novel_elements(args, genome)
    truth_fields = ["barcode", "read_id", "category", "clone_id", "element", "element_family",
                    "ref_pos", "orientation", "tsd_len", "tsd_seq", "insert_len",
                    "genome_contig", "genome_start", "genome_end", "read_strand", "fragmented", "read_len"]
    clone_fields = ["barcode", "clone_id", "category", "element", "element_family", "ref_pos",
                    "orientation", "tsd_len", "tsd_seq", "insert_len", "genome_contig",
                    "genome_start", "genome_end", "weight", "n_reads"]
    truth_fh = open(os.path.join(args.out_dir, "truth.tsv"), "w", newline="")
    clone_fh = open(os.path.join(args.out_dir, "clones.tsv"), "w", newline="")
    truth = csv.DictWriter(truth_fh, truth_fields, delimiter="\t")
    clones_out = csv.DictWriter(clone_fh, clone_fields, delimiter="\t")
    truth.writeheader()
    clones_out.writeheader()

    for b in range(1, args.barcodes + 1):
        bc = f"barcode{b:02d}"
        rate = rates[0] if len(rates) == 1 else rates[b - 1]
        bc_dir = os.path.join(args.out_dir, bc)
        os.makedirs(bc_dir, exist_ok=True)

        # Insertion clones for this barcode
        clones = []
        for k in range(args.clones):
            pos = pick_position(rng, plen, cut, margin=50)
            orient = rng.choice("+-")
            cid = f"{bc}_c{k:03d}"
            is_novel = rng.random() < args.novel_fraction
            if is_novel and elements:
                # element mode: one of the synthetic elements, TSD = its DR length
                el = rng.choice(elements)
                element, name, fam = el["seq"], el["element_id"], "novel_element"
                tsd_len = el["genome_dr_len"]
                category = "novel_insertion"
                gcoords = (el["genome_contig"], el["genome_insert_pos"], el["genome_end"])
            elif is_novel:
                contig, gseq = rng.choice(genome)
                ilen = rng.randint(*args.novel_len)
                gstart = rng.randrange(0, len(gseq) - ilen)
                element, name, fam = gseq[gstart:gstart + ilen], "novel", "novel"
                tsd_len = rng.randint(2, 9)
                category = "novel_insertion"
                gcoords = (contig, gstart, gstart + ilen)
            else:
                name, element = rng.choice(is_db)
                fam = is_family(name)
                tsd_len = tsd_table.get(fam, FALLBACK_TSD)
                category = "is_insertion"
                gcoords = ("", "", "")
            construct, tsd = build_construct(plasmid, pos, element, orient, tsd_len)
            clones.append({
                "barcode": bc, "clone_id": cid, "category": category, "element": name,
                "element_family": fam, "ref_pos": pos, "orientation": orient,
                "tsd_len": tsd_len, "tsd_seq": tsd, "insert_len": len(element),
                "genome_contig": gcoords[0], "genome_start": gcoords[1], "genome_end": gcoords[2],
                "weight": 0.0, "n_reads": 0, "_construct": construct,
            })
        weights = skewed_weights(len(clones), rng, alpha=1.2)
        total_w = sum(weights)
        for c, w in zip(clones, weights):
            c["weight"] = round(w / total_w, 6)

        with open(os.path.join(bc_dir, f"{bc}_raw.fastq"), "w") as fq:
            for r in range(args.reads):
                rid = f"{bc}_r{r:06d}"
                rec = {f: "" for f in truth_fields}
                rec.update(barcode=bc, read_id=rid, fragmented=False)
                x = rng.random()
                if x < args.contamination_rate:
                    contig, gseq = rng.choice(genome)
                    length = min(len(gseq), max(500, int(rng.lognormvariate(8.5, 0.6))))
                    start = rng.randrange(0, len(gseq) - length + 1)
                    mol = gseq[start:start + length]
                    rec.update(category="contamination", genome_contig=contig,
                               genome_start=start, genome_end=start + length)
                elif x < args.contamination_rate + args.deletion_rate:
                    dlen = rng.randint(50, 800)
                    dpos = pick_position(rng, plen, cut, margin=dlen + 50)
                    deleted = plasmid[:dpos] + plasmid[dpos + dlen:]
                    c = cut if cut <= dpos else cut - dlen
                    mol = deleted[c:] + deleted[:c]
                    rec.update(category="deletion", ref_pos=dpos, insert_len=-dlen)
                elif rng.random() < rate:
                    clone = rng.choices(clones, weights=weights, k=1)[0]
                    clone["n_reads"] += 1
                    extra = clone["insert_len"] + clone["tsd_len"]
                    mol = linearise(clone["_construct"], cut, clone["ref_pos"] + clone["tsd_len"] - 1, extra)
                    rec.update({k: clone[k] for k in ("category", "clone_id", "element_family", "ref_pos",
                                                       "orientation", "tsd_len", "tsd_seq", "insert_len",
                                                       "genome_contig", "genome_start", "genome_end")})
                    rec["element"] = clone["element"]
                else:
                    mol = linearise(plasmid, cut, None, 0)
                    rec.update(category="empty_plasmid")

                if rng.random() < args.fragment_rate and len(mol) > 1000:
                    a = rng.randrange(0, len(mol) - 500)
                    b_ = rng.randrange(a + 500, len(mol) + 1)
                    mol = mol[a:b_]
                    rec["fragmented"] = True
                strand = rng.choice("+-")
                if strand == "-":
                    mol = revcomp(mol)
                seq = add_errors(mol, rng, args.sub, args.ins, args.dele)
                rec.update(read_strand=strand, read_len=len(seq))
                fq.write(f"@{rid}\n{seq}\n+\n{quality_string(len(seq), rng, args.mean_q)}\n")
                truth.writerow(rec)

        for c in clones:
            clones_out.writerow({k: v for k, v in c.items() if not k.startswith("_")})

    truth_fh.close()
    clone_fh.close()
    with open(os.path.join(args.out_dir, "simulation.json"), "w") as fh:
        params = dict(vars(args))
        if args.novel_mode == "segment":
            # keep simulation.json identical to earlier versions in the default mode
            for k in ("novel_mode", "n_novel_elements", "element_len", "tir_len", "tir_mismatch",
                      "min_orf", "element_dr"):
                params.pop(k)
        else:
            params["elements"] = [{k: e[k] for k in ELEMENT_FIELDS} for e in elements]
        json.dump({**params, "cut_pos_used": cut, "plasmid_name": plasmid_name,
                   "plasmid_len": plen, "tsd_table": tsd_table}, fh, indent=2)
    print(f"Simulated {args.barcodes} barcode(s) x {args.reads} reads into {args.out_dir} (cut site {cut})")


if __name__ == "__main__":
    main()
