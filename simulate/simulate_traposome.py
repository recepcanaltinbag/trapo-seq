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
            if rng.random() < args.novel_fraction:
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
        json.dump({**vars(args), "cut_pos_used": cut, "plasmid_name": plasmid_name,
                   "plasmid_len": plen, "tsd_table": tsd_table}, fh, indent=2)
    print(f"Simulated {args.barcodes} barcode(s) x {args.reads} reads into {args.out_dir} (cut site {cut})")


if __name__ == "__main__":
    main()
