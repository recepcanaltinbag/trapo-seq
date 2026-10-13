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

Sequencing errors (--error_model)
---------------------------------
  uniform (default)  Independent per-base substitution / insertion / deletion
                     errors (--sub, --ins, --del) and Gaussian qualities around
                     --mean_q. Output is byte-identical to earlier versions for
                     the same --seed and options.
  ont                Oxford Nanopore-like reads:
                       * per-read accuracy ~ Beta(--read_accuracy_mean,
                         --read_accuracy_sd), clipped at 0.5 (left-skewed: a
                         tail of poor reads);
                       * error composition --ont_sub_frac / --ont_ins_frac /
                         --ont_del_frac (indel-dominant by default);
                       * homopolymer runs >= 3 bp get extra length-1 indels
                         (mostly deletions), probability growing with run
                         length (--homopolymer_factor);
                       * Phred qualities follow the read's accuracy, vary along
                         the read and dip at and next to errors.
                     --ont_profile r9.4.1 | r10.4.1 sets all of these. The
                     presets (r9.4.1: mean accuracy 0.92, sd 0.04; r10.4.1:
                     0.98, sd 0.015; see ONT_PROFILES) are ROUGH APPROXIMATIONS
                     chosen to look qualitatively like ONT data, not measured
                     values for any chemistry or basecaller; every value can be
                     overridden. Errors come from a separate random stream
                     (seeded from --seed); read structure still comes from the
                     main stream, but because the uniform model draws from the
                     main stream per base, the reads of an ont run differ from
                     those of a uniform run with the same seed. truth.tsv gains
                     read_accuracy_target, read_accuracy (1 - errors / source
                     length), n_errors, n_sub, n_ins, n_del, n_hp_indel.
                     The ont path is O(errors) per read in Python, so it is
                     much faster than the uniform per-base loop.
"""

import argparse
import bisect
import csv
import json
import math
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


# ---------------------------------------------------------------------------
# Oxford Nanopore-like error model (--error_model ont)
#
# All numbers below are rough approximations chosen to give reads that look
# qualitatively like ONT data (indel-dominant errors, homopolymer length
# errors, per-read accuracy spread with a tail of poor reads, error-correlated
# qualities). They are NOT measured values for any flow cell, chemistry or
# basecaller; override them with the CLI options when a specific profile is
# needed.
# ---------------------------------------------------------------------------
ONT_PROFILES = {
    #             read accuracy        error composition (fractions of all errors)
    "r9.4.1": {"read_accuracy_mean": 0.92, "read_accuracy_sd": 0.04,
               "ont_sub_frac": 0.40, "ont_ins_frac": 0.25, "ont_del_frac": 0.35,
               "homopolymer_factor": 4.0},
    "r10.4.1": {"read_accuracy_mean": 0.98, "read_accuracy_sd": 0.015,
                "ont_sub_frac": 0.40, "ont_ins_frac": 0.20, "ont_del_frac": 0.40,
                "homopolymer_factor": 5.0},
}
ONT_PARAM_KEYS = ("read_accuracy_mean", "read_accuracy_sd", "ont_sub_frac", "ont_ins_frac",
                  "ont_del_frac", "homopolymer_factor")
ONT_NOTE = ("ONT error-model parameters are rough approximations for simulation, "
            "not measured values for any flow cell, chemistry or basecaller.")
ONT_TRUTH_FIELDS = ["read_accuracy_target", "read_accuracy", "n_errors", "n_sub", "n_ins",
                    "n_del", "n_hp_indel"]
HP_MIN_RUN = 3          # homopolymer runs of at least this length get extra indels
HP_DEL_SHARE = 0.75     # share of homopolymer errors that are deletions (rest insertions)
HP_MAX_SHARE = 0.8      # homopolymer errors use at most this share of a read's error budget
MIN_READ_ACCURACY = 0.5
_I, _D, _HI, _HD = 4, 8, 12, 16       # event codes (substitution = 0-2, base in the low bits)
_KIND_NAMES = ("S", "I", "D", "hI", "hD")
_HP_RE = None


def homopolymer_runs(seq, min_run=HP_MIN_RUN):
    """{run length: sorted list of run start positions} for runs >= min_run."""
    global _HP_RE
    if _HP_RE is None:
        import re
        _HP_RE = re.compile(r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % ((min_run,) * 4))
    runs = {}
    for m in _HP_RE.finditer(seq):
        runs.setdefault(m.end() - m.start(), []).append(m.start())
    return runs


class OntErrorModel:
    """Fast ONT-like read generator.

    Work per read is O(errors + distinct homopolymer lengths) in Python; every
    per-base operation (slicing, joining, quality strings via str.translate)
    runs in C. Errors are placed by sampling the gaps between events
    (geometric skips) instead of a per-base loop.

    Model per read:
      * target accuracy a ~ Beta(mean, sd) clipped to [0.5, 1]; expected
        errors N = (1 - a) * source length. The Beta is left-skewed, giving a
        tail of poor reads.
      * homopolymer runs of length L >= 3 get a length-1 indel (deletion with
        probability 0.75, else insertion) with probability
        min(0.95, (1 - a) * homopolymer_factor * (L - 2)), drawn from the same
        error budget (scaled down so homopolymer deletions/insertions never
        exceed the configured deletion/insertion share);
      * the rest of the budget is spread uniformly as substitutions,
        insertions and deletions so that the read's overall composition
        matches ont_sub_frac / ont_ins_frac / ont_del_frac in expectation;
      * qualities: per-read level -10 log10(1 - a), varying in segments along
        the read, per-base noise, and low values at and next to errors.
    """

    QUAL_POOL = 1 << 17
    QUAL_SD = 3.0

    def __init__(self, rng, accuracy_mean, accuracy_sd, sub_frac, ins_frac, del_frac,
                 homopolymer_factor):
        if not 0 < accuracy_mean <= 1:
            raise ValueError("read accuracy mean must be in (0, 1]")
        if accuracy_sd < 0:
            raise ValueError("read accuracy sd must be >= 0")
        if accuracy_mean < 1 and accuracy_sd > 0 and accuracy_sd ** 2 >= accuracy_mean * (1 - accuracy_mean):
            raise ValueError("read accuracy sd too large for its mean (need sd^2 < mean*(1-mean))")
        total = sub_frac + ins_frac + del_frac
        if min(sub_frac, ins_frac, del_frac) < 0 or total <= 0:
            raise ValueError("error fractions must be >= 0 and not all 0")
        if homopolymer_factor < 0:
            raise ValueError("homopolymer factor must be >= 0")
        self.rng = rng
        self.mean, self.sd = accuracy_mean, accuracy_sd
        self.fs, self.fi, self.fd = sub_frac / total, ins_frac / total, del_frac / total
        self.hp = homopolymer_factor
        if accuracy_mean < 1 and accuracy_sd > 0:
            k = accuracy_mean * (1 - accuracy_mean) / accuracy_sd ** 2 - 1
            self.beta = (accuracy_mean * k, (1 - accuracy_mean) * k)
        else:
            self.beta = None
        self._runs = {}
        self._tables = {}
        # quality noise pool: chars encode an offset o in [-12, 12] as chr(o + 64)
        g = rng.gauss
        pool = "".join(chr(64 + max(-12, min(12, int(round(g(0, self.QUAL_SD))))))
                       for _ in range(self.QUAL_POOL))
        self._pool = pool + pool
        self._subs = {b: "".join(x for x in "ACGT" if x != b) for b in "ACGT"}

    def _table(self, q):
        t = self._tables.get(q)
        if t is None:
            t = {64 + o: chr(33 + max(2, min(40, q + o))) for o in range(-12, 13)}
            self._tables[q] = t
        return t

    def runs_for(self, mol):
        """Homopolymer runs of a molecule, cached (molecules repeat across reads)."""
        r = self._runs.get(mol)
        if r is None:
            if len(self._runs) > 512:
                self._runs.clear()
            r = self._runs[mol] = homopolymer_runs(mol)
        return r

    def draw_accuracy(self):
        if self.beta is None:
            return min(1.0, self.mean)
        return max(MIN_READ_ACCURACY, self.rng.betavariate(*self.beta))

    def _skips(self, p, n):
        """Sorted positions in [0, n) hit by independent Bernoulli(p) events
        (geometric gaps, so the cost is proportional to the number of hits)."""
        if p <= 0 or n <= 0:
            return []
        if p >= 1:
            return list(range(n))
        rnd, log, lp = self.rng.random, math.log, math.log1p(-p)
        out, i = [], -1
        while True:
            i += 1 + int(log(1.0 - rnd()) / lp)
            if i >= n:
                return out
            out.append(i)

    def _quality_base(self, length, qread):
        """Background qualities: a level per segment (~2 kb) around qread plus
        per-base noise sliced from a pre-drawn pool (str.translate, no per-base
        Python work)."""
        rng, rnd = self.rng, self.rng.random
        pool, plen = self._pool, self.QUAL_POOL
        parts, i = [], 0
        while i < length:
            L = min(length - i, 1 + int(rng.expovariate(1 / 2000.0)))
            off = int(rnd() * plen)
            noise = pool[off:off + L] if L <= plen else (pool * (L // plen + 1))[off:off + L]
            parts.append(noise.translate(self._table(int(round(qread + rng.gauss(0, 1.5))))))
            i += L
        return "".join(parts)

    def read(self, mol, start=0, end=None, strand="+", accuracy=None, log=None):
        """Sequence mol[start:end] with ONT-like errors.

        Returns (read, quality, stats). Errors are applied to the forward
        window; for strand '-' the read is reverse-complemented and the
        qualities reversed. If `log` is a list, (source_pos, kind, read_pos)
        tuples are appended: source_pos relative to the window, read_pos in the
        forward read (before any reverse complement) of the substituted base,
        the inserted base (insertions go after source_pos) or the base right
        after a deletion; kind is S, I, D, hI or hD (h = homopolymer).
        """
        rng = self.rng
        rnd = rng.random
        end = len(mol) if end is None else end
        seg = mol[start:end]
        n = len(seg)
        acc = self.draw_accuracy() if accuracy is None else accuracy
        e = 1.0 - acc
        budget = e * n
        # events are ints: (source position << 5) | code, code = kind << 2 | base
        events = []
        add = events.append
        hp_ins = hp_del = 0.0
        if e > 0 and self.hp > 0 and n:
            probs = []
            for L, starts in self.runs_for(mol).items():
                i0 = bisect.bisect_left(starts, start)
                i1 = bisect.bisect_right(starts, end - L)
                if i1 > i0:
                    probs.append((starts, i0, i1, min(0.95, e * self.hp * (L - 2))))
            h = sum(p * (i1 - i0) for _, i0, i1, p in probs)
            if h > 0:
                scale = min(1.0, HP_MAX_SHARE * budget / h,
                            self.fi * budget / (h * (1 - HP_DEL_SHARE)),
                            self.fd * budget / (h * HP_DEL_SHARE))
                for starts, i0, i1, p in probs:
                    for k in self._skips(p * scale, i1 - i0):
                        add(((starts[i0 + k] - start) << 5) | (_HD if rnd() < HP_DEL_SHARE else _HI))
                hp_ins = h * scale * (1 - HP_DEL_SHARE)
                hp_del = h * scale * HP_DEL_SHARE
        if e > 0 and n:
            bs = self.fs * budget
            bi = max(0.0, self.fi * budget - hp_ins)
            bd = max(0.0, self.fd * budget - hp_del)
            bg = bs + bi + bd
            p = bg / n
            if p > 0:
                ts, ti = bs / bg, (bs + bi) / bg
                tsi = ti - ts
                log_ = math.log
                lp = math.log1p(-p) if p < 1 else None
                i = -1
                while True:
                    i += 1 if lp is None else 1 + int(log_(1.0 - rnd()) / lp)
                    if i >= n:
                        break
                    # one draw picks the error type and, within its band, the base
                    u = rnd()
                    if u < ts:
                        add((i << 5) | int(u / ts * 3))
                    elif u < ti:
                        add((i << 5) | _I | int((u - ts) / tsi * 4))
                    else:
                        add((i << 5) | _D)
        events.sort()

        qread = 40.0 if e <= 1e-4 else -10.0 * math.log10(e)
        qbase = self._quality_base(n + len(events) + 3, qread)
        qn = chr(33 + max(2, int(qread) // 2))
        lows = [qn + chr(35 + k) + qn for k in range(6)]   # Q2-7 at the error, qread/2 beside it
        rp, qp = [], []
        prev = qprev = shift = 0
        ns = ni = nd = nh = 0
        subs = self._subs
        for ev in events:
            pos = ev >> 5
            if pos < prev:          # one event per source position
                continue
            code = ev & 31
            o = pos + shift
            if code < 4:                            # substitution
                rp.append(seg[prev:pos])
                rp.append(subs.get(seg[pos], "ACG")[code])
                c = o
                ns += 1
            elif code < 8 or code == _HI:           # insertion after seg[pos]
                rp.append(seg[prev:pos + 1])
                rp.append(seg[pos] if code == _HI else "ACGT"[code & 3])
                c = o + 1
                shift += 1
                ni += 1
                nh += code == _HI
            else:                                   # deletion of seg[pos]
                rp.append(seg[prev:pos])
                c = o
                shift -= 1
                nd += 1
                nh += code == _HD
            prev = pos + 1
            if c > qprev:
                qp.append(qbase[qprev:c - 1])
                qp.append(lows[pos % 6])
                qprev = c + 2
            elif c == qprev:
                qp.append(lows[pos % 6][1:])
                qprev = c + 2
            if log is not None:
                log.append((pos, _KIND_NAMES[code >> 2], c))
        rp.append(seg[prev:])
        read = "".join(rp)
        m = len(read)
        if qprev < m:
            qp.append(qbase[qprev:m])
        qual = "".join(qp)[:m]
        if strand == "-":
            read = revcomp(read)
            qual = qual[::-1]
        nerr = ns + ni + nd
        stats = {"read_accuracy_target": round(acc, 5),
                 "read_accuracy": round(1.0 - nerr / n, 5) if n else 1.0,
                 "n_errors": nerr, "n_sub": ns, "n_ins": ni, "n_del": nd, "n_hp_indel": nh}
        return read, qual, stats


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
    ap.add_argument("--error_model", choices=["uniform", "ont"], default="uniform",
                    help="uniform: independent per-base errors from --sub/--ins/--del and --mean_q "
                         "(default, original behaviour); ont: Oxford Nanopore-like model (per-read "
                         "accuracy, indel-dominant, homopolymer indels, error-correlated qualities). "
                         "--sub/--ins/--del/--mean_q are ignored in ont mode")
    ap.add_argument("--ont_profile", choices=sorted(ONT_PROFILES), default="r9.4.1",
                    help="ont mode: preset for the parameters below (default r9.4.1). Presets are "
                         "rough approximations, not measured values")
    ap.add_argument("--read_accuracy_mean", type=float, default=None,
                    help="ont mode: mean per-read accuracy (r9.4.1 0.92, r10.4.1 0.98)")
    ap.add_argument("--read_accuracy_sd", type=float, default=None,
                    help="ont mode: SD of per-read accuracy, Beta-distributed (r9.4.1 0.04, r10.4.1 0.015)")
    ap.add_argument("--ont_sub_frac", type=float, default=None,
                    help="ont mode: fraction of errors that are substitutions (r9.4.1 0.40, r10.4.1 0.40)")
    ap.add_argument("--ont_ins_frac", type=float, default=None,
                    help="ont mode: fraction of errors that are insertions (r9.4.1 0.25, r10.4.1 0.20)")
    ap.add_argument("--ont_del_frac", type=float, default=None,
                    help="ont mode: fraction of errors that are deletions (r9.4.1 0.35, r10.4.1 0.40)")
    ap.add_argument("--homopolymer_factor", type=float, default=None,
                    help="ont mode: indel probability of a homopolymer run of length L>=3 is "
                         "(1-accuracy)*factor*(L-2) (r9.4.1 4, r10.4.1 5; 0 disables)")
    ap.add_argument("--tsd", nargs="*", default=[], metavar="FAMILY=LEN", help="Override TSD lengths, e.g. IS4=10 Tn3=5")
    ap.add_argument("--seed", type=int, default=1, help="Random seed (default 1)")
    args = ap.parse_args(argv)
    if args.error_model == "ont":
        for k, v in ONT_PROFILES[args.ont_profile].items():
            if getattr(args, k) is None:
                setattr(args, k, v)
    return args


def make_ont_model(args):
    """ONT model with its own random stream, so --seed still fixes the run and
    the main stream (read structure) is untouched by error-model settings."""
    return OntErrorModel(random.Random(f"ont-errors-{args.seed}"), args.read_accuracy_mean,
                         args.read_accuracy_sd, args.ont_sub_frac, args.ont_ins_frac,
                         args.ont_del_frac, args.homopolymer_factor)


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

    try:
        ont = make_ont_model(args) if args.error_model == "ont" else None
    except ValueError as exc:
        sys.exit(f"ONT error model: {exc}")
    os.makedirs(args.out_dir, exist_ok=True)
    elements = []
    if args.novel_mode == "element":
        genome, elements = build_novel_elements(args, genome)
    truth_fields = ["barcode", "read_id", "category", "clone_id", "element", "element_family",
                    "ref_pos", "orientation", "tsd_len", "tsd_seq", "insert_len",
                    "genome_contig", "genome_start", "genome_end", "read_strand", "fragmented", "read_len"]
    if ont is not None:
        truth_fields += ONT_TRUTH_FIELDS
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

        with open(os.path.join(bc_dir, f"{bc}_raw.fastq"), "w", buffering=1 << 22) as fq:
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

                a, b_ = 0, len(mol)
                if rng.random() < args.fragment_rate and len(mol) > 1000:
                    a = rng.randrange(0, len(mol) - 500)
                    b_ = rng.randrange(a + 500, len(mol) + 1)
                    rec["fragmented"] = True
                strand = rng.choice("+-")
                if ont is not None:
                    # errors from the separate ONT stream; the full molecule is
                    # passed so its homopolymer runs can be cached
                    seq, qual, stats = ont.read(mol, a, b_, strand)
                    rec.update(stats)
                    rec.update(read_strand=strand, read_len=len(seq))
                    fq.write(f"@{rid}\n{seq}\n+\n{qual}\n")
                    truth.writerow(rec)
                    continue
                if rec["fragmented"]:
                    mol = mol[a:b_]
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
        if ont is None:
            # keep simulation.json identical to earlier versions in uniform mode
            for k in ("error_model", "ont_profile") + ONT_PARAM_KEYS:
                params.pop(k)
        else:
            for k in ("sub", "ins", "dele", "mean_q"):
                params.pop(k)
            params["ont_model_note"] = ONT_NOTE
        json.dump({**params, "cut_pos_used": cut, "plasmid_name": plasmid_name,
                   "plasmid_len": plen, "tsd_table": tsd_table}, fh, indent=2)
    print(f"Simulated {args.barcodes} barcode(s) x {args.reads} reads into {args.out_dir} (cut site {cut})")


if __name__ == "__main__":
    main()
