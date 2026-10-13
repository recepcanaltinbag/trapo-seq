#!/usr/bin/env python
"""Build a consolidated, non-redundant mobile-genetic-element FASTA for trapo-seq.

Merges third-party element databases (ISfinder IS.fna + taxonomy, TnCentral /
INTEGRALL / ISfinder exports, PICOTA) into one FASTA whose headers follow the
convention of `blast_annot --is_fasta`:

    >{Family}_{Name}

The pipeline splits the first header token at its FIRST underscore into family
and name (src/h_is_stats.py), so Family never contains an underscore (e.g.
ISfinder "IS200/IS605" becomes "IS200-IS605"); Name may contain underscores.
The first header token is unique and there is no description after it.

The source databases are licensed third-party data: the outputs go to a local,
git-ignored directory (default IS_Db/trapo_db/) and must not be committed.

INPUTS (auto-discovered under --src, default IS_Db/DBIS, or given explicitly)
  --isfinder   ISfinder nucleotide FASTA, headers Name_Family_Group
               (e.g. ISPs1_IS5_IS427; names may contain "_", so the header is
               split from the right).
  --taxonomy   ISfinder taxonomy TSV (IS_Name, IS_Family, IS_Group); its group
               replaces "unknown" groups of the FASTA headers and gives the
               family of ISfinder elements found only in the other files.
  --extra      any number of TnCentral / INTEGRALL / PICOTA FASTAs. Header
               styles: "ISFinder_<name>" (ISfinder element), "PICOTA_<name>",
               "INTEGRALL-<cassette>-<accession>" (gene cassette) and
               "<name>-<accession>" (TnCentral: IS variants, Tn, In, ISCR).
               Malformed records - a header glued to the end of the previous
               sequence line, possibly wrapped over two lines and glued to the
               next sequence ("...gaaaacgaca>In1223-KX784502tgtcattt") - are
               repaired; records that cannot be repaired are skipped and
               reported.

OUTPUTS (in --out)
  elements.fasta          selected element types (--types), uppercase.
  elements_all.fasta      every element type.
  elements_is_only.fasta  elements present in the ISfinder FASTA (IS + MITE),
                          deduplicated, with the same IDs as elements.fasta:
                          comparable with runs on ISfinder alone.
  elements.tsv            one row per kept record (all types) with metadata,
                          aliases (identical / reverse-complement-identical
                          sequences merged into it) and notes.
  summary.txt             counts per source / type / family, duplicates,
                          repaired / skipped records, unknown families.

FAMILY MAPPING (header family; the original ISfinder family is kept in the TSV)
  IS from ISfinder / TnCentral IS variants -> ISfinder family, "/" and "_" -> "-"
                                             (IS200/IS605 -> IS200-IS605);
                                             variants (IS26.var1, IS10_p,
                                             ISGNB1.1, ...) inherit the family
                                             of their root name.
  MITE (MITE* names in ISfinder)           -> MITE (--mite_family parent: the
                                             related IS family instead).
  ISCR*                                    -> ISCR
  Tn* distributed by ISfinder (ISFinder_Tn*) -> Tn3 (ISfinder only hosts Tn3-
                                             family transposons)
  Tn* from TnCentral (and other Tn names)  -> Tn
  In*, In_Tn*, InS* (integrons)            -> Integron
  INTEGRALL-* (gene cassettes)             -> Cassette
  compound TnCentral entries whose name
  joins several elements (IS1326_IS1353,
  IS1133_Tn10_IS903B, ISPa12_ISPpu17)      -> Compound
  anything else / IS without taxonomy      -> unknown
  --family_map TSV (name<TAB>family) overrides any of the above.

DEDUPLICATION
  Sequences are uppercased (U->T); whitespace and characters that are not
  IUPAC nucleotide codes are removed (noted). Records with identical sequence
  or identical reverse complement are merged; the kept record is the best
  annotated one: ISfinder FASTA > ISfinder (via other files) > TnCentral >
  PICOTA > INTEGRALL, then known family, known group, non-MITE, shorter name.
  The others are listed as aliases. An "ISFinder_<name>" record whose name is
  in the ISfinder FASTA with a different sequence is an older copy of that
  element and is dropped ("superseded", listed in the summary).

ELEMENT TYPES (--types; default IS,MITE,ISCR)
  IS, MITE, ISCR, transposon, integron, cassette, compound, other  ("all": every type)
  Large elements are not in the default because an IS-sized insert is fully
  covered by every Tn/integron/compound that carries a copy of that IS, and
  blast_annot ranks subjects by covered insert bases, so the label of such an
  insert would depend on BLAST tie order. Subjects > partial_len (8 kb) are
  only accepted when >= partial_threshold (80 %) of them is aligned, which
  protects against long Tn but not against the many shorter ones.

Usage:
  python scripts/build_element_db.py                       # IS_Db/DBIS -> IS_Db/trapo_db
  python scripts/build_element_db.py --types all --out IS_Db/trapo_db_all
  python scripts/build_element_db.py --isfinder IS.fna --taxonomy tax.tsv \\
         --extra tncentral.fa picota.fasta --out mydb

Standard library only.
"""

import argparse
import collections
import datetime
import glob
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SRC = os.path.join(REPO, "IS_Db", "DBIS")
DEFAULT_OUT = os.path.join(REPO, "IS_Db", "trapo_db")

TYPES = ("IS", "MITE", "ISCR", "transposon", "integron", "cassette", "compound", "other")
DEFAULT_TYPES = ("IS", "MITE", "ISCR")
TYPE_ORDER = {t: i for i, t in enumerate(TYPES)}

# Origin of a record, best annotated first.
ORIGINS = ("ISfinder", "ISfinder-copy", "TnCentral", "PICOTA", "INTEGRALL")
ORIGIN_RANK = {o: i for i, o in enumerate(ORIGINS)}

IUPAC = set("ACGTRYKMSWBDHVN")
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")

# Accession-like suffix of TnCentral/INTEGRALL headers: AF071413, NC_002525.1,
# X03216.1.1, AJ6276431, AZ907J2 ...
ACCESSION_RE = re.compile(r"^[A-Z]{1,6}_?[A-Z0-9]*\d[A-Z0-9]*(\.\d+)*$")
PSEUDO_ACCESSIONS = {"ISfinder", "ND"}
# Variant suffixes of TnCentral IS names: IS26.var1, IS10_p, IS1133_p2, ISCfr1.1, ISPme1_surrogate.1
VARIANT_RE = re.compile(r"(\.var\d+|_p\d*|_surrogate(\.\d+)?|\.\d+)$")
ELEMENT_TOKEN_RE = re.compile(r"^(IS|Tn|In|MITE)[A-Za-z0-9]")


#------------------------------------------------------------ sequences -----

def revcomp(seq):
    return seq.translate(COMP)[::-1]


def canonical(seq):
    """Strand-independent key: the smaller of the sequence and its reverse complement."""
    rc = revcomp(seq)
    return seq if seq <= rc else rc


def clean_sequence(raw):
    """Uppercase, U->T, drop whitespace and non-IUPAC characters.

    Returns (seq, notes) where notes lists removed characters and ambiguity
    codes (non-ACGT IUPAC letters, kept)."""
    s = "".join(raw.split()).upper().replace("U", "T")
    removed = collections.Counter(c for c in s if c not in IUPAC)
    if removed:
        s = "".join(c for c in s if c in IUPAC)
    notes = []
    if removed:
        notes.append("removed non-IUPAC characters " +
                     ",".join(f"{c!r}x{n}" for c, n in sorted(removed.items())))
    amb = collections.Counter(c for c in s if c not in "ACGT")
    if amb:
        notes.append("ambiguity codes " + ",".join(f"{c}x{n}" for c, n in sorted(amb.items())))
    return s, notes


#------------------------------------------------------------- FASTA I/O ----

Raw = collections.namedtuple("Raw", "header seq repaired")


def _split_glued_header(text):
    """Split 'In1223-KX784502tgtcattt' into (header, sequence).

    The sequence of the malformed source files is lowercase, so the header
    ends at the last character that is not a lowercase a/c/g/t/n. Returns
    (header, seq) or (None, text) if no plausible header is found."""
    m = re.match(r"^(.*[^acgtn])([acgtn]*)$", text, re.S)
    if not m:
        return None, text
    header, seq = m.group(1), m.group(2)
    if len(header) > 200 or not _plausible_header(header):
        return None, text
    return header, seq


def _plausible_header(header):
    h = header.strip()
    if not h or any(c in h for c in "\n>"):
        return False
    name, acc = _split_accession(h)
    return acc is not None or h.startswith(("ISFinder_", "PICOTA_", "INTEGRALL"))


def read_fasta_repair(path):
    """Parse a FASTA, repairing headers glued into sequence lines.

    Returns (records, problems) where records is a list of Raw and problems a
    list of messages for records that could not be repaired (skipped)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    # Wrap width: most common length of non-final sequence lines (None for
    # single-line sequences). A glued header at the end of a full-width line
    # may continue on the next line.
    widths = collections.Counter()
    for a, b in zip(lines, lines[1:]):
        if a and not a.startswith(">") and ">" not in a and b and not b.startswith(">"):
            widths[len(a)] += 1
    wrap = widths.most_common(1)[0][0] if widths else None

    records, problems = [], []
    header, seq, repaired = None, [], False

    def flush():
        if header is not None:
            records.append(Raw(header, "".join(seq), repaired))

    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\r")
        i += 1
        if line.startswith(">"):
            flush()
            header, seq, repaired = line[1:].strip(), [], False
            continue
        if ">" not in line:
            seq.append(line)
            continue
        # Header(s) glued into a sequence line.
        parts = line.split(">")
        seq.append(parts[0])
        for k, tail in enumerate(parts[1:], 1):
            flush()
            last = k == len(parts) - 1
            text = tail
            if last and wrap and len(line) == wrap and i < len(lines) \
                    and not lines[i].startswith(">") and not re.search(r"[acgtn]$", tail):
                # header wrapped onto the next line: 'In_Tn6234-HG9' + '34082tgtc...'
                text = tail + lines[i].rstrip("\r")
                i += 1
            if last and not re.search(r"[acgtn]$", text):
                # header alone at the end of the line (single-line FASTA)
                h, s = (text.strip(), "") if _plausible_header(text) else (None, text)
            else:
                h, s = _split_glued_header(text)
            if h is None:
                problems.append(f"{os.path.basename(path)} line {i}: cannot repair glued header "
                                f"{text[:60]!r}; following sequence skipped")
                header, seq, repaired = None, [], False
                # skip until the next real header
                while i < len(lines) and not lines[i].startswith(">"):
                    i += 1
                break
            header, seq, repaired = h.strip(), [s], True
    flush()
    return records, problems


def write_fasta(path, records, width=80):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for j in range(0, len(seq), width):
                fh.write(seq[j:j + width] + "\n")


#------------------------------------------------------------ names etc. ----

def sanitize_family(family):
    f = re.sub(r"[/_\s]+", "-", str(family).strip())
    f = re.sub(r"[^A-Za-z0-9.-]+", "", f).strip("-")
    return f or "unknown"


def sanitize_name(name):
    n = re.sub(r"\s+", "", str(name))
    n = re.sub(r"[^A-Za-z0-9._+-]+", "-", n).strip("-")
    return n or "unnamed"


def _split_accession(text):
    """'IS26.var1-GU562437.2' -> ('IS26.var1', 'GU562437.2'); 'ISAba19-ISfinder'
    -> ('ISAba19', ''); no accession -> (text, None)."""
    t = re.sub(r"\s+", "", text)
    if "-" not in t:
        return t, None
    name, acc = t.rsplit("-", 1)
    if not name:
        return t, None
    if acc in PSEUDO_ACCESSIONS:
        return name, ""
    if ACCESSION_RE.match(acc):
        return name, acc
    return t, None


def element_type(name):
    """Element type from a (TnCentral / ISfinder) element name."""
    if name.startswith("MITE"):
        return "MITE"
    if re.match(r"^ISCR\d", name):
        return "ISCR"
    if re.match(r"^In(\d|_|S\d|[A-Z])", name) or name == "In":
        return "integron"
    if name.startswith("Tn"):
        return "transposon"
    if name.startswith("IS"):
        return "IS"
    return "other"


def is_compound(name):
    """True for TnCentral names joining several elements (IS1326_IS1353)."""
    tokens = name.split("_")
    return len(tokens) > 1 and any(ELEMENT_TOKEN_RE.match(t) for t in tokens[1:])


def root_names(name):
    """Candidate ISfinder names of a variant: IS26.var1 -> IS26; ISGNB1.1 -> ISGNB1-1."""
    out = [name]
    n = name
    while True:
        m = VARIANT_RE.search(n)
        if not m or m.start() == 0:
            break
        n = n[:m.start()]
        out.append(n)
    out += [name.replace(".", "-")]
    # left/right copies of composite-transposon ISs share the family: IS10 -> IS10R, IS50L -> IS50R
    for n in list(out):
        if n[-1:] in ("L", "R"):
            out.append(n[:-1] + ("R" if n[-1] == "L" else "L"))
        else:
            out += [n + "R", n + "L"]
    return list(dict.fromkeys(out))


#------------------------------------------------------------- records ------

class Rec:
    __slots__ = ("name", "family", "isf_family", "group", "etype", "origin", "accession",
                 "seq", "sources", "notes", "aliases", "raw_id")

    def __init__(self, name, family, isf_family, group, etype, origin, accession, seq,
                 source, notes, raw_id):
        self.name, self.family, self.isf_family, self.group = name, family, isf_family, group
        self.etype, self.origin, self.accession, self.seq = etype, origin, accession, seq
        self.sources = [source]
        self.notes = list(notes)
        self.aliases = []
        self.raw_id = raw_id

    def rank(self):
        return (ORIGIN_RANK[self.origin], self.family == "unknown", not self.group,
                self.etype == "MITE", len(self.name), self.name)


def load_taxonomy(path):
    tax = {}
    if not path:
        return tax
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or f[0] == "IS_Name":
                continue
            fam = f[1].strip()
            grp = f[2].strip() if len(f) > 2 else ""
            tax[f[0].strip()] = (fam, "" if grp.lower() == "unknown" else grp)
    return tax


def parse_isfinder(path, tax, mite_family, stats):
    """ISfinder IS.fna: Name_Family_Group."""
    recs, problems = read_fasta_repair(path)
    out = []
    src = os.path.basename(path)
    for r in recs:
        stats["parsed", src] += 1
        parts = r.header.split()[0].rsplit("_", 2)
        if len(parts) != 3:
            name, fam, grp = r.header.split()[0], "unknown", ""
        else:
            name, fam, grp = parts
        grp = "" if grp.lower() == "unknown" else grp
        if name in tax:
            tfam, tgrp = tax[name]
            if tfam and tfam != fam:
                stats["taxonomy_family_conflict"] += 1
            grp = grp or tgrp
        seq, notes = clean_sequence(r.seq)
        if r.repaired:
            notes.append("repaired malformed record")
            stats["repaired", src] += 1
        etype = "MITE" if name.startswith("MITE") else "IS"
        family = sanitize_family(fam)
        if etype == "MITE" and mite_family == "MITE":
            family = "MITE"
        out.append(Rec(sanitize_name(name), family, fam, grp, etype, "ISfinder", "", seq, src,
                       notes, r.header))
    return out, problems


def parse_extra(path, tax, isf_names, mite_family, stats):
    """TnCentral / INTEGRALL / PICOTA style FASTA."""
    recs, problems = read_fasta_repair(path)
    src = os.path.basename(path)
    out = []
    for r in recs:
        stats["parsed", src] += 1
        h = r.header.strip()
        notes = []
        accession = ""
        if h.startswith("ISFinder_"):
            name, origin = h[len("ISFinder_"):].split()[0], "ISfinder-copy"
        elif h.startswith("PICOTA_"):
            name, origin = h[len("PICOTA_"):].split()[0], "PICOTA"
        elif h.startswith("INTEGRALL"):
            name, acc = _split_accession(re.sub(r"^INTEGRALL[-_]?", "", h))
            accession, origin = acc or "", "INTEGRALL"
        else:
            name, acc = _split_accession(h)
            accession = acc or ""
            origin = "TnCentral" if "tncentral" in src.lower() else "PICOTA"
        if origin == "INTEGRALL":
            etype = "cassette"
        else:
            etype = element_type(name)
            if etype == "IS" and origin in ("TnCentral", "PICOTA") and is_compound(name):
                etype = "compound"
        isf_family, grp = "", ""
        for cand in root_names(name):
            if cand in tax:
                isf_family, grp = tax[cand]
                break
        if etype == "MITE":
            family = "MITE" if mite_family == "MITE" else sanitize_family(isf_family or "unknown")
        elif etype == "ISCR":
            family = "ISCR"
        elif etype == "transposon":
            family = "Tn3" if origin == "ISfinder-copy" else "Tn"
            if origin == "ISfinder-copy":
                notes.append("family Tn3 inferred (ISfinder hosts Tn3-family transposons)")
        elif etype == "integron":
            family = "Integron"
        elif etype == "cassette":
            family = "Cassette"
        elif etype == "compound":
            family = "Compound"
        else:
            family = sanitize_family(isf_family) if isf_family else "unknown"
        if etype == "IS" and not isf_family:
            notes.append("not in ISfinder taxonomy")
        seq, cn = clean_sequence(r.seq)
        notes += cn
        if r.repaired:
            notes.append("repaired malformed record")
            stats["repaired", src] += 1
        rec = Rec(sanitize_name(name), family, isf_family, grp, etype, origin, accession, seq,
                  src, notes, h)
        if origin == "ISfinder-copy" and name in isf_names:
            rec.notes.append("superseded")
        out.append(rec)
    return out, problems


#--------------------------------------------------------------- build ------

def discover(src):
    """Default input files under `src` (and src/source)."""
    def first(*patterns):
        for p in patterns:
            hits = sorted(glob.glob(os.path.join(src, p)))
            if hits:
                return hits[0]
        return None
    isf = first("IS.fna", "source/IS.fna")
    tax = first("isfinder_taxonomy.tsv", "source/isfinder_taxonomy.tsv")
    extras = []
    for p in ("*.fa", "*.fasta", "source/*.fa", "source/*.fasta"):
        for f in sorted(glob.glob(os.path.join(src, p))):
            real = os.path.realpath(f)
            if real == os.path.realpath(isf or "") or real in {os.path.realpath(e) for e in extras}:
                continue
            extras.append(f)
    return isf, tax, extras


def build(isfinder=None, taxonomy=None, extras=(), family_map=None, mite_family="MITE",
          min_len=0):
    """Parse, clean and deduplicate. Returns (kept records, report dict)."""
    stats = collections.Counter()
    problems = []
    tax = load_taxonomy(taxonomy)
    fmap = {}
    if family_map:
        with open(family_map) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0] and not f[0].startswith("#"):
                    fmap[f[0].strip()] = sanitize_family(f[1])

    all_recs = []
    isf_names = set()
    if isfinder:
        recs, pr = parse_isfinder(isfinder, tax, mite_family, stats)
        problems += pr
        all_recs += recs
        isf_names = {r.raw_id.rsplit("_", 2)[0] for r in recs}
    for path in extras:
        recs, pr = parse_extra(path, tax, isf_names, mite_family, stats)
        problems += pr
        all_recs += recs

    superseded, short, empty = [], [], []
    by_key = {}
    order = []
    for r in all_recs:
        if fmap and r.name in fmap:
            r.family = fmap[r.name]
            r.notes.append("family from --family_map")
        if not r.seq:
            empty.append(r)
            continue
        if len(r.seq) < min_len:
            short.append(r)
            continue
        key = canonical(r.seq)
        by_key.setdefault(key, []).append(r)
        if len(by_key[key]) == 1:
            order.append(key)

    kept = []
    dup_identical = dup_rc = 0
    for key in order:
        group = by_key[key]
        # An ISfinder copy that is identical to an IS.fna record is a plain duplicate;
        # one whose name is in IS.fna with a different sequence is an old version.
        live = [r for r in group if "superseded" not in r.notes]
        if not live:
            superseded.extend(group)
            continue
        group_sup = [r for r in group if "superseded" in r.notes]
        best = min(live, key=Rec.rank)
        best.sources = sorted({s for r in group for s in r.sources})
        for r in group:
            if r is best:
                continue
            if r.seq == best.seq:
                dup_identical += 1
            else:
                dup_rc += 1
            if r in group_sup:
                r.notes.remove("superseded")
            if r.name != best.name and r.name not in best.aliases:
                best.aliases.append(r.name)
            if r.accession and not best.accession:
                best.accession = r.accession
            if r.family != best.family and r.etype == best.etype and r.family != "unknown" \
                    and best.family != "unknown":
                best.notes.append(f"identical to {r.name} ({r.family})")
        if any(r.seq != best.seq for r in group):
            best.notes.append("reverse-complement duplicate(s) merged")
        kept.append(best)

    # Unique IDs: Family_Name; a clash gets the accession, then a counter.
    by_id = {}
    kept.sort(key=lambda r: (TYPE_ORDER[r.etype], ORIGIN_RANK[r.origin], r.family, r.name))
    for r in kept:
        base = f"{r.family}_{r.name}"
        if base in by_id and r.accession:
            base = f"{base}-{sanitize_name(r.accession)}"
        rid, n = base, 1
        while rid in by_id:
            n += 1
            rid = f"{base}.{n}"
        if rid != f"{r.family}_{r.name}":
            r.notes.append(f"name clash with {r.family}_{r.name}; renamed")
        by_id[rid] = r
        r.raw_id = rid  # final id
    report = dict(stats=stats, problems=problems, superseded=superseded, short=short,
                  empty=empty, dup_identical=dup_identical, dup_rc=dup_rc, n_input=len(all_recs),
                  inputs=dict(isfinder=isfinder, taxonomy=taxonomy, extras=list(extras),
                              family_map=family_map))
    return kept, report


def parse_types(text):
    if text is None or text.strip().lower() == "default":
        return set(DEFAULT_TYPES)
    if text.strip().lower() == "all":
        return set(TYPES)
    lookup = {t.lower(): t for t in TYPES}
    out = set()
    for t in text.split(","):
        t = t.strip()
        if not t:
            continue
        if t.lower() not in lookup:
            raise SystemExit(f"--types: unknown type {t!r}; choose from {', '.join(TYPES)} or all")
        out.add(lookup[t.lower()])
    return out


TSV_COLUMNS = ["id", "name", "family", "isfinder_family", "group", "element_type", "origin",
               "sources", "accession", "length", "selected", "aliases", "notes"]


def write_outputs(kept, report, out_dir, types):
    os.makedirs(out_dir, exist_ok=True)
    selected = [r for r in kept if r.etype in types]
    write_fasta(os.path.join(out_dir, "elements.fasta"), [(r.raw_id, r.seq) for r in selected])
    write_fasta(os.path.join(out_dir, "elements_all.fasta"), [(r.raw_id, r.seq) for r in kept])
    isf_src = os.path.basename(report["inputs"]["isfinder"]) if report["inputs"]["isfinder"] else None
    is_only = [r for r in kept if isf_src and isf_src in r.sources and r.etype in ("IS", "MITE")]
    write_fasta(os.path.join(out_dir, "elements_is_only.fasta"), [(r.raw_id, r.seq) for r in is_only])
    with open(os.path.join(out_dir, "elements.tsv"), "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in kept:
            vals = [r.raw_id, r.name, r.family, r.isf_family, r.group, r.etype, r.origin,
                    ";".join(r.sources), r.accession, len(r.seq),
                    "yes" if r.etype in types else "no", ";".join(r.aliases),
                    "; ".join(dict.fromkeys(r.notes))]
            fh.write("\t".join(str(v).replace("\t", " ") for v in vals) + "\n")
    summary = format_summary(kept, report, types, selected, is_only)
    with open(os.path.join(out_dir, "summary.txt"), "w") as fh:
        fh.write(summary)
    return summary


def format_summary(kept, report, types, selected, is_only):
    st = report["stats"]
    L = [f"trapo-seq element DB  ({datetime.datetime.now():%Y-%m-%d %H:%M})", ""]
    L.append("Inputs:")
    inp = report["inputs"]
    L.append(f"  ISfinder FASTA : {inp['isfinder']}")
    L.append(f"  taxonomy       : {inp['taxonomy']}")
    for e in inp["extras"]:
        L.append(f"  extra          : {e}")
    if inp["family_map"]:
        L.append(f"  family map     : {inp['family_map']}")
    L.append("")
    L.append("Records parsed per file (repaired malformed records):")
    for (k, src), n in sorted((k, n) for k, n in st.items() if isinstance(k, tuple) and k[0] == "parsed"):
        L.append(f"  {src:50s} {n:6d}  (repaired {st['repaired', src]})")
    L.append(f"  total input records: {report['n_input']}")
    L.append("")
    L.append(f"Skipped: {len(report['problems'])} unrepairable, {len(report['empty'])} empty, "
             f"{len(report['short'])} shorter than --min_len")
    for p in report["problems"]:
        L.append(f"  {p}")
    for r in report["empty"] + report["short"]:
        L.append(f"  {r.sources[0]}: {r.name} ({len(r.seq)} bp)")
    L.append(f"Duplicates merged: {report['dup_identical']} identical + {report['dup_rc']} "
             f"reverse-complement identical")
    sup = report["superseded"]
    L.append(f"Superseded ISfinder copies dropped (name in ISfinder FASTA, different sequence): {len(sup)}")
    for r in sup:
        L.append(f"  {r.sources[0]}: {r.name} ({len(r.seq)} bp)")
    L.append("")
    L.append(f"Kept (non-redundant) records: {len(kept)}")
    L.append(f"Selected types for elements.fasta: {','.join(t for t in TYPES if t in types)} "
             f"-> {len(selected)} records, {sum(len(r.seq) for r in selected):,} bp")
    L.append(f"elements_is_only.fasta (in ISfinder FASTA, IS+MITE): {len(is_only)} records")
    L.append("")
    L.append("Per type (kept / selected; length min-median-max):")
    for t in TYPES:
        rs = [r for r in kept if r.etype == t]
        if not rs:
            continue
        ls = sorted(len(r.seq) for r in rs)
        L.append(f"  {t:11s} {len(rs):6d} / {sum(1 for r in rs if r.etype in types):6d}   "
                 f"{ls[0]}-{ls[len(ls) // 2]}-{ls[-1]} bp; > 8 kb: {sum(l > 8000 for l in ls)}")
    L.append("")
    L.append("Per origin of the kept record:")
    for o in ORIGINS:
        n = sum(1 for r in kept if r.origin == o)
        if n:
            L.append(f"  {o:14s} {n:6d}")
    L.append("")
    L.append("Records per source file (after merging; a record can come from several):")
    srcs = collections.Counter(s for r in kept for s in r.sources)
    for s, n in srcs.most_common():
        L.append(f"  {s:50s} {n:6d}")
    L.append("")
    L.append("Per family (all kept records):")
    fams = collections.Counter(r.family for r in kept)
    for f, n in sorted(fams.items(), key=lambda x: (-x[1], x[0])):
        L.append(f"  {f:14s} {n:6d}")
    unk = [r for r in kept if r.family == "unknown"]
    L.append("")
    L.append(f"Family unknown ({len(unk)}): " + ", ".join(sorted(r.name for r in unk)))
    L.append("")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="See the module docstring for the family mapping and dedup rules.")
    ap.add_argument("--src", default=DEFAULT_SRC,
                    help="directory to auto-discover inputs in (default IS_Db/DBIS)")
    ap.add_argument("--isfinder", help="ISfinder nucleotide FASTA (Name_Family_Group headers)")
    ap.add_argument("--taxonomy", help="ISfinder taxonomy TSV (IS_Name, IS_Family, IS_Group)")
    ap.add_argument("--extra", nargs="*", default=None,
                    help="TnCentral / INTEGRALL / PICOTA FASTA files (default: all other *.fa/*.fasta in --src)")
    ap.add_argument("-o", "--out", default=DEFAULT_OUT, help="output directory (default IS_Db/trapo_db)")
    ap.add_argument("--types", default="default",
                    help=f"element types in elements.fasta, comma-separated from {','.join(TYPES)}, "
                         f"or 'all' (default: {','.join(DEFAULT_TYPES)})")
    ap.add_argument("--family_map", help="TSV name<TAB>family overriding the family of named elements")
    ap.add_argument("--mite_family", choices=("MITE", "parent"), default="MITE",
                    help="header family of MITEs: MITE (default) or the related IS family")
    ap.add_argument("--min_len", type=int, default=0, help="drop sequences shorter than this (default 0)")
    a = ap.parse_args(argv)

    types = parse_types(a.types)
    isf, tax, extras = discover(a.src) if os.path.isdir(a.src) else (None, None, [])
    isf = a.isfinder or isf
    tax = a.taxonomy or tax
    extras = a.extra if a.extra is not None else extras
    if not isf and not extras:
        ap.error(f"no input FASTA found (looked in {a.src}); give --isfinder and/or --extra")
    kept, report = build(isf, tax, extras, a.family_map, a.mite_family, a.min_len)
    summary = write_outputs(kept, report, a.out, types)
    print(summary)
    print(f"Outputs written to {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
