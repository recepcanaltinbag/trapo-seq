#-------------------------------------
# Batch version of g_blast_annot: runs main_annot on every sample folder of
# `data_dir` (folder <sample>/ with <sample>.fasta and insertion_from_bam.tab)
# and writes <sample>/best_alignment.tab. BLAST DBs are built once and shared.
#
# After all samples, sites (label + plasmid position) are compared across
# samples: <data_dir>/sites_across_barcodes.tsv lists the reads per sample of
# every site, and a site of a sample with < bleed_pct % of the reads it has in
# another sample is a possible barcode bleed-through (demultiplexing error,
# ARTEFACTS A3); its calls get the caution flag possible_bleed_through in the
# Artefact column (the Note is never changed: hotspots shared by samples are
# real too).
#
# Recep Can Altınbağ, 07 11 2024, v0.2; v0.3 reuses g_blast_annot (audit B19, B38)
#-------------------------------------

import csv
import os
import shutil

from src import artefacts as art
from src.g_blast_annot import DEFAULT_EVALUE, main_annot, sites_path_for

CROSS_SITES_NAME = "sites_across_barcodes.tsv"
DEFAULT_BLEED_PCT = 5.0

TAB_FILE_NAME = "insertion_from_bam.tab"
OUTPUT_NAME = "best_alignment.tab"


def find_file_with_prefix(directory, prefix):
    """Return os.path.join(directory, prefix) if that file exists, else None."""
    path = os.path.join(directory, prefix)
    return path if os.path.isfile(path) else None


def main_annot_batch(data_dir, genom_fasta, is_fasta=None, temp_dir=None,
                     threshold=70, threads=2, debug=False, keep_temp=False,
                     partial_threshold=80, partial_len=8000,
                     evalue=DEFAULT_EVALUE, task="megablast", homology_action="relabel",
                     artefact_action="flag", bleed_pct=DEFAULT_BLEED_PCT):
    """Annotate every sample folder in data_dir. temp_dir defaults to
    <data_dir>/temp. Returns {sample: {Query ID: Note}}.

    artefact_action: see g_blast_annot.main_annot. bleed_pct: a site with less
    than this % of the reads it has in another sample is flagged
    possible_bleed_through (sites_across_barcodes.tsv)."""
    if temp_dir is None:
        temp_dir = os.path.join(data_dir, "temp")
    temp_real = os.path.realpath(temp_dir)
    temp_existed = os.path.isdir(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)
    db_dir = os.path.join(temp_dir, "blastdb")

    folders = sorted(f for f in os.listdir(data_dir)
                     if os.path.isdir(os.path.join(data_dir, f)))
    print('\n')
    print(folders)
    results = {}
    try:
        for folder in folders:
            folder_path = os.path.join(data_dir, folder)
            if os.path.realpath(folder_path) == temp_real:
                continue
            print('\nProcessing: ', folder)
            output_file = os.path.join(folder_path, OUTPUT_NAME)
            mapped_fasta = os.path.join(folder_path, f"{folder}.fasta")
            tab_file = find_file_with_prefix(folder_path, TAB_FILE_NAME)
            if tab_file is None:
                print("No tab insertions file in ", folder, ", skipping!")
                continue
            print('Corresponding tab file: ', tab_file)
            if not os.path.exists(mapped_fasta):
                print("No mapped fasta file in ", folder, ", skipping!")
                continue
            results[folder] = main_annot(
                tab_file, mapped_fasta, genom_fasta, is_fasta,
                os.path.join(temp_dir, folder), output_file, threshold, threads,
                debug, keep_temp, partial_threshold, partial_len,
                evalue=evalue, db_dir=db_dir, task=task, homology_action=homology_action,
                artefact_action=artefact_action)
    finally:
        if not keep_temp:
            shutil.rmtree(db_dir, ignore_errors=True)
            if not temp_existed:
                shutil.rmtree(temp_dir, ignore_errors=True)
    cross_barcode(data_dir, sorted(results), bleed_pct)
    return results


def _read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cross_barcode(data_dir, samples, bleed_pct=DEFAULT_BLEED_PCT):
    """Write <data_dir>/sites_across_barcodes.tsv and add possible_bleed_through
    to the Artefact column of the affected calls. Returns the table rows."""
    per_sample, plen = {}, None
    for sample in samples:
        path = sites_path_for(os.path.join(data_dir, sample, OUTPUT_NAME))
        if not os.path.isfile(path):
            continue
        rows = _read_tsv(path)
        for r in rows:
            for k in ("position", "n_reads"):
                r[k] = int(r[k])
        per_sample[sample] = rows
        cut = art.cut_site_path_for(os.path.join(data_dir, sample, TAB_FILE_NAME))
        if plen is None and os.path.isfile(cut):
            for line in open(cut):
                k, _, v = line.rstrip("\n").partition("\t")
                if k == "plasmid_len" and v.strip().isdigit():
                    plen = int(v)
    out = os.path.join(data_dir, CROSS_SITES_NAME)
    if not per_sample:
        return []
    table = art.cross_barcode_sites(per_sample, plen, bleed_pct=bleed_pct)
    cols = ["label", "position", "n_samples", "total_reads", "possible_bleed_through"] + \
           [f"reads_{s}" for s in sorted(per_sample)] + ["site_ids"]
    art.write_rows(table, out, cols)
    # flag the calls of bleed-through sites
    bleed = {}
    for row in table:
        if row["possible_bleed_through"] == "-":
            continue
        for item in row["site_ids"].split(","):
            sample, sid = item.split(":", 1)
            if sample in row["possible_bleed_through"].split(","):
                bleed.setdefault(sample, set()).add(sid)
    n = 0
    for sample, sids in bleed.items():
        qids = set()
        for r in per_sample[sample]:
            if r["site_id"] in sids:
                qids.update(r["query_ids"].split(","))
        path = os.path.join(data_dir, sample, OUTPUT_NAME)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        if "Artefact" not in fields:
            continue
        for r in rows:
            if r["Query ID"] in qids:
                v = r["Artefact"]
                if v in ("", "-"):
                    r["Artefact"] = "possible_bleed_through"
                elif "possible_bleed_through" not in v:
                    r["Artefact"] = v + ",possible_bleed_through"
                n += 1
        art.write_rows(rows, path, fields)
    if n:
        print(f"Cross-barcode sites: {n} call(s) at sites with < {bleed_pct} % of the reads of another "
              f"sample flagged possible_bleed_through ({out})")
    return table
