#-------------------------------------
# Batch version of g_blast_annot: runs main_annot on every sample folder of
# `data_dir` (folder <sample>/ with <sample>.fasta and insertion_from_bam.tab)
# and writes <sample>/best_alignment.tab. BLAST DBs are built once and shared.
#
# Recep Can Altınbağ, 07 11 2024, v0.2; v0.3 reuses g_blast_annot (audit B19, B38)
#-------------------------------------

import os
import shutil

from src.g_blast_annot import main_annot, DEFAULT_EVALUE

TAB_FILE_NAME = "insertion_from_bam.tab"
OUTPUT_NAME = "best_alignment.tab"


def find_file_with_prefix(directory, prefix):
    """Return os.path.join(directory, prefix) if that file exists, else None."""
    path = os.path.join(directory, prefix)
    return path if os.path.isfile(path) else None


def main_annot_batch(data_dir, genom_fasta, is_fasta=None, temp_dir=None,
                     threshold=70, threads=2, debug=False, keep_temp=False,
                     partial_threshold=80, partial_len=8000,
                     evalue=DEFAULT_EVALUE, task="megablast"):
    """Annotate every sample folder in data_dir. temp_dir defaults to
    <data_dir>/temp. Returns {sample: {Query ID: Note}}."""
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
                evalue=evalue, db_dir=db_dir, task=task)
    finally:
        if not keep_temp:
            shutil.rmtree(db_dir, ignore_errors=True)
            if not temp_existed:
                shutil.rmtree(temp_dir, ignore_errors=True)
    return results
