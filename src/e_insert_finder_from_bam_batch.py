"""Batch version of e_insert_finder_from_bam: run the insertion finder on the
``03_*.bam`` file of every sample folder in a data directory and write
``<folder>/insertion_from_bam.tab``.

See e_insert_finder_from_bam for the algorithm and output columns.

Recep Can Altınbağ, 12 11 2024, v0.3
"""

import logging
import os

from .e_insert_finder_from_bam import OUTPUT_NAME, main_insert_finder_from_bam

logger = logging.getLogger(__name__)


def find_file_with_prefix(directory, prefix, enfix):
    """First file (sorted by name) in ``directory`` that starts with ``prefix``
    and ends with ``enfix``, or None."""
    for filename in sorted(os.listdir(directory)):
        if filename.startswith(prefix) and filename.endswith(enfix):
            return os.path.join(directory, filename)
    return None


def main_insert_finder_from_bam_batch(data_dir, insertion_threshold=500, plasmid_fasta=None,
                                      enzyme=None, cut_site=None, overhang=None, adapter_fasta=None):
    """Returns {folder: list of Insertion records} for the processed folders.

    plasmid_fasta / enzyme / cut_site / overhang / adapter_fasta: artefact
    checks (see e_insert_finder_from_bam.main_insert_finder_from_bam)."""
    folders = sorted(f for f in os.listdir(data_dir)
                     if os.path.isdir(os.path.join(data_dir, f)))
    results = {}
    for folder in folders:
        folder_path = os.path.join(data_dir, folder)
        bamfile = find_file_with_prefix(folder_path, "03_", ".bam")
        if bamfile is None:
            logger.warning("No 03_*.bam file in %s, skipping.", folder)
            continue
        logger.info("Processing %s (%s)", folder, bamfile)
        output_file = os.path.join(folder_path, OUTPUT_NAME)
        results[folder] = main_insert_finder_from_bam(
            bamfile, output_file, insertion_threshold, plasmid_fasta=plasmid_fasta, enzyme=enzyme,
            cut_site=cut_site, overhang=overhang, adapter_fasta=adapter_fasta)
    return results
