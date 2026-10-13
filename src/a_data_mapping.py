"""Map the reads of every sample folder to the plasmid and the genome.

The work is done by ``scripts/01_data_prep_from_py.sh`` (minimap2 + samtools),
run once with ``data_dir`` as working directory. Each sample folder
``<name>/`` must contain ``<name>.fastq``.

Outputs per folder: 02_mapped_to_plasmid.bam,
03_sorted_mapped_to_plasmid.bam(.bai), 04_mapped_to_plasmid.fastq,
07_sorted_aligned_to_other_genome.bam(.bai), 08_mapping_stats.txt, <name>.fasta

The plasmid is treated as circular (AUDIT B09): reads are also mapped to the
plasmid rotated by half its length and, per read, the mapping that covers more
of the read is kept, converted back to the original plasmid coordinates (see
src/circular_mapping.py). The BAM files therefore stay ordinary BAMs against
the plasmid FASTA as given. ``circular=False`` restores the linear mapping.

Plasmid-genome shared segments (GitHub issue #5): the plasmid is compared with
the genome first; the segments are written to
``<data_dir>/plasmid_genome_homology.tsv`` and a warning is logged for
segments of at least ``insertion_threshold`` bp, which can create false
insertion calls (see src/plasmid_genome_homology.py). 08_mapping_stats.txt is
corrected for reads that align to the genome only through a shared segment.

Requires minimap2 and samtools in PATH (and pysam for the circular merge,
run with the current Python interpreter).

Recep Can Altınbağ, 11 11 2024, v0.1
"""

import collections
import logging
import os
import subprocess
import sys

from . import plasmid_genome_homology as homology

logger = logging.getLogger(__name__)

SCRIPT_NAME = "01_data_prep_from_py.sh"


def script_path():
    """Absolute path of the mapping shell script (independent of the cwd)."""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(repo, "scripts", SCRIPT_NAME)


class MappingError(subprocess.CalledProcessError):
    """The mapping script exited with a non-zero status."""

    def __str__(self):
        msg = super().__str__()
        if self.stderr:
            msg += "\nLast lines of stderr:\n" + self.stderr
        return msg


def main_mapping(data_dir, query_fasta, genome_fasta, force=False, threads=4, circular=True,
                 homology_min_len=homology.DEFAULT_MIN_LEN,
                 homology_min_identity=homology.DEFAULT_MIN_IDENTITY,
                 insertion_threshold=homology.DEFAULT_WARN_LEN):
    """Run the mapping script on every sample folder in ``data_dir``.

    ``homology_min_len`` / ``homology_min_identity``: shortest / least similar
    plasmid-genome shared segment listed in plasmid_genome_homology.tsv;
    ``insertion_threshold``: segments at least this long are reported with a
    warning (use the insert_finder threshold).

    Returns the list of shared segments (dicts, see plasmid_genome_homology),
    or None when the table was not written.

    Raises FileNotFoundError for missing inputs and MappingError if any
    mapping step fails (no partial outputs are left behind in that case).
    """
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"input directory not found: {data_dir}")
    for path in (query_fasta, genome_fasta):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"FASTA file not found: {path}")

    script = script_path()
    args = ["bash", script, os.path.abspath(query_fasta), os.path.abspath(genome_fasta)]
    if force:
        args.append("--force")
    args += ["--threads", str(max(1, int(threads)))]
    if not circular:
        args.append("--linear")
    args += ["--homology-min-len", str(int(homology_min_len)),
             "--homology-min-identity", str(float(homology_min_identity)),
             "--insertion-threshold", str(int(insertion_threshold))]
    logger.info("Running %s in %s", script, data_dir)

    tail = collections.deque(maxlen=30)
    with subprocess.Popen(
        args,
        cwd=data_dir,
        text=True,
        stderr=subprocess.PIPE,
        env={**os.environ, "PYTHONUNBUFFERED": "1", "TRAPOSEQ_PYTHON": sys.executable},
    ) as process:
        for line in process.stderr:
            line = line.rstrip("\n")
            tail.append(line)
            logger.info("mapping: %s", line)

    if process.returncode != 0:
        raise MappingError(process.returncode, args, stderr="\n".join(tail))
    logger.info("Mapping finished successfully.")

    table = os.path.join(data_dir, homology.HOMOLOGY_NAME)
    if not os.path.isfile(table):
        return None
    segments = homology.read_homology(table)
    warnings, infos = homology.homology_messages(segments, insertion_threshold)
    for w in warnings:
        logger.warning("Plasmid-genome homology: %s (table: %s)", w, table)
    for i in infos:
        logger.info("Plasmid-genome homology: %s", i)
    return segments
