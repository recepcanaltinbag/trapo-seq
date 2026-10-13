"""Map the reads of every sample folder to the plasmid and the genome.

The work is done by ``scripts/01_data_prep_from_py.sh`` (minimap2 + samtools),
run once with ``data_dir`` as working directory. Each sample folder
``<name>/`` must contain ``<name>.fastq``.

Outputs per folder: 01_aligned_to_plasmid.sam, 02_mapped_to_plasmid.bam,
03_sorted_mapped_to_plasmid.bam(.bai), 04_mapped_to_plasmid.fastq,
05_aligned_to_other_genome.sam, 06_aligned_to_other_genome.bam,
07_sorted_aligned_to_other_genome.bam(.bai), 08_mapping_stats.txt, <name>.fasta

Requires minimap2 and samtools in PATH.

Recep Can Altınbağ, 11 11 2024, v0.1
"""

import collections
import logging
import os
import subprocess

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


def main_mapping(data_dir, query_fasta, genome_fasta, force=False, threads=4):
    """Run the mapping script on every sample folder in ``data_dir``.

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
    logger.info("Running %s in %s", script, data_dir)

    tail = collections.deque(maxlen=30)
    with subprocess.Popen(
        args,
        cwd=data_dir,
        text=True,
        stderr=subprocess.PIPE,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    ) as process:
        for line in process.stderr:
            line = line.rstrip("\n")
            tail.append(line)
            logger.info("mapping: %s", line)

    if process.returncode != 0:
        raise MappingError(process.returncode, args, stderr="\n".join(tail))
    logger.info("Mapping finished successfully.")
