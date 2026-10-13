"""Shared pytest fixtures for trapo-seq."""
import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DATA = os.path.join(REPO, "test-data")
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "simulate"))


def have(tool):
    return shutil.which(tool) is not None


requires_tools = pytest.mark.skipif(
    not all(have(t) for t in ("minimap2", "samtools", "blastn", "makeblastdb")),
    reason="external tools (minimap2, samtools, blastn, makeblastdb) not in PATH",
)


@pytest.fixture(scope="session")
def sim_dir(tmp_path_factory):
    """A small simulated traposome dataset with ground truth.

    The directory name contains a space on purpose: paths with spaces must work.
    """
    out = tmp_path_factory.mktemp("sim data")
    subprocess.run(
        [sys.executable, os.path.join(REPO, "simulate", "simulate_traposome.py"),
         "-p", os.path.join(TEST_DATA, "plasmid.fasta"),
         "-g", os.path.join(TEST_DATA, "genome.fasta"),
         "-i", os.path.join(TEST_DATA, "IS_curated.fasta"),
         "-o", str(out), "--barcodes", "2", "--reads", "150", "--seed", "11"],
        check=True,
    )
    return out
