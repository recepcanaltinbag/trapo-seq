"""Smoke tests for the command-line interface."""
import os
import subprocess
import sys

import pytest

from conftest import REPO

CLI = os.path.join(REPO, "trapo-seq.py")
SUBCOMMANDS = ["read_histogram", "filter", "map", "insert_finder", "insert_finder_batch",
               "blast_annot", "blast_annot_batch", "is_stat", "heatmap", "dr_finder",
               "dr_logo", "in_del_plot", "kde_mobile", "map_dist", "test"]


def run(args, cwd):
    return subprocess.run([sys.executable, CLI] + args, cwd=cwd, capture_output=True, text=True)


def test_version_from_any_directory(tmp_path):
    from src import __version__
    res = run(["--version"], cwd=tmp_path)
    assert res.returncode == 0
    assert __version__ in res.stdout


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_subcommand_help(sub, tmp_path):
    res = run([sub, "--help"], cwd=tmp_path)
    assert res.returncode == 0, res.stderr
    assert "usage" in res.stdout


def test_no_invalid_escape_warnings(tmp_path):
    res = subprocess.run([sys.executable, "-W", "error::SyntaxWarning", "-c",
                          f"import ast; ast.parse(open({CLI!r}).read())"],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
