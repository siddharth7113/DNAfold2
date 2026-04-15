"""Parity checks between Python stage tools and legacy C tools."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from dnafold2.stage_tools import (
    convert_conf_to_pdb,
    extract_min_conformations,
    generate_initial_ch_dat,
)


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "stage_tools"
REPO_ROOT = Path(__file__).resolve().parent.parent


def _have_compiler(name: str) -> bool:
    return shutil.which(name) is not None


def _drop_legacy_seq_initial_trailing_garbage(text: str) -> str:
    """Drop known trailing garbage line emitted by legacy seq_initial.c."""
    lines = text.splitlines()
    if lines and lines[-1].startswith("0 0 "):
        lines = lines[:-1]
    return "\n".join(lines) + ("\n" if lines else "")


@pytest.mark.contract
def test_seq_initial_python_matches_c(tmp_path: Path) -> None:
    if not _have_compiler("g++"):
        pytest.skip("g++ compiler not available")

    src = FIXTURE_ROOT / "seq_initial"
    c_file = REPO_ROOT / "src" / "initial" / "seq_initial.c"

    # Python output
    py_out = tmp_path / "py_ch.dat"
    generate_initial_ch_dat(src / "seq.dat", src / "initial.dat", py_out)

    # C output
    c_dir = tmp_path / "c_seq_initial"
    c_dir.mkdir()
    shutil.copy2(c_file, c_dir / "seq_initial.c")
    shutil.copy2(src / "seq.dat", c_dir / "seq.dat")
    shutil.copy2(src / "initial.dat", c_dir / "initial.dat")

    subprocess.run(["g++", "seq_initial.c", "-o", "seq_initial"], cwd=c_dir, check=True)
    subprocess.run(["./seq_initial"], cwd=c_dir, check=True)

    c_out = c_dir / "ch.dat"
    assert c_out.exists()
    c_text = _drop_legacy_seq_initial_trailing_garbage(c_out.read_text())
    assert py_out.read_text() == c_text


@pytest.mark.contract
def test_a_state_python_matches_c(tmp_path: Path) -> None:
    if not _have_compiler("gcc"):
        pytest.skip("gcc compiler not available")

    src = FIXTURE_ROOT / "min_extract"
    c_file = REPO_ROOT / "src" / "scoring" / "A_state.c"

    py_out = tmp_path / "py_min_conf.dat"
    extract_min_conformations(src / "ch.dat", src / "min.dat", src / "conf_0.dat", py_out)

    c_dir = tmp_path / "c_a_state"
    c_dir.mkdir()
    shutil.copy2(c_file, c_dir / "A_state.c")
    shutil.copy2(src / "ch.dat", c_dir / "ch.dat")
    shutil.copy2(src / "min.dat", c_dir / "min.dat")
    shutil.copy2(src / "conf_0.dat", c_dir / "conf_0.dat")

    subprocess.run(["gcc", "-Wall", "A_state.c", "-o", "A_state", "-lm"], cwd=c_dir, check=True)
    subprocess.run(["./A_state"], cwd=c_dir, check=True)

    c_out = c_dir / "min_conf.dat"
    assert c_out.exists()
    assert py_out.read_text() == c_out.read_text()


@pytest.mark.contract
def test_tc_python_matches_c(tmp_path: Path) -> None:
    if not _have_compiler("g++"):
        pytest.skip("g++ compiler not available")

    src = FIXTURE_ROOT / "pdb_convert"
    c_file = REPO_ROOT / "src" / "scoring" / "tc.c"

    py_out = tmp_path / "py_cf.pdb"
    convert_conf_to_pdb(src / "min_conf.dat", py_out)

    c_dir = tmp_path / "c_tc"
    c_dir.mkdir()
    shutil.copy2(c_file, c_dir / "tc.c")
    shutil.copy2(src / "min_conf.dat", c_dir / "min_conf.dat")

    subprocess.run(["g++", "tc.c", "-o", "tc"], cwd=c_dir, check=True)
    subprocess.run(["./tc"], cwd=c_dir, check=True)

    c_out = c_dir / "cf.pdb"
    assert c_out.exists()
    assert py_out.read_text().strip() == c_out.read_text().strip()
