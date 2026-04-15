"""Deterministic fixture-based contract tests for stage tools."""

from __future__ import annotations

from pathlib import Path

import pytest

from dnafold2.stage_tools import (
    convert_conf_to_pdb,
    extract_min_conformations,
    generate_initial_ch_dat,
)


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "stage_tools"


@pytest.mark.contract
def test_seq_initial_contract_fixture(tmp_path: Path) -> None:
    src = FIXTURE_ROOT / "seq_initial"
    out_file = tmp_path / "ch.dat"

    generate_initial_ch_dat(src / "seq.dat", src / "initial.dat", out_file)

    expected = (src / "expected_ch.dat").read_text()
    assert out_file.read_text() == expected


@pytest.mark.contract
def test_min_extract_contract_fixture(tmp_path: Path) -> None:
    src = FIXTURE_ROOT / "min_extract"
    out_file = tmp_path / "min_conf.dat"

    extract_min_conformations(src / "ch.dat", src / "min.dat", src / "conf_0.dat", out_file)

    expected = (src / "expected_min_conf.dat").read_text()
    assert out_file.read_text() == expected


@pytest.mark.contract
def test_pdb_convert_contract_fixture(tmp_path: Path) -> None:
    src = FIXTURE_ROOT / "pdb_convert"
    out_file = tmp_path / "cf.pdb"

    convert_conf_to_pdb(src / "min_conf.dat", out_file)

    expected = (src / "expected_cf.pdb").read_text()
    assert out_file.read_text() == expected
