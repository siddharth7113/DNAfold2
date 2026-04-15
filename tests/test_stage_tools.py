"""Tests for Python stage helper replacements."""

from __future__ import annotations

from pathlib import Path

import pytest

from dnafold2.stage_tools import (
    convert_conf_to_pdb,
    extract_min_conformations,
    generate_initial_ch_dat,
)


@pytest.mark.unit
def test_generate_initial_ch_dat_relabels_base_beads(tmp_path: Path) -> None:
    seq_file = tmp_path / "seq.dat"
    template_file = tmp_path / "initial.dat"
    output_file = tmp_path / "ch.dat"

    seq_file.write_text("AT\n")
    template_file.write_text(
        "\n".join(
            [
                "1 1 P 0 0 0 1 2 3",
                "2 1 S 0 0 0 1 2 3",
                "3 1 X 0 0 0 1 2 3",
                "4 2 P 0 0 0 1 2 3",
                "5 2 S 0 0 0 1 2 3",
                "6 2 X 0 0 0 1 2 3",
            ]
        )
        + "\n"
    )

    generate_initial_ch_dat(seq_file, template_file, output_file)

    lines = output_file.read_text().strip().splitlines()
    assert len(lines) == 6
    assert lines[2].split()[2] == "A"
    assert lines[5].split()[2] == "T"


@pytest.mark.contract
def test_extract_min_conformations_selects_requested_blocks(tmp_path: Path) -> None:
    ch_file = tmp_path / "ch.dat"
    min_file = tmp_path / "min.dat"
    conf_file = tmp_path / "conf_0.dat"
    output_file = tmp_path / "min_conf.dat"

    ch_file.write_text("1 1 P 0 0 0 1 2 3\n2 1 S 0 0 0 1 2 3\n3 1 A 0 0 0 1 2 3\n")
    min_file.write_text("2\n")

    conf_file.write_text(
        "\n".join(
            [
                # block 1
                "1 1 P 0 0 0 1 2 3",
                "2 1 S 0 0 0 1 2 3",
                "3 1 A 0 0 0 1 2 3",
                # block 2
                "4 2 P 1 1 1 1 2 3",
                "5 2 S 1 1 1 1 2 3",
                "6 2 T 1 1 1 1 2 3",
                # block 3
                "7 3 P 2 2 2 1 2 3",
                "8 3 S 2 2 2 1 2 3",
                "9 3 C 2 2 2 1 2 3",
            ]
        )
        + "\n"
    )

    extract_min_conformations(ch_file, min_file, conf_file, output_file)

    selected = output_file.read_text().strip().splitlines()
    assert len(selected) == 3
    assert selected[0].startswith("4 2 P")
    assert selected[2].startswith("6 2 T")


@pytest.mark.contract
def test_convert_conf_to_pdb_writes_multimodel_output(tmp_path: Path) -> None:
    conf_file = tmp_path / "min_conf.dat"
    out_file = tmp_path / "cf.pdb"

    conf_file.write_text(
        "\n".join(
            [
                "1 1 P 0 0 0 1 2 3",
                "1 2 S 0 1 0 1 2 3",
                "1 3 A 0 2 0 1 2 3",
                "2 4 P 1 0 0 1 2 3",
                "2 5 S 1 1 0 1 2 3",
                "2 6 T 1 2 0 1 2 3",
            ]
        )
        + "\n"
    )

    convert_conf_to_pdb(conf_file, out_file)

    content = out_file.read_text()
    assert "CRYST1" in content
    assert content.count("END") == 2
    assert "ATOM" in content
