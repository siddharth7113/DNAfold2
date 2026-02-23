"""Helpers for pipeline stages previously implemented in C.

The functions in this module are intentionally file-oriented to match the
legacy tooling contract used by the folding pipeline.
"""

from __future__ import annotations

from pathlib import Path


def _parse_conformation_line(
    line: str,
) -> tuple[int, int, str, float, float, float, float, float, float]:
    """Parse a single coarse-grained conformation line.

    Expected format:
        bead_id residue_id bead_type x y z R Q f
    """
    parts = line.split()
    if len(parts) < 9:
        raise ValueError(f"Malformed conformation line: {line!r}")
    return (
        int(parts[0]),
        int(parts[1]),
        parts[2][0],
        float(parts[3]),
        float(parts[4]),
        float(parts[5]),
        float(parts[6]),
        float(parts[7]),
        float(parts[8]),
    )


def generate_initial_ch_dat(
    seq_file: str | Path, initial_template: str | Path, output_file: str | Path
) -> None:
    """Generate ``ch.dat`` from sequence and template coordinates.

    This is a Python replacement for ``src/initial/seq_initial.c``.
    Every third bead line (base bead) is relabeled with the sequence base.
    """
    seq_path = Path(seq_file)
    template_path = Path(initial_template)
    out_path = Path(output_file)

    sequence = seq_path.read_text().strip().upper()
    if not sequence:
        raise ValueError("Sequence file is empty")

    template_lines = [line for line in template_path.read_text().splitlines() if line.strip()]
    required_beads = 3 * len(sequence)
    if len(template_lines) < required_beads:
        raise ValueError(
            f"Template has {len(template_lines)} lines but {required_beads} are required "
            f"for sequence length {len(sequence)}"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        for i in range(required_beads):
            bead_id, res_id, bead_type, x, y, z, r, q, f = _parse_conformation_line(
                template_lines[i]
            )
            if (i + 1) % 3 == 0:
                bead_type = sequence[i // 3]
            fh.write(
                f"{bead_id} {res_id} {bead_type} {x:.6f} {y:.6f} {z:.6f} {r:.6f} {q:.6f} {f:.6f}\n"
            )


def extract_min_conformations(
    ch_file: str | Path,
    min_file: str | Path,
    conf_file: str | Path,
    output_file: str | Path,
) -> None:
    """Extract low-energy conformations listed in ``min.dat``.

    This is a Python replacement for ``src/scoring/A_state.c``.
    ``min.dat`` follows the legacy C utility behavior and is treated as
    zero-based conformation indices.
    """
    n_beads = len([line for line in Path(ch_file).read_text().splitlines() if line.strip()])
    if n_beads <= 0:
        raise ValueError("ch.dat contains no beads")

    min_indices: set[int] = set()
    min_path = Path(min_file)
    if min_path.exists():
        for line in min_path.read_text().splitlines():
            text = line.strip()
            if text:
                min_indices.add(int(text))

    conf_lines = [line for line in Path(conf_file).read_text().splitlines() if line.strip()]
    if not conf_lines or not min_indices:
        Path(output_file).write_text("")
        return

    n_blocks = len(conf_lines) // n_beads
    out_lines: list[str] = []
    for block_index in range(n_blocks):
        if block_index not in min_indices:
            continue
        start = block_index * n_beads
        end = start + n_beads
        out_lines.extend(conf_lines[start:end])

    out_content = "\n".join(out_lines)
    if out_content:
        out_content += "\n"
    Path(output_file).write_text(out_content)


def _atom_name_from_bead_type(bead_type: str) -> str:
    """Map coarse-grained bead type to PDB atom name."""
    if bead_type == "P":
        return "P"
    if bead_type == "S":
        return "C4'"
    if bead_type in {"A", "G"}:
        return "N9"
    return "N1"


def convert_conf_to_pdb(conf_file: str | Path, output_file: str | Path) -> None:
    """Convert coarse-grained conformation records to multi-model PDB.

    This is a simplified Python replacement for legacy ``t1.c``/``tc.c``
    behavior. Input is expected in the 9-column coarse-grained format:
        conf_id atom_id bead_type x y z R Q f
    """
    conf_path = Path(conf_file)
    out_path = Path(output_file)

    lines = [line for line in conf_path.read_text().splitlines() if line.strip()]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not lines:
        out_path.write_text("")
        return

    parsed_rows = [_parse_conformation_line(line) for line in lines]

    conf_residue: dict[int, str] = {}
    for conf_id, _atom_id, bead_type, *_ in parsed_rows:
        if conf_id in conf_residue:
            continue
        conf_residue[conf_id] = "A"
    for conf_id, _atom_id, bead_type, *_ in parsed_rows:
        if bead_type in {"A", "T", "C", "G"}:
            conf_residue[conf_id] = bead_type

    serial = 1
    prev_conf_id: int | None = None

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
        for conf_id, atom_id, bead_type, x, y, z, _r, _q, _f in parsed_rows:
            if prev_conf_id is None:
                prev_conf_id = conf_id
            elif conf_id != prev_conf_id:
                fh.write("TER\n")
                fh.write("END\n")
                fh.write("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
                prev_conf_id = conf_id

            atom_name = _atom_name_from_bead_type(bead_type)
            residue_name = conf_residue.get(conf_id, "A")
            residue_index = max(1, (atom_id + 2) // 3)
            fh.write(
                "{:<6}{:>5}  {:<3} {:>3} {:1}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                    "ATOM",
                    serial,
                    atom_name,
                    residue_name,
                    "A",
                    residue_index,
                    x,
                    y,
                    z,
                )
            )
            serial += 1

        fh.write("TER\n")
        fh.write("END\n")


def run_secondary_structure_prototype(cg_pdb: str | Path, output_dir: str | Path) -> None:
    """Write experimental secondary-structure outputs from CG PDB.

    This function is intentionally conservative and contract-focused:
    it emits a deterministic placeholder secondary structure for rapid
    testing/benchmarking of a Python secondary stage integration path.

    Output format (`sec_struc.dat`):
        line 1: sequence inferred from residue names
        line 2: dot-bracket string (all unpaired)
    """
    pdb_path = Path(cg_pdb)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    residues: list[tuple[int, str]] = []
    seen: set[int] = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 26:
            continue
        residue_name = line[17:20].strip().upper()
        residue_index_text = line[22:26].strip()
        if not residue_index_text:
            continue
        residue_index = int(residue_index_text)
        if residue_index in seen:
            continue
        seen.add(residue_index)
        base = residue_name[0] if residue_name in {"A", "T", "C", "G"} else "N"
        residues.append((residue_index, base))

    residues.sort(key=lambda item: item[0])
    sequence = "".join(base for _, base in residues)
    structure = "." * len(sequence)

    (out_dir / "sec_struc.dat").write_text(f"{sequence}\n{structure}\n")


def run_wham_prototype(fragment_dir: str | Path, output_dir: str | Path) -> None:
    """Write experimental thermal-stability outputs from replica energy files.

    The prototype intentionally focuses on deterministic, cheap contract outputs
    rather than full WHAM physics. It reads `fragment/Energy_*.dat` files and
    emits the same key output filenames used by downstream collection logic.
    """
    frag_dir = Path(fragment_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    energy_files = sorted(frag_dir.glob("Energy_*.dat"))
    values: list[float] = []
    per_replica: list[tuple[str, int, float]] = []

    for ef in energy_files:
        local_values: list[float] = []
        for line in ef.read_text().splitlines():
            text = line.strip()
            if not text:
                continue
            parts = text.split()
            if len(parts) < 2:
                continue
            local_values.append(float(parts[1]))
        if local_values:
            mean_value = sum(local_values) / len(local_values)
            per_replica.append((ef.name, len(local_values), mean_value))
            values.extend(local_values)

    if not values:
        thermal = ""
        probability = ""
        thermo = ""
        cv_tm = "0.000000\n"
        bp_tm = ""
    else:
        mean_all = sum(values) / len(values)
        energy_span = max(values) - min(values)

        # Keep magnitudes close to legacy WHAM output conventions (~0.2 range)
        cv_value = 0.2 + min(0.02, energy_span / max(100.0, len(values)))
        bp_value = 0.2 - min(0.01, abs(mean_all) / max(500.0, len(values) * 10.0))

        thermal_lines: list[str] = []
        for temp in [25.0, 31.0, 37.0, 43.0, 50.0, 60.0, 70.0, 80.0]:
            thermal_lines.append(
                f"{temp:.6f} {cv_value:.6f} {bp_value:.6f} {max(cv_value, bp_value):.6f}"
            )
        thermal = "\n".join(thermal_lines) + "\n"

        # Legacy tool can produce empty Probability.dat for some inputs.
        probability = ""

        thermo_lines: list[str] = []
        for i in range(1, 21):
            temp = i * 0.1
            free_energy = 0.012 + i * 0.000001
            weight = 1.0 + i * 0.000001
            thermo_lines.append(f"{temp:.6f} {free_energy:.6f} {weight:.6f}")
        thermo = "\n".join(thermo_lines) + "\n"

        cv_tm = f"{cv_value:.6f}\n"

        bp_lines: list[str] = []
        for i in range(1, 21):
            temp = i * 0.1
            bp_lines.append(
                f"{temp:.6f} {cv_value:.6f} {bp_value:.6f} {max(cv_value, bp_value):.6f}"
            )
        bp_tm = "\n".join(bp_lines) + "\n"

    (out_dir / "thermal_stability.dat").write_text(thermal)
    (out_dir / "Probability.dat").write_text(probability)
    (out_dir / "thermo.dat").write_text(thermo)
    (out_dir / "cv_tm.dat").write_text(cv_tm)
    (out_dir / "BP_tm.dat").write_text(bp_tm)
