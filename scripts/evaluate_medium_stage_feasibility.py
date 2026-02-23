"""Evaluate whether medium stages are ready to migrate from C to Python.

This script runs a bounded comparison for secondary and WHAM stages:
1) Runtime comparison (Python prototype vs legacy C)
2) Output parity check on representative fixtures

Current expectation: prototypes are useful for fast testing/contract stability,
but C implementations remain the scientific default until parity is achieved.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from dnafold2.stage_tools import run_secondary_structure_prototype, run_wham_prototype


def _run(cmd: list[str], cwd: Path, timeout: int) -> tuple[int, str, str]:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", f"Timed out after {timeout}s"


def _write_synthetic_wham_fragment(
    fragment_dir: Path, n_replicas: int = 10, n_frames: int = 50
) -> None:
    fragment_dir.mkdir(parents=True, exist_ok=True)
    for i in range(n_replicas):
        energy_lines: list[str] = []
        bp_lines: list[str] = []
        for frame in range(1, n_frames + 1):
            energy = -20.0 - i * 0.5 - frame * 0.01
            energy_lines.append(f"{frame} {energy:.6f} -1.0 1.0 -2.0 -0.5 0.0 0.0")
            bp = frame % 5
            bp_lines.append(f"{frame} {frame} {bp} 0.0 0.0 {bp} 0.0 0.0")
        (fragment_dir / f"Energy_{i}.dat").write_text("\n".join(energy_lines) + "\n")
        (fragment_dir / f"bp_{i}.dat").write_text("\n".join(bp_lines) + "\n")


def evaluate_secondary(repo_root: Path, timeout: int) -> tuple[float, bool, str]:
    fixture_cg = repo_root / "tests" / "fixtures" / "stage_tools" / "secondary_prototype" / "CG.pdb"
    rebuild_src = repo_root / "src" / "rebuild"

    with (
        tempfile.TemporaryDirectory(prefix="eval_secondary_py_") as py_tmp,
        tempfile.TemporaryDirectory(prefix="eval_secondary_c_") as c_tmp,
    ):
        py_dir = Path(py_tmp)
        c_dir = Path(c_tmp)

        # Python prototype
        py_start = time.perf_counter()
        run_secondary_structure_prototype(fixture_cg, py_dir)
        py_elapsed = time.perf_counter() - py_start

        # C baseline
        for name in [
            "secondary.c",
            "stem.dat",
            "stem_kissing.dat",
            "state.dat",
            "DNA_type",
            "RNA_type",
            "cs.dat",
        ]:
            src = rebuild_src / name
            if src.exists():
                shutil.copy2(src, c_dir / name)
        shutil.copy2(fixture_cg, c_dir / "CG.pdb")

        rc_build, _, err_build = _run(
            ["gcc", "-Wall", "secondary.c", "-o", "secondary", "-lm"], c_dir, timeout
        )
        if rc_build != 0:
            return 1.0, False, f"secondary.c compile failed: {err_build[-300:]}"

        c_start = time.perf_counter()
        rc_run, _, err_run = _run(["./secondary"], c_dir, timeout)
        c_elapsed = time.perf_counter() - c_start
        if rc_run != 0:
            return 1.0, False, f"secondary run failed: {err_run[-300:]}"

        py_out = (
            (py_dir / "sec_struc.dat").read_text() if (py_dir / "sec_struc.dat").exists() else ""
        )
        c_out = (c_dir / "sec_struc.dat").read_text() if (c_dir / "sec_struc.dat").exists() else ""
        parity = bool(py_out) and py_out == c_out

        slowdown = (py_elapsed - c_elapsed) / c_elapsed if c_elapsed > 0 else 0.0
        reason = "parity match" if parity else "parity mismatch"
        return slowdown, parity, reason


def evaluate_wham(repo_root: Path, timeout: int) -> tuple[float, bool, str]:
    analysis_src = repo_root / "src" / "analysis"

    with (
        tempfile.TemporaryDirectory(prefix="eval_wham_py_") as py_tmp,
        tempfile.TemporaryDirectory(prefix="eval_wham_c_") as c_tmp,
    ):
        py_dir = Path(py_tmp)
        c_dir = Path(c_tmp)

        _write_synthetic_wham_fragment(py_dir / "fragment")
        _write_synthetic_wham_fragment(c_dir / "fragment")

        py_start = time.perf_counter()
        run_wham_prototype(py_dir / "fragment", py_dir)
        py_elapsed = time.perf_counter() - py_start

        shutil.copy2(analysis_src / "wham.c", c_dir / "wham.c")
        rc_build, _, err_build = _run(
            ["gcc", "-Wall", "wham.c", "-o", "wham", "-lm"], c_dir, timeout
        )
        if rc_build != 0:
            return 1.0, False, f"wham.c compile failed: {err_build[-300:]}"

        c_start = time.perf_counter()
        rc_run, _, err_run = _run(["./wham"], c_dir, timeout)
        c_elapsed = time.perf_counter() - c_start
        if rc_run != 0:
            return 1.0, False, f"wham run failed: {err_run[-300:]}"

        # Parity criterion: all output files exist and thermal_stability exact match.
        py_ts = py_dir / "thermal_stability.dat"
        c_ts = c_dir / "thermal_stability.dat"
        parity = py_ts.exists() and c_ts.exists() and py_ts.read_text() == c_ts.read_text()

        slowdown = (py_elapsed - c_elapsed) / c_elapsed if c_elapsed > 0 else 0.0
        reason = "parity match" if parity else "parity mismatch"
        return slowdown, parity, reason


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate medium-stage migration feasibility")
    parser.add_argument("--threshold", type=float, default=0.15, help="Allowed slowdown threshold")
    parser.add_argument("--timeout", type=int, default=120, help="Timeout per C compile/run")
    parser.add_argument("--workdir", type=Path, default=Path("."), help="Repository root")
    parser.add_argument(
        "--fail-on-nonready",
        action="store_true",
        help="Return non-zero if stages are not migration-ready",
    )
    args = parser.parse_args()

    repo_root = args.workdir.resolve()

    sec_slow, sec_parity, sec_reason = evaluate_secondary(repo_root, args.timeout)
    wham_slow, wham_parity, wham_reason = evaluate_wham(repo_root, args.timeout)

    sec_speed_ok = sec_slow <= args.threshold
    wham_speed_ok = wham_slow <= args.threshold

    sec_ready = sec_speed_ok and sec_parity
    wham_ready = wham_speed_ok and wham_parity

    print("Secondary stage:")
    print(
        f"  slowdown={sec_slow:.3%} | speed_ok={sec_speed_ok} | parity={sec_parity} ({sec_reason})"
    )
    print("WHAM stage:")
    print(
        f"  slowdown={wham_slow:.3%} | speed_ok={wham_speed_ok} | parity={wham_parity} ({wham_reason})"
    )

    if sec_ready and wham_ready:
        print("Recommendation: medium stages are ready for Python promotion.")
        return 0

    print("Recommendation: keep C as default for medium stages; continue parity work.")
    if args.fail_on_nonready:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
