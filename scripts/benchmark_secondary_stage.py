"""Benchmark experimental Python secondary prototype against legacy C tool."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from dnafold2.stage_tools import run_secondary_structure_prototype


def run_python(cg_pdb: Path, out_dir: Path) -> float:
    start = time.perf_counter()
    run_secondary_structure_prototype(cg_pdb, out_dir)
    return time.perf_counter() - start


def run_c(cg_pdb: Path, rebuild_src: Path, out_dir: Path, timeout: int) -> tuple[int, float, str]:
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
            shutil.copy2(src, out_dir / name)
    shutil.copy2(cg_pdb, out_dir / "CG.pdb")

    start = time.perf_counter()
    try:
        build = subprocess.run(
            ["gcc", "-Wall", "secondary.c", "-o", "secondary", "-lm"],
            cwd=out_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if build.returncode != 0:
            return build.returncode, time.perf_counter() - start, build.stderr

        run = subprocess.run(
            ["./secondary"],
            cwd=out_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return run.returncode, time.perf_counter() - start, run.stderr
    except subprocess.TimeoutExpired:
        return 124, time.perf_counter() - start, f"Timed out after {timeout}s"


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Python secondary prototype")
    parser.add_argument("--threshold", type=float, default=0.15)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--cg-pdb",
        type=Path,
        default=Path("tests/fixtures/stage_tools/secondary_prototype/CG.pdb"),
        help="Input CG.pdb fixture",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("."),
        help="Repository root",
    )
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    cg_pdb = (workdir / args.cg_pdb).resolve()
    rebuild_src = workdir / "src" / "rebuild"

    with (
        tempfile.TemporaryDirectory(prefix="bench_secondary_py_") as py_tmp,
        tempfile.TemporaryDirectory(prefix="bench_secondary_c_") as c_tmp,
    ):
        py_time = run_python(cg_pdb, Path(py_tmp))
        c_rc, c_time, c_err = run_c(cg_pdb, rebuild_src, Path(c_tmp), args.timeout)

    print(f"Python secondary prototype: time={py_time:.6f}s")
    print(f"C secondary tool: rc={c_rc}, time={c_time:.6f}s")

    if c_rc != 0:
        print("C secondary run failed:")
        print(c_err[-1000:])
        return c_rc

    slowdown = (py_time - c_time) / c_time if c_time > 0 else 0.0
    print(f"Slowdown ratio: {slowdown:.3%}")
    if slowdown > args.threshold:
        print(f"FAIL: slowdown {slowdown:.3%} exceeds threshold {args.threshold:.3%}")
        return 3

    print(f"PASS: slowdown {slowdown:.3%} within threshold {args.threshold:.3%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
