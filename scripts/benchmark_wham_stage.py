"""Benchmark experimental Python WHAM prototype against legacy C wham tool."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from dnafold2.stage_tools import run_wham_prototype


def write_synthetic_fragment(fragment_dir: Path, n_replicas: int, n_frames: int) -> None:
    """Write synthetic Energy_* and bp_* files for benchmark runs."""
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


def run_python(fragment_dir: Path, out_dir: Path) -> float:
    start = time.perf_counter()
    run_wham_prototype(fragment_dir, out_dir)
    return time.perf_counter() - start


def run_c(analysis_src: Path, out_dir: Path, timeout: int) -> tuple[int, float, str]:
    shutil.copy2(analysis_src / "wham.c", out_dir / "wham.c")
    start = time.perf_counter()
    try:
        build = subprocess.run(
            ["gcc", "-Wall", "wham.c", "-o", "wham", "-lm"],
            cwd=out_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if build.returncode != 0:
            return build.returncode, time.perf_counter() - start, build.stderr

        run = subprocess.run(
            ["./wham"],
            cwd=out_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return run.returncode, time.perf_counter() - start, run.stderr
    except subprocess.TimeoutExpired:
        return 124, time.perf_counter() - start, f"Timed out after {timeout}s"


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Python WHAM prototype")
    parser.add_argument("--threshold", type=float, default=0.15)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--frames", type=int, default=50)
    parser.add_argument("--workdir", type=Path, default=Path("."))
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    analysis_src = workdir / "src" / "analysis"

    with (
        tempfile.TemporaryDirectory(prefix="bench_wham_py_") as py_tmp,
        tempfile.TemporaryDirectory(prefix="bench_wham_c_") as c_tmp,
    ):
        py_dir = Path(py_tmp)
        c_dir = Path(c_tmp)
        write_synthetic_fragment(py_dir / "fragment", n_replicas=10, n_frames=args.frames)
        write_synthetic_fragment(c_dir / "fragment", n_replicas=10, n_frames=args.frames)

        py_time = run_python(py_dir / "fragment", py_dir)
        c_rc, c_time, c_err = run_c(analysis_src, c_dir, args.timeout)

    print(f"Python WHAM prototype: time={py_time:.6f}s")
    print(f"C wham tool: rc={c_rc}, time={c_time:.6f}s")

    if c_rc != 0:
        print("C wham run failed:")
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
