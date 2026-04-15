"""Benchmark Python stage tools against legacy C tools.

This script runs two quick pipeline invocations up to a chosen stage and
compares wall-clock runtime:
1) Python stage tools (default path)
2) Legacy C stage tools (`--use-c-stage-tools`)

The benchmark is intended for developer validation and does not replace full
scientific integration tests.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path


def run_case(command: list[str], cwd: Path, timeout: int) -> tuple[int, float, str]:
    """Run one benchmark case and return (returncode, seconds, stderr)."""
    start = time.perf_counter()
    try:
        proc = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        elapsed = time.perf_counter() - start
        return 124, elapsed, f"Timed out after {timeout}s"
    elapsed = time.perf_counter() - start
    return proc.returncode, elapsed, proc.stderr


def build_base_command(output_dir: Path, stop_after: str, sequence: str) -> list[str]:
    """Build base CLI command for quick benchmark run."""
    return [
        "dnafold2",
        "fold",
        "--sequence",
        sequence,
        "--output",
        str(output_dir),
        "--quick-test",
        "--stop-after",
        stop_after,
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Python vs C stage tools")
    parser.add_argument(
        "--sequence",
        default="ATCCTAGTTATAGGAT",
        help="DNA sequence for benchmark run",
    )
    parser.add_argument(
        "--stop-after",
        choices=["initial", "folding", "scoring"],
        default="initial",
        help="Pipeline stage at which to stop benchmark runs",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.15,
        help="Maximum allowed slowdown ratio for Python tools (default: 0.15 = 15%%)",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("."),
        help="Repository working directory",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=180,
        help="Timeout per benchmark run in seconds",
    )
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    out_py = workdir / "benchmark_py_tools"
    out_c = workdir / "benchmark_c_tools"

    for out in (out_py, out_c):
        if out.exists():
            shutil.rmtree(out)

    cmd_py = build_base_command(out_py, args.stop_after, args.sequence)
    cmd_c = build_base_command(out_c, args.stop_after, args.sequence) + ["--use-c-stage-tools"]

    rc_py, t_py, err_py = run_case(cmd_py, workdir, args.timeout)
    rc_c, t_c, err_c = run_case(cmd_c, workdir, args.timeout)

    print(f"Python stage tools: rc={rc_py}, time={t_py:.3f}s")
    print(f"C stage tools:      rc={rc_c}, time={t_c:.3f}s")

    if rc_py != 0:
        print("Python stage tool run failed:")
        print(err_py[-1000:])
        return rc_py
    if rc_c != 0:
        print("C stage tool run failed:")
        print(err_c[-1000:])
        return rc_c

    if t_c <= 0:
        print("Invalid C baseline runtime")
        return 2

    slowdown = (t_py - t_c) / t_c
    print(f"Slowdown ratio: {slowdown:.3%}")

    if slowdown > args.threshold:
        print(f"FAIL: slowdown {slowdown:.3%} exceeds threshold {args.threshold:.3%}")
        return 3

    print(f"PASS: slowdown {slowdown:.3%} within threshold {args.threshold:.3%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
