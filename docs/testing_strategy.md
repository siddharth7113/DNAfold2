# Testing Strategy for Fast Iteration

This project has expensive end-to-end stages, so tests are split into tiers.

## Test tiers

- `unit`: fast deterministic tests for pure Python logic.
- `contract`: fixture-driven stage tests that validate file contracts and output shape.
- `smoke`: short pipeline checks (opt-in with `--run-smoke`).
- `slow`: long-running integration/scientific tests (opt-in with `--run-slow`).

## Commands

```bash
# Fast default tiers
pytest -m "unit or contract" -v

# Include smoke checks
pytest --run-smoke -m "unit or contract or smoke" -v

# Run slow tests only
pytest --run-slow -m "slow" -v
```

Quick CLI smoke profile:

```bash
dnafold2 fold --sequence ATCCTAGTTATAGGAT --output results_quick --quick-test --stop-after scoring
```

Benchmark Python vs C stage tools (15% threshold example):

```bash
python scripts/benchmark_stage_tools.py --stop-after initial --threshold 0.15 --timeout 180
```

## Regression policy

- Deterministic stages: exact output comparisons are preferred.
- Stochastic stages: use invariant checks and tolerances, not bitwise equality.
- Full long runs should be reserved for nightly/manual validation.
