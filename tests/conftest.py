"""Pytest configuration for tiered test execution."""

from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add CLI flags for optional slow/smoke test tiers."""
    parser.addoption(
        "--run-smoke",
        action="store_true",
        default=False,
        help="Run tests marked as smoke",
    )
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run tests marked as slow",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip optional test tiers unless explicitly requested."""
    run_smoke = config.getoption("--run-smoke")
    run_slow = config.getoption("--run-slow")

    skip_smoke = pytest.mark.skip(reason="need --run-smoke option to run")
    skip_slow = pytest.mark.skip(reason="need --run-slow option to run")

    for item in items:
        if "smoke" in item.keywords and not run_smoke:
            item.add_marker(skip_smoke)
        if "slow" in item.keywords and not run_slow:
            item.add_marker(skip_slow)
