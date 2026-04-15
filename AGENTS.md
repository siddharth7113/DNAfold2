# AGENTS.md
Guidance for coding agents working in `DNAfold2`.

Repository map:
- `dnafold2/`: Python package (CLI + orchestration)
- `src/`: C/C++ scientific kernels and support files
- `tests/`: pytest tests
- `Makefile`: binary build and housekeeping targets
- `pyproject.toml`: Python packaging + tool config

Run all commands from `/home/siddharth/work/DNAfold2`.

## Build, lint, test commands

### Setup
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
pip install -e ".[dev]"
```

### Build
Build all binaries:
```bash
make all
```
Build subsets:
```bash
make core
make analysis
make utils
make rebuild
make initial
```
Install helpers:
```bash
make install
make install-dev
```
Clean:
```bash
make clean
make distclean
```

### Tests (pytest)
Full suite:
```bash
pytest
```
Common explicit form:
```bash
pytest tests/ -v
make test
```
Coverage:
```bash
pytest tests/ --cov=dnafold2
```
Single test file:
```bash
pytest tests/test_config.py -v
```
Single class:
```bash
pytest tests/test_config.py::TestFoldingConfig -v
```
Single test function:
```bash
pytest tests/test_config.py::TestFoldingConfig::test_default_config -v
```
Keyword expression:
```bash
pytest -k "validate and not integration" -v
```
Debug with live output:
```bash
pytest -s tests/test_cli.py::TestValidateCommand::test_invalid_sequence -v
```

Tiered test execution:
```bash
# Fast default tiers
pytest -m "unit or contract" -v

# Include smoke tests
pytest --run-smoke -m "unit or contract or smoke" -v

# Include slow tests (opt-in)
pytest --run-slow -m "slow" -v
```

Quick pipeline smoke invocation (CLI profile):
```bash
dnafold2 fold --sequence ATCCTAGTTATAGGAT --output results_quick --quick-test --stop-after scoring
```

Stage tool benchmark (Python vs legacy C, 15% threshold):
```bash
python scripts/benchmark_stage_tools.py --stop-after initial --threshold 0.15 --timeout 180
```

### Lint, format, type-check
Black (`line-length = 100`):
```bash
black dnafold2 tests
black --check dnafold2 tests
```
Ruff:
```bash
ruff check dnafold2 tests
ruff check --fix dnafold2 tests
```
Mypy:
```bash
mypy dnafold2
```
Recommended pre-PR run:
```bash
ruff check dnafold2 tests && black --check dnafold2 tests && mypy dnafold2 && pytest -v
```

## Code style guidelines

### Source of truth
- Tool behavior comes from `pyproject.toml`.
- Follow established patterns in nearby files unless task requirements say otherwise.

### Imports
- Order imports: standard library, third-party, local package.
- Keep imports explicit; avoid wildcard imports.
- Remove unused imports (keep Ruff clean).
- Use relative imports for internal modules (example: `from .config import FoldingConfig`).

### Formatting
- Use Black-compatible formatting.
- Keep max line length at 100.
- Write docstrings for public modules, classes, and functions.
- Prefer extracting helpers instead of adding deeply nested logic.

### Types
- Add type hints for new public APIs.
- Add return annotations for non-trivial internal helpers.
- Prefer `pathlib.Path` for filesystem paths.
- Use existing union style (`str | Path`) where appropriate.
- Use dataclasses for structured config/result data.
- Keep code mypy-friendly; avoid expanding untyped public surfaces.

### Naming conventions
- Classes: `PascalCase`.
- Functions/variables/methods: `snake_case`.
- Constants: `UPPER_SNAKE_CASE`.
- Tests: `test_*.py`, `Test*`, and `test_*` names.
- CLI handlers: `cmd_<subcommand>` naming.

### Error handling
- Validate input early and fail fast.
- Raise specific exceptions (`ValueError`, `FileNotFoundError`, etc.).
- Catch broad exceptions only at CLI boundaries where exit codes are returned.
- Avoid silent failures in library code.
- Include actionable context in errors (step, file, command).

### Logging
- Use module logger: `logger = logging.getLogger(__name__)`.
- Use lazy formatting in logs (`logger.info("x=%s", x)`).
- Log major pipeline stages at `info`, detailed diagnostics at `debug`.

### Subprocess and filesystem patterns
- Prefer `subprocess.run(..., capture_output=True, text=True, timeout=...)`.
- Check `returncode` and report stderr/stdout on failure.
- Use `Path` APIs (`exists`, `mkdir`, `glob`, `read_text`, `write_text`).
- Ensure executable permissions for binaries that must run (`chmod 0o755`).

### Testing expectations
- Add/update tests for any behavior change.
- Keep unit tests deterministic and fast.
- Mark/skip integration-heavy flows that require compiled binaries or long runs.
- Add regression tests for bug fixes.

### Scope boundaries
- Prefer changes in `dnafold2/` for most tasks.
- Touch legacy `src/` code only when required.
- Preserve CLI argument names, defaults, and exit code behavior.

## Cursor/Copilot instructions
No additional instruction files were found:
- `.cursorrules`
- `.cursor/rules/`
- `.github/copilot-instructions.md`

If these files appear later, treat them as higher-priority instructions and update this document.
