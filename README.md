# DNAfold2

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An improved coarse-grained (CG) model for **ab initio prediction of complex DNA structure folding**.

## Features

- 🧬 **3D Structure Prediction** - Predict DNA tertiary structures from sequence
- 🔬 **Coarse-Grained Model** - Three-bead CG model with sequence-dependent potentials
- 🌡️ **Thermal Stability** - Compute melting temperatures and stability profiles
- ⚡ **Multiple Algorithms** - Replica Exchange Monte Carlo (REMC) and Simulated Annealing (SA)
- 🐍 **Python API** - Modern Python bindings for programmatic access
- 🖥️ **Cluster Ready** - SLURM templates for HPC environments

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/DNAfold2.git
cd DNAfold2

# Install Python package
pip install -e .

# Compile C binaries (requires GCC with OpenMP)
make all
```

### Basic Usage

**Command Line:**
```bash
# Fold a DNA sequence
dnafold2 fold --sequence ATCCTAGTTATAGGAT --output results/

# From a sequence file
dnafold2 fold --input seq.dat --config config.dat --output results/

# Validate a sequence
dnafold2 validate ATCCTAGTTATAGGAT
```

**Python API:**
```python
from dnafold2 import DNAFolder, FoldingConfig

# Configure the simulation
config = FoldingConfig(
    sampling_method="remc",      # or "sa" for Simulated Annealing
    folding_steps=500000,        # REMC/SA simulation steps
    optimizing_steps=100000,     # Energy optimization steps
    na_concentration=1000,       # Na+ concentration (mM)
    mg_concentration=0,          # Mg2+ concentration (mM)
    n_structures=10              # Number of structures to predict
)

# Fold a sequence
folder = DNAFolder(config)
result = folder.fold("ATCCTAGTTATAGGAT", output_dir="results/")

# Access results
print(f"CG structures: {result.cg_structures}")
print(f"Elapsed time: {result.elapsed_time:.1f}s")
```

## Requirements

| Component | Requirement |
|-----------|-------------|
| OS | Linux or macOS |
| CPU | ≥ 8 cores recommended |
| Compiler | GCC ≥ 7.5 (with C++11 and OpenMP support) |
| Python | ≥ 3.9 |
| Dependencies | Biopython ≥ 1.81 |

## Project Structure

```
DNAfold2/
├── dnafold2/           # Python package
│   ├── core.py         # DNAFolder class
│   ├── config.py       # Configuration handling
│   ├── cli.py          # Command-line interface
│   └── utils.py        # Utility functions
├── src/                # C source code
│   ├── core/           # Main algorithms (REMC, SA, optimization)
│   ├── analysis/       # Secondary structure, WHAM analysis
│   ├── rebuild/        # All-atom reconstruction
│   └── utils/          # Utility programs
├── bin/                # Compiled binaries
├── examples/           # Example inputs
├── scripts/slurm/      # SLURM job templates
├── tests/              # Unit tests
└── docs/               # Documentation
```

## Configuration

Create a `config.dat` file or use `FoldingConfig` in Python:

```
Sampling 1                 # 1=REMC, 2=SA
Folding_steps 500000       # Simulation steps (min 500000 for structure)
Optimizing_steps 100000    # Optimization steps (min 100000)
CNa 1000                   # Na+ concentration (mM)
CMg 0                      # Mg2+ concentration (mM)
Ncout 10                   # Number of output structures
```

### Recommended Settings

| Use Case | Folding Steps | Optimizing Steps |
|----------|---------------|------------------|
| Quick test | 150,000 | 50,000 |
| 3D structure | ≥ 500,000 | ≥ 100,000 |
| Thermal stability | ≥ 2,000,000 | ≥ 100,000 |
| Publication quality | ≥ 4,000,000 | ≥ 500,000 |

## Output Files

Results are organized in the output directory:

| Directory | Contents |
|-----------|----------|
| `CG_structure/` | Coarse-grained PDB structures |
| `All_atom_structure/` | All-atom PDB models |
| `Secondary_structure/` | Dot-bracket notation files |
| `Folding_trajectory/` | Simulation trajectories |
| `Thermal_Stability/` | Stability data and melting curves |

## Cluster Usage (SLURM)

Submit jobs on HPC clusters:

```bash
# Single sequence
sbatch scripts/slurm/submit_folding.sh

# Batch processing (array job)
sbatch scripts/slurm/submit_batch.sh
```

Edit the SLURM scripts to configure:
- Sequence input
- Resource allocation (CPUs, memory, time)
- Email notifications

## Compiling from Source

If you need to recompile the C binaries:

```bash
# Build all binaries
make all

# Or build specific components
make core       # REMC/SA algorithms
make analysis   # Analysis tools
make utils      # Utility programs

# Clean build
make clean
```

### Compiler Requirements

```bash
# Check GCC version (needs ≥ 7.5)
gcc --version

# Check OpenMP support
echo | gcc -fopenmp -E -x c - >/dev/null 2>&1 && echo "OpenMP OK" || echo "OpenMP missing"
```

## Running Tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=dnafold2
```

## References

1. Z.-C. Mu, Y.-L. Tan, B.-G. Zhang, J. Liu, Y.-Z. Shi. **Ab initio predictions for 3D structure and stability of single- and double-stranded DNAs in ion solutions.** *PLoS Comput. Biol.* 18, e1010501 (2022).

2. X. Wang, Y.-Z. Shi. **3D structure and stability prediction of DNA with multi-way junctions in ionic solutions.** arXiv:2501.11891

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Authors

- **Xunxun Wang** - Guizhou Medical University
- **Ya-Zhou Shi** - Wuhan Textile University

Developed at Tan's Group (Wuhan University) and collaborating institutions.
