# Installation Guide

This guide covers the installation of DNAfold2 on Linux and macOS systems.

## Prerequisites

### System Requirements

- **Operating System**: Linux (recommended) or macOS
- **CPU**: ≥ 8 cores with 2 threads recommended
- **Memory**: ≥ 8 GB RAM
- **Disk**: ≥ 1 GB free space

### Required Software

| Software | Version | Purpose |
|----------|---------|---------|
| GCC/G++ | ≥ 7.5 | Compile C source code |
| Python | ≥ 3.9 | Python bindings and CLI |
| OpenMP | (included with GCC) | Parallel processing |

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/your-org/DNAfold2.git
cd DNAfold2
```

### 2. Install Python Package

```bash
# Standard installation
pip install -e .

# With development dependencies (for testing)
pip install -e ".[dev]"
```

### 3. Compile C Binaries

```bash
# Build all binaries
make all

# Or build specific components
make core       # TiRNA_sa, TiRNA_optimize
make analysis   # secondary, wham, A_state
make utils      # center, tc, cat, t1
```

### 4. Verify Installation

```bash
# Check CLI is working
dnafold2 --version

# Check installation status
dnafold2 info

# Validate a sequence
dnafold2 validate ATCGATCG
```

## Platform-Specific Notes

### Ubuntu/Debian

```bash
# Install build tools
sudo apt update
sudo apt install build-essential python3-dev

# Verify OpenMP support
apt list --installed | grep libomp
```

### CentOS/RHEL

```bash
# Install build tools
sudo yum groupinstall "Development Tools"
sudo yum install python3-devel

# For newer GCC (if needed)
sudo yum install devtoolset-9
source /opt/rh/devtoolset-9/enable
```

### macOS

```bash
# Install Xcode command line tools
xcode-select --install

# Install GCC via Homebrew (for OpenMP support)
brew install gcc

# Use Homebrew GCC for compilation
export CC=/usr/local/bin/gcc-13
export CXX=/usr/local/bin/g++-13
make all
```

## Compiling Binaries Manually

If the Makefile doesn't work, compile manually:

```bash
cd src/core

# REMC sampler (if source is available)
g++ -O3 -fopenmp TiRNA_remc.c -o ../../bin/TiRNA_remc -lm

# SA sampler
gcc -O3 -Wall TiRNA_sa.c -o ../../bin/TiRNA_sa -lm

# Optimizer
gcc -O3 -Wall -fopenmp TiRNA_optimize.c -o ../../bin/TiRNA_optimize -lm
```

## Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install DNAfold2
pip install -e .
```

## Conda Environment

```bash
# Create conda environment
conda create -n dnafold2 python=3.11 gcc_linux-64 gxx_linux-64
conda activate dnafold2

# Install package
pip install -e .

# Compile binaries
make all
```

## Troubleshooting

### OpenMP Not Found

```bash
# Check OpenMP support
echo | gcc -fopenmp -E -x c - >/dev/null 2>&1 && echo "OK" || echo "Missing"

# On Ubuntu, install:
sudo apt install libomp-dev

# On macOS with Homebrew:
brew install libomp
```

### Python Module Not Found

```bash
# Ensure package is installed
pip show dnafold2

# If using virtual environment, ensure it's activated
which python  # Should point to venv/bin/python
```

### Binaries Not Executable

```bash
# Make binaries executable
chmod +x bin/*
```

## Updating

```bash
# Pull latest changes
git pull origin main

# Reinstall Python package
pip install -e .

# Rebuild binaries (if source changed)
make clean
make all
```
