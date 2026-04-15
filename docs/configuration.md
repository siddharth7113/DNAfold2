# Configuration Guide

DNAfold2 uses a configuration file (`config.dat`) to control simulation parameters. This guide explains all available options.

## Configuration File Format

The configuration file uses a simple key-value format:

```
Sampling 1              # 1=REMC, 2=SA
Folding_steps 500000    # Number of simulation steps
Optimizing_steps 100000 # Number of optimization steps
CNa 1000                # Na+ concentration (mM)
CMg 0                   # Mg2+ concentration (mM)
Ncout 10                # Number of output structures
```

Lines starting with `#` after a value are comments.

## Parameters

### Sampling Method

```
Sampling 1
```

| Value | Method | Description |
|-------|--------|-------------|
| 1 | REMC | Replica Exchange Monte Carlo (recommended) |
| 2 | SA | Simulated Annealing |

**Recommendation**: Use REMC (`1`) for most applications. SA is faster but may miss global minima.

### Folding Steps

```
Folding_steps 500000
```

Number of Monte Carlo steps for structure sampling.

| Use Case | Recommended Steps |
|----------|-------------------|
| Quick test | 100,000 - 150,000 |
| Structure prediction | 500,000 - 750,000 |
| Thermal stability | 2,000,000 - 4,000,000 |
| Publication quality | ≥ 4,000,000 |

### Optimization Steps

```
Optimizing_steps 100000
```

Number of energy minimization steps after sampling.

| Use Case | Recommended Steps |
|----------|-------------------|
| Quick test | 50,000 |
| Standard | 100,000 - 200,000 |
| High quality | ≥ 500,000 |

### Ion Concentrations

```
CNa 1000    # Na+ in mM
CMg 0       # Mg2+ in mM
```

Ionic conditions affect DNA stability and structure.

| Condition | CNa (mM) | CMg (mM) |
|-----------|----------|----------|
| Physiological | 150 | 1-2 |
| Standard buffer | 100-200 | 0-10 |
| High salt | 1000 | 0 |
| Mg2+ stabilization | 0 | 10-20 |

### Number of Structures

```
Ncout 10
```

Number of lowest-energy structures to output (1-20 recommended).

## Python Configuration

Use `FoldingConfig` for programmatic configuration:

```python
from dnafold2 import FoldingConfig

# Create configuration
config = FoldingConfig(
    sampling_method="remc",    # "remc" or "sa"
    folding_steps=750000,
    optimizing_steps=200000,
    na_concentration=1000.0,
    mg_concentration=0.0,
    n_structures=10
)

# Save to file
config.to_file("my_config.dat")

# Load from file
config = FoldingConfig.from_file("config.dat")

# View as dictionary
print(config.to_dict())
```

## CLI Configuration Options

Override configuration via command line:

```bash
dnafold2 fold \
    --sequence ATCCTAGTTATAGGAT \
    --method remc \
    --steps 500000 \
    --na-concentration 1000 \
    --mg-concentration 0 \
    --n-structures 10 \
    --output results/
```

## Example Configurations

### Quick Test

```
Sampling 1
Folding_steps 100000
Optimizing_steps 50000
CNa 1000
CMg 0
Ncout 5
```

### 3D Structure Prediction

```
Sampling 1
Folding_steps 750000
Optimizing_steps 200000
CNa 1000
CMg 0
Ncout 10
```

### Thermal Stability Analysis

```
Sampling 1
Folding_steps 4000000
Optimizing_steps 500000
CNa 1000
CMg 0
Ncout 10
```

### Physiological Conditions

```
Sampling 1
Folding_steps 750000
Optimizing_steps 200000
CNa 150
CMg 2
Ncout 10
```
