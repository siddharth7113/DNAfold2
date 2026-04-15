# API Reference

Python API documentation for DNAfold2.

## Core Classes

### DNAFolder

Main class for DNA structure prediction.

```python
from dnafold2 import DNAFolder

folder = DNAFolder(config=None, bin_dir=None, data_dir=None)
```

**Parameters:**
- `config` (FoldingConfig, optional): Folding configuration. Uses defaults if not provided.
- `bin_dir` (str or Path, optional): Directory containing compiled binaries.
- `data_dir` (str or Path, optional): Directory containing data files.

#### Methods

##### fold()

```python
result = folder.fold(sequence, output_dir=None, verbose=False)
```

Fold a DNA sequence to predict its 3D structure.

**Parameters:**
- `sequence` (str): DNA sequence (only A, T, C, G)
- `output_dir` (str or Path, optional): Output directory
- `verbose` (bool): Print progress information

**Returns:** `FoldingResult`

##### fold_batch()

```python
results = folder.fold_batch(sequences, output_dir, verbose=False)
```

Fold multiple DNA sequences.

**Parameters:**
- `sequences` (List[str]): List of DNA sequences
- `output_dir` (str or Path): Base output directory
- `verbose` (bool): Print progress

**Returns:** `List[FoldingResult]`

---

### FoldingConfig

Configuration dataclass for folding parameters.

```python
from dnafold2 import FoldingConfig

config = FoldingConfig(
    sampling_method="remc",
    folding_steps=500000,
    optimizing_steps=100000,
    na_concentration=1000.0,
    mg_concentration=0.0,
    n_structures=10
)
```

**Attributes:**
| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `sampling_method` | str | "remc" | "remc" or "sa" |
| `folding_steps` | int | 150000 | Simulation steps |
| `optimizing_steps` | int | 50000 | Optimization steps |
| `na_concentration` | float | 1000.0 | Na+ concentration (mM) |
| `mg_concentration` | float | 0.0 | Mg2+ concentration (mM) |
| `n_structures` | int | 10 | Number of output structures |

#### Methods

##### from_file()

```python
config = FoldingConfig.from_file("config.dat")
```

Load configuration from a file.

##### to_file()

```python
config.to_file("my_config.dat")
```

Save configuration to a file.

##### to_dict()

```python
d = config.to_dict()
```

Convert to dictionary.

---

### FoldingResult

Dataclass containing folding results.

**Attributes:**
| Attribute | Type | Description |
|-----------|------|-------------|
| `sequence` | str | Input DNA sequence |
| `cg_structures` | List[Path] | Coarse-grained PDB files |
| `all_atom_structures` | List[Path] | All-atom PDB files |
| `secondary_structures` | List[Path] | Secondary structure files |
| `folding_trajectories` | List[Path] | Trajectory files |
| `thermal_stability` | Path or None | Thermal stability data |
| `output_dir` | Path or None | Output directory |
| `elapsed_time` | float | Runtime in seconds |

---

## Convenience Functions

### fold()

```python
from dnafold2 import fold

result = fold(
    sequence="ATCCTAGTTATAGGAT",
    config=None,
    output_dir="results/",
    verbose=True
)
```

Convenience function for quick folding.

---

## Utility Functions

### validate_sequence()

```python
from dnafold2.config import validate_sequence

sequence = validate_sequence("atcg atcg")  # Returns "ATCGATCG"
```

Validate and normalize DNA sequences.

**Raises:** `ValueError` if sequence is invalid.

---

## CLI Entry Points

After installation, the `dnafold2` command is available:

```bash
dnafold2 fold --sequence SEQ [OPTIONS]
dnafold2 config [OPTIONS]
dnafold2 validate SEQUENCE
dnafold2 info
```

See `dnafold2 --help` for full documentation.
