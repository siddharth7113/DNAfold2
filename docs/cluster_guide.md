# Cluster Guide (SLURM)

Running DNAfold2 on HPC clusters using SLURM.

## Quick Start

```bash
# Single sequence
sbatch scripts/slurm/submit_folding.sh

# Batch processing
sbatch scripts/slurm/submit_batch.sh
```

## Single Sequence Job

### Configuration

Edit `scripts/slurm/submit_folding.sh`:

```bash
# Resource allocation
#SBATCH --cpus-per-task=16    # Number of CPU cores
#SBATCH --mem=8G              # Memory
#SBATCH --time=24:00:00       # Wall time limit
#SBATCH --partition=standard  # Partition name

# Job parameters
SEQUENCE="ATCCTAGTTATAGGAT"   # Your sequence
METHOD="remc"                  # remc or sa
FOLDING_STEPS=750000          # Simulation steps
```

### Submission

```bash
# Submit job
sbatch scripts/slurm/submit_folding.sh

# Check status
squeue -u $USER

# View output
cat dnafold2_<jobid>.out
```

## Batch Processing (Array Jobs)

For multiple sequences:

### Setup

1. Create sequence files:
```bash
mkdir sequences
echo "ATCCTAGTTATAGGAT" > sequences/seq_001.dat
echo "GCTAGCTAGCTAGCTA" > sequences/seq_002.dat
# ... etc
```

2. Edit `scripts/slurm/submit_batch.sh`:
```bash
#SBATCH --array=1-10%5   # 10 sequences, 5 concurrent

SEQUENCES_DIR="${SLURM_SUBMIT_DIR}/sequences"
```

### Submission

```bash
sbatch scripts/slurm/submit_batch.sh
```

## Resource Guidelines

| Sequence Length | CPUs | Memory | Time (REMC 500k steps) |
|-----------------|------|--------|------------------------|
| < 20 nt | 8 | 4 GB | 2-4 hours |
| 20-50 nt | 16 | 8 GB | 4-12 hours |
| 50-100 nt | 16-32 | 16 GB | 12-24 hours |
| > 100 nt | 32+ | 32 GB | 24-48 hours |

## Environment Setup

Add to your SLURM script:

```bash
# Load modules (adjust for your cluster)
module load gcc/11.2.0
module load python/3.11

# Or activate conda
source /path/to/conda/etc/profile.d/conda.sh
conda activate dnafold2

# Set OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
```

## Email Notifications

```bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@example.com
```

## Monitoring Jobs

```bash
# Queue status
squeue -u $USER

# Job details
scontrol show job <jobid>

# Cancel job
scancel <jobid>

# View output in real-time
tail -f dnafold2_<jobid>.out
```

## Collecting Results

After batch jobs complete:

```bash
# Results are in batch_results_<jobid>/
ls batch_results_*/seq_*/CG_structure/

# Combine all CG structures
cat batch_results_*/seq_*/CG_structure/*.pdb > all_structures.pdb
```

## Troubleshooting

### Job Fails Immediately

```bash
# Check error file
cat dnafold2_<jobid>.err

# Common issues:
# - Module not found → Check module load commands
# - Python not found → Activate conda/venv
# - Binary not found → Check DNAfold2 installation
```

### Out of Memory

Increase memory allocation:
```bash
#SBATCH --mem=16G
```

### Timeout

Increase wall time or reduce folding steps:
```bash
#SBATCH --time=48:00:00
```
