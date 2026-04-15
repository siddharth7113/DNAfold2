#!/bin/bash
#===============================================================================
# DNAfold2 SLURM Array Job Script
# For batch processing multiple DNA sequences
#===============================================================================

#SBATCH --job-name=dnafold2_batch
#SBATCH --output=logs/dnafold2_%A_%a.out
#SBATCH --error=logs/dnafold2_%A_%a.err
#SBATCH --array=1-10%5              # Process 10 sequences, 5 at a time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=standard

#===============================================================================
# Configuration
#===============================================================================

# Directory containing sequence files (seq_001.dat, seq_002.dat, ...)
SEQUENCES_DIR="${SLURM_SUBMIT_DIR}/sequences"

# Base output directory
OUTPUT_BASE="${SLURM_SUBMIT_DIR}/batch_results_${SLURM_ARRAY_JOB_ID}"

# Sampling parameters
METHOD="remc"
FOLDING_STEPS=750000
NA_CONCENTRATION=1000
MG_CONCENTRATION=0
N_STRUCTURES=10

#===============================================================================
# Environment Setup
#===============================================================================

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Create logs directory
mkdir -p "${SLURM_SUBMIT_DIR}/logs"

#===============================================================================
# Process Current Array Task
#===============================================================================

# Get sequence file for this array task
TASK_ID=${SLURM_ARRAY_TASK_ID}
SEQ_FILE=$(printf "%s/seq_%03d.dat" "${SEQUENCES_DIR}" "${TASK_ID}")

if [ ! -f "${SEQ_FILE}" ]; then
    echo "Error: Sequence file not found: ${SEQ_FILE}"
    exit 1
fi

SEQUENCE=$(cat "${SEQ_FILE}" | tr -d '\n')
OUTPUT_DIR="${OUTPUT_BASE}/seq_$(printf '%03d' ${TASK_ID})"

echo "=============================================="
echo "DNAfold2 Batch Processing"
echo "=============================================="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${TASK_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Sequence file: ${SEQ_FILE}"
echo "Sequence: ${SEQUENCE:0:50}..."
echo "Output: ${OUTPUT_DIR}"
echo "Start time: $(date)"
echo "=============================================="

mkdir -p "${OUTPUT_DIR}"

# Run DNAfold2
dnafold2 fold \
    --sequence "${SEQUENCE}" \
    --output "${OUTPUT_DIR}" \
    --method ${METHOD} \
    --steps ${FOLDING_STEPS} \
    --na-concentration ${NA_CONCENTRATION} \
    --mg-concentration ${MG_CONCENTRATION} \
    --n-structures ${N_STRUCTURES} \
    --verbose

EXIT_CODE=$?

echo ""
echo "Task ${TASK_ID} completed with exit code: ${EXIT_CODE}"
echo "End time: $(date)"

exit ${EXIT_CODE}
