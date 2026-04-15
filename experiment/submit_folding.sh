#!/bin/bash
#===============================================================================
# DNAfold2 SLURM Job Submission Script
# For single DNA sequence structure prediction
#===============================================================================

#SBATCH --job-name=dnafold2_2rqj
#SBATCH --output=dnafold2_%j.out
#SBATCH --error=dnafold2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=main


#===============================================================================
# Configuration - Edit these parameters
#===============================================================================

# Input sequence (either inline or from file)
SEQUENCE="GGAGGAGGAGGA "
# Or use: SEQUENCE=$(cat /path/to/seq.dat)

# Configuration file (optional - uses defaults if not specified)
CONFIG_FILE=""

# Output directory
OUTPUT_DIR="${SLURM_SUBMIT_DIR}/results_${SLURM_JOB_ID}"

# Sampling method: "remc" or "sa"
METHOD="remc"

# Number of folding steps (recommended: >=500000 for structure, >=2000000 for stability)
FOLDING_STEPS=100000

# Ion concentrations (mM)
NA_CONCENTRATION=137
MG_CONCENTRATION=0

# Number of structures to predict
N_STRUCTURES=10

# Number of temperature replicas (should match cpus-per-task for REMC)
N_THREADS=${SLURM_CPUS_PER_TASK:-4}

#===============================================================================
# Environment Setup
#===============================================================================

# Activate conda/virtual environment if needed
source /home/siddharth/work/DNAfold2/.venv/bin/activate
# conda activate dnafold2

# Set OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Path to DNAfold2 installation
DNAFOLD2_DIR="${SLURM_SUBMIT_DIR}/.."

#===============================================================================
# Run Folding
#===============================================================================

echo "=============================================="
echo "DNAfold2 Structure Prediction"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start time: $(date)"
echo ""
echo "Sequence: ${SEQUENCE:0:50}..."
echo "Method: ${METHOD}"
echo "Folding steps: ${FOLDING_STEPS}"
echo "Output: ${OUTPUT_DIR}"
echo "=============================================="

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Build command
CMD="dnafold2 fold"
CMD="${CMD} --sequence ${SEQUENCE}"
CMD="${CMD} --output ${OUTPUT_DIR}"
CMD="${CMD} --method ${METHOD}"
CMD="${CMD} --steps ${FOLDING_STEPS}"
CMD="${CMD} --na-concentration ${NA_CONCENTRATION}"
CMD="${CMD} --mg-concentration ${MG_CONCENTRATION}"
CMD="${CMD} --n-structures ${N_STRUCTURES}"
CMD="${CMD} --n-threads ${N_THREADS}"
CMD="${CMD} --verbose"

if [ -n "${CONFIG_FILE}" ] && [ -f "${CONFIG_FILE}" ]; then
    CMD="${CMD} --config ${CONFIG_FILE}"
fi

# Run DNAfold2
echo ""
echo "Running: ${CMD}"
echo ""

${CMD}
EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Job completed with exit code: ${EXIT_CODE}"
echo "End time: $(date)"
echo "=============================================="

exit ${EXIT_CODE}
