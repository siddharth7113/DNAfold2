import numpy as np
import sys

# Read energy data file
try:
    data = np.loadtxt('Energy_0.dat')
except Exception as e:
    print(f"Warning: Could not load Energy_0.dat: {e}")
    np.savetxt('min.dat', [], fmt='%d')
    sys.exit(0)

# Handle empty or insufficient data
if data.size == 0:
    print("Warning: Energy_0.dat is empty — no conformations to score")
    np.savetxt('min.dat', [], fmt='%d')
    sys.exit(0)

# Handle single-row data (1D array)
if data.ndim == 1:
    data = data.reshape(1, -1)

if data.shape[1] < 2:
    print("Warning: Energy_0.dat has fewer than 2 columns")
    np.savetxt('min.dat', [], fmt='%d')
    sys.exit(0)

# Get second column (energy) and sort to find lowest-energy conformations
second_column = data[:, 1]
indices = np.argsort(second_column)

# Take up to 500 lowest-energy conformations
n = min(500, len(indices))
top_indices = indices[:n]

# Save indices to file (1-based line numbers)
np.savetxt('min.dat', top_indices + 1, fmt='%d')
print(f"Selected {n} lowest-energy conformations")
