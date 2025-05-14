# File: utils/delft3d_rgh.py
# Author: Shabnam
# Modified: April 2025

import numpy as np

def create_rgh_from_dep(dep_filename, rgh_filename='output.rgh', roughness_value=22.0):
    with open(dep_filename, 'r') as f:
        lines = f.readlines()

    # Extract all numbers into a flat list
    numbers = []
    for line in lines:
        if line.strip() and not line.startswith('*'):
            parts = line.strip().split()
            numbers.extend([float(p) for p in parts])

    total_values = len(numbers)
    grid_size = int(np.sqrt(total_values))  # Assuming square grid

    if grid_size ** 2 != total_values:
        raise ValueError(f"DEP file does not contain a square grid: {total_values} values found.")

    # Reshape into 2D array
    grid_array = np.array(numbers).reshape((grid_size, grid_size))

    N, M = grid_array.shape  # Rows (N), Columns (M)

    # Create roughness arrays
    rgh_M = np.full((N, M), roughness_value)
    rgh_N = np.full((N, M), roughness_value)

    # Stack for Delft3D format: M-direction first, then N-direction
    rgh_combined = np.vstack([rgh_M, rgh_N])

    # Write .rgh file
    with open(rgh_filename, 'w') as f:
        for row in rgh_combined:
            f.write(' '.join(f'{val:.4f}' for val in row) + '\n')

    print(f"Roughness file '{rgh_filename}' created based on '{dep_filename}' with grid size {M}x{N}.")

# Example usage
create_rgh_from_dep('45x45.dep', rgh_filename='45x45.rgh', roughness_value=22.0)
