# File: src_delft3d4/delft3d4_rgh.py
# Author: Shabnam
# Modified: May 2025
# Purpose: Create a Manning's n .rgh file using .grd dimensions (for non-square grids)


import numpy as np
import os

def read_grid_dimensions(grd_file):
    """
    Reads NX and NY from line 7 of the .grd file.
    Returns number of nodes in X and Y directions (NX, NY).
    """
    with open(grd_file, 'r') as f:
        lines = f.readlines()
    parts = lines[6].strip().split()
    NX, NY = int(parts[0]), int(parts[1])
    return NX, NY

def create_rgh_from_grd(grd_file, roughness_value=0.02):
    """
    Creates a .rgh file using constant Manning’s n based on grid dimensions from .grd file.
    
    Parameters:
    - grd_file: str, path to input .grd file
    - roughness_value: float, constant Manning's n value to apply (default = 0.02)
    """
    NX, NY = read_grid_dimensions(grd_file)

    # Create output filename by replacing .grd with .rgh
    base, _ = os.path.splitext(grd_file)
    output_rgh_file = f"{base}.rgh"

    # Each .rgh block is (NX+1) x (NY+1)
    ncols = NX + 1
    nrows = NY + 1
    
    # Delft3D expects .rgh files to contain two blocks:
    #   - The first block (M-direction roughness) and 
    #   - The second block (N-direction roughness)
    # Each block has shape (NY+1, NX+1), so total rows = 2 * (NY+1)
    # Even if the values are the same, both blocks must be present or the model will throw an error.

    rgh_M = np.full((nrows, ncols), roughness_value)
    rgh_N = np.full((nrows, ncols), roughness_value)
    rgh_combined = np.vstack([rgh_M, rgh_N])

    with open(output_rgh_file, 'w') as f:
        for row in rgh_combined:
            f.write(" ".join(f"{val:.4f}" for val in row) + "\n")

    print(f"? .rgh file written: {output_rgh_file}  ?  shape: {rgh_combined.shape} (M + N directions)")

# Only runs when the script is executed directly
if __name__ == "__main__":
    create_rgh_from_grd("45x45.grd", roughness_value=0.02)
