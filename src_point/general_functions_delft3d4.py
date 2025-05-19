# File: src_point/general_functions_delft3d4.py
# Author: Shabnam
# Modified: July 2025

import numpy as np
import math

def read_dep_file(dep_file):
    """
    Reads a Delft3D .dep file and returns a flat numpy array of bathymetry values.
    Skips empty lines and handles multiple values per line.
    """
    values = []
    with open(dep_file, 'r') as f:
        for line in f:
            if line.strip():
                values.extend([float(val) for val in line.strip().split()])
    return np.array(values)


def infer_grid_shape_from_dep(dep_file):
    """
    Infers the shape (nx, ny) of the grid from the .dep file, assuming either a square
    or rectangular layout. Tries to guess the shape based on total number of values.
    """
    values = read_dep_file(dep_file)
    total = len(values)

    # First try perfect square
    side = math.isqrt(total)
    if side * side == total:
        return side, side

    # Then try rectangular factors
    for ny in range(1, total + 1):
        if total % ny == 0:
            nx = total // ny
            if nx * ny == total:
                return nx, ny

    raise ValueError(f"Cannot infer grid shape from {dep_file} with {total} values.")
