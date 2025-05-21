# File: src_point/general_functions_delft3d4.py
# Author: Shabnam
# Modified: May 2025

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
    
def read_rgh_file(rgh_file):
    """
    Reads a Delft3D .rgh file and returns a 1D array of Manning's n values.
    Assumes the file has 2*N rows (M-direction + N-direction), where only one is needed.
    """
    values = []
    with open(rgh_file, 'r') as f:
        for line in f:
            if line.strip():
                values.extend([float(val) for val in line.strip().split()])
    half = len(values) // 2
    return np.array(values[:half])  # Return only M-direction values (first half)


def read_grid_shape_from_grd(grd_file):
    """
    Reads the number of grid lines in x (NX) and y (NY) from line 7 of a Delft3D .grd file.
    Returns (NX+1, NY+1) for full node dimensions including boundary rows/cols.
    """
    with open(grd_file, 'r') as f:
        lines = f.readlines()

    if len(lines) < 7:
        raise ValueError(f"Grid file {grd_file} is too short to read NX and NY.")

    tokens = lines[6].strip().split()
    if len(tokens) < 2:
        raise ValueError(f"Line 7 of {grd_file} does not contain two integers.")

    nx, ny = int(tokens[0]), int(tokens[1])
    return nx + 1, ny + 1  # because .dep includes boundary nodes

def infer_grid_shape_from_grd(grd_file):
    """
    Wrapper to read grid shape using the .grd file for use in .dep logic.
    """
    return read_grid_shape_from_grd(grd_file)
