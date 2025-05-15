# File: src_point/general_functions_delft3d4.py
# Author: Shabnam
# Modified: May 2025

import numpy as np

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
