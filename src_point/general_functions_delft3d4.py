# File: src_point/general_functions_delft3d4.py
# Author: Shabnam and Jin Ikeda
# Modified: May 2025

import numpy as np
import pandas as pd
import math
import re
import logging


def read_grid_file(grd_file, output_Flag=False):
    """
    Reads a .grd file and returns:
    - df: DataFrame with columns ['node_id', 'x', 'y'] (includes ghost cells if enabled)
    - nx, ny: number of grid nodes in x and y directions

    Parameters:
    - grd_file (str): Path to the .grd file
    - output_flag (bool): Whether to output a CSV file with the parsed data
    """
    # Read the grid file and extract metadata  
    with open(grd_file, 'r') as f:
        lines = f.readlines()

    if len(lines) < 7:
        raise ValueError(f"Grid file '{grd_file}' is too short to contain valid metadata.")

    # Extract grid dimensions from line 7 (index 6)
    try:
        nx, ny = map(int, lines[6].strip().split())
    except Exception as e:
        raise ValueError(f"Failed to parse grid dimensions from line 7: {e}")

    node_count = nx * ny
    data_lines = lines[8:]  # Skip to data section

    # Initialize lists to store x and y data
    x_data = []
    y_data = []
    current_block = []
    blocks = []
    eta_ids = []

    data_lines = lines[8:]  # Data starts at line 9

    for line in data_lines:
        line = line.strip()
        if line.startswith("ETA="):
            if current_block:
                blocks.append(current_block)
                current_block = []

            try:
                eta_id = int(line.split()[1])
                eta_ids.append(eta_id)
            except Exception:
                eta_ids.append("??")

            parts = line.split()[2:]
            values = [float(v) for v in parts if any(c in v for c in "E.-0123456789")]
            current_block.extend(values)
        else:
            values = [float(v) for v in line.split() if any(c in v for c in "E.-0123456789")]
            current_block.extend(values)

    if current_block:
        blocks.append(current_block)

    expected_blocks = ny + nx
    if len(blocks) != expected_blocks:
        print(f"Found ETA block indices: {eta_ids}")
        raise ValueError(f"Expected {expected_blocks} blocks (X rows + Y cols), but got {len(blocks)}")

    # First ny blocks → x_data, each block = 1 row with nx values
    for b in blocks[:ny]:
        if len(b) != nx:
            raise ValueError(f"Expected {nx} values in X-row block, got {len(b)}")
        x_data.extend(b)

    # Next nx blocks → y_data, each block = 1 column with ny values
    for b in blocks[ny:]:
        if len(b) != ny:
            raise ValueError(f"Expected {ny} values in Y-col block, got {len(b)}")
        y_data.extend(b)

    if len(x_data) != node_count:
        raise ValueError(f"x_data length mismatch: got {len(x_data)}, expected {node_count}")
    if len(y_data) != node_count:
        raise ValueError(f"y_data length mismatch: got {len(y_data)}, expected {node_count}")

    x_array = np.array(x_data).reshape((ny, nx))
    y_array = np.array(y_data).reshape((ny, nx))

    # Reshape x_array and y_array to cordinate system (nx * ny, 2)
    xy_array = np.column_stack(( x_array.ravel(), y_array.ravel()))  # x first then y
    logging.debug(f"Preview of xy_array (first 5 rows): {xy_array[:5]}")
   
    # Reshape to 3D: (ny, nx, 2)
    xy_rows = xy_array.reshape(ny, nx, 2)

    # Prepare final list to collect rows + ghost rows
    rows_with_ghosts = []
    # add ghost cell in the first row
    ghost = np.array([[0.0, 0.0]])  # shape (1, 2)

    for i in range(ny):
        # Add ghost cell to the end of each row
        row = xy_rows[i]                 # shape (nx, 2)
        rows_with_ghosts.append(np.vstack([row, ghost]))  # shape (nx+1, 2)        

    # Add the last full row of ghost cells (nx+1 times [0, 0])
    ghost_row = np.tile(np.array([[0.0, 0.0]]), (nx + 1, 1))
    rows_with_ghosts.append(ghost_row)

        # Stack all into final array
    xy_final = np.vstack(rows_with_ghosts)
    assert xy_final.shape[0] == (ny + 1) * (nx + 1), f"Final shape mismatch: {xy_final.shape[0]} != {(ny + 1) * (nx + 1)}"

    # Create DataFrame with columns 'x' and 'y'
    df = pd.DataFrame(xy_final, columns=['x', 'y'])
    df["node_id"] = df.index + 1

    # Reorder columns to put node_id first
    df = df[['node_id', 'x', 'y']]

    if output_Flag:
        # Save DataFrame to CSV
        output_file = grd_file.replace(".grd", "_xyz.csv")
        df.to_csv(output_file, index=False)
        print(f"Grid data saved to {output_file}")

    return df, nx, ny

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
    Reads a Delft3D .rgh file and returns a 1D array of Manning's n values for the M-direction only.

    Notes:
    - The .rgh file has 2 blocks of data: M-direction (top) and N-direction (bottom).
    - Each block has shape (NY+1, NX+1), so total number of values = 2 * grid_size.
    - For preprocessing, we only use the M-direction block (first half), assuming symmetry.
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

