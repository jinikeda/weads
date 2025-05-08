# File: src_point/postprocessing_delft3d4.py
# Author: Shabnam
# Updated: July 2025 (Fix .dep update mismatch and clean logic)

import pandas as pd
import numpy as np
import shutil
import os
from .basics_delft3d4 import fileexists

def infer_grid_shape_from_dep(dep_file):
    with open(dep_file, 'r', encoding='utf-8') as f:
        values = [float(val) for line in f if line.strip() for val in line.strip().split()]

    total = len(values)
    side = int(np.sqrt(total))
    if side * side != total:
        raise ValueError(f"Grid not square: {total} values → {side}x{side} is invalid")

    return side, side

def update_delft3d_dep(original_dep_file, output_dep_file, tb_update_array, values_per_line=12):
    fileexists(original_dep_file)

    with open(original_dep_file, 'r', encoding='utf-8') as f:
        values = [float(val) for line in f if line.strip() for val in line.strip().split()]

    grid_size = int(np.sqrt(len(values)))
    z_matrix = np.array(values).reshape((grid_size, grid_size))
    z_flat = z_matrix.flatten()

    valid_mask = z_flat != -999.0

    if tb_update_array.size != np.count_nonzero(valid_mask):
        raise ValueError(f"Mismatch: tb_update ({tb_update_array.size}) vs valid .dep cells ({np.count_nonzero(valid_mask)})")

    z_flat[valid_mask] = tb_update_array
    updated_matrix = z_flat.reshape(z_matrix.shape)

    with open(output_dep_file, 'w', encoding='utf-8') as f:
        for row in updated_matrix:
            for j, val in enumerate(row):
                val_str = f'{val:.7E}'
                spacing = '  ' if val < 0 else '   '
                f.write(spacing + val_str)
                if (j + 1) % values_per_line == 0 and j != len(row) - 1:
                    f.write('\n')
            f.write('\n')

    print(f"✅ Updated .dep file written in Quickin style: {output_dep_file}")

def update_delft3d_rgh(outputMEMFile, grid_shape, output_rgh_file):
    df = pd.read_csv(outputMEMFile)

    class_to_chezy = {
        55: 65.0,   # land / upland
        40: 20.0,   # subtidal or water
        16: 45.0,   # low productivity
        23: 35.0,   # medium productivity
        32: 25.0,   # high productivity
        128: 65.0   # nodata default
    }

    if 'P' not in df.columns:
        raise KeyError("Column 'P' not found in ecology file.")

    chezy_vals = df['P'].map(class_to_chezy).fillna(65.0).values

    expected_size = grid_shape[0] * grid_shape[1]
    if chezy_vals.size != expected_size:
        raise ValueError(f"Mismatch: {chezy_vals.size} values vs grid size {expected_size}.")

    chezy_grid = chezy_vals.reshape(grid_shape)
    rgh_combined = np.vstack([chezy_grid, chezy_grid])

    with open(output_rgh_file, 'w') as f:
        for row in rgh_combined:
            f.write("   ".join(f"{val:.4f}" for val in row) + "\n")

    print(f"✅ .rgh file created: {output_rgh_file} (shape: {rgh_combined.shape})")

def update_delft3d_mdf(mdf_file, new_dep_file, new_rgh_file):
    fileexists(mdf_file)
    with open(mdf_file, 'r') as f:
        lines = f.readlines()

    found_dep = found_rgh = found_roumet = False
    updated = []

    for line in lines:
        if line.lower().startswith("fildep"):
            updated.append(f"Fildep = #{new_dep_file}#\n"); found_dep = True
        elif line.lower().startswith("filrgh"):
            updated.append(f"Filrgh = #{new_rgh_file}#\n"); found_rgh = True
        elif line.lower().startswith("roumet"):
            updated.append("Roumet = #C#\n"); found_roumet = True
        else:
            updated.append(line)

    if not found_dep:
        updated.append(f"Fildep = #{new_dep_file}#\n")
    if not found_rgh:
        updated.append(f"Filrgh = #{new_rgh_file}#\n")
    if not found_roumet:
        updated.append("Roumet = #C#\n")

    with open(mdf_file, 'w') as f:
        f.writelines(updated)

    print(f"✔ Updated .mdf file to point to: {new_dep_file} and {new_rgh_file}")

def postprocessing_delft(inputDepFile, inputMdfFile, outputMEMFile,
                         dep_out='updated.dep', rgh_out='updated.rgh'):
    print("\nPOSTPROCESSING Delft3D-MEM Coupling...")
    fileexists(outputMEMFile)
    df = pd.read_csv(outputMEMFile)
    print(f"✔ Loaded MEM results: {outputMEMFile}  →  shape {df.shape}")

    tb_all = df['tb_update'].values
    grid_shape = infer_grid_shape_from_dep(inputDepFile)
    shutil.copy(inputDepFile, dep_out)

    # Load .dep values and create valid mask
    with open(inputDepFile, 'r', encoding='utf-8') as f:
        dep_vals = [float(val) for line in f if line.strip() for val in line.strip().split()]
    dep_array = np.array(dep_vals)
    valid_mask = dep_array != -999.0

    # Apply mask to tb_update
    tb_valid = tb_all[valid_mask]

    update_delft3d_dep(inputDepFile, dep_out, tb_valid)
    update_delft3d_rgh(outputMEMFile, grid_shape, rgh_out)
    update_delft3d_mdf(inputMdfFile, dep_out, rgh_out)

    print("✔ Postprocessing complete.\n")
