# File: src_point/preprocessing_delft3d4.py
# Author: Shabnam
# Modified: June 2025

import time
import numpy as np
import pandas as pd
from .basics import fileexists
import xarray as xr
from scipy.spatial import cKDTree
from src_point.general_functions_delft3d4 import read_dep_file


def idw_interpolate(x_known, y_known, values, x_target, y_target, k=6, power=2):
    """
    Interpolate missing values using Inverse Distance Weighting (IDW).
    """
    tree = cKDTree(np.column_stack((x_known, y_known)))
    distances, idxs = tree.query(np.column_stack((x_target, y_target)), k=k)

    weights = 1.0 / (distances ** power + 1e-12)
    interpolated = np.sum(weights * values[idxs], axis=1) / np.sum(weights, axis=1)
    return interpolated

def preprocessing_Delft3D(
    inputBathymetryFile,
    inputWaterLevelCSV,
    inputShapeFile,         # not used yet
    domainIOFile,
    inEPSG,                 # not used yet
    outEPSG                 # not used yet
):
    start_time = time.time()
    print("\nLAUNCH: Running Delft3D Preprocessing Script with WEADS-style HydroClass + IDW Interpolation\n")

    # --- Check required files ---
    fileexists(inputBathymetryFile)
    fileexists(inputWaterLevelCSV)

    # --- Load .dep file into flat z array ---
    from .general_functions_delft3d4 import read_dep_file
    z = read_dep_file(inputBathymetryFile)
    print(f"✔ Read {len(z)} bathymetry points from {inputBathymetryFile}")

    # --- Load tidal metrics from NetCDF or CSV ---
    if inputWaterLevelCSV.endswith('.nc'):
        print("✔ Detected NetCDF file input for water levels...")
        ds = xr.open_dataset(inputWaterLevelCSV)
        mhw_array = ds['MHW'].values
        mlw_array = ds['MLW'].values
        x = ds['x'].values
        y = ds['y'].values
    else:
        print("Detected CSV file input for water levels...")
        df = pd.read_csv(inputWaterLevelCSV)
        mhw_array = df['MHW'].values
        mlw_array = df['MLW'].values
        x = df['x'].values
        y = df['y'].values

    # --- Interpolate missing MHW/MLW using IDW ---
    valid_mask = np.isfinite(mhw_array) & np.isfinite(mlw_array)
    x_known = x[valid_mask]
    y_known = y[valid_mask]
    mhw_known = mhw_array[valid_mask]
    mlw_known = mlw_array[valid_mask]
    x_missing = x[~valid_mask]
    y_missing = y[~valid_mask]

    if len(x_missing) > 0:
        print(f"Interpolating {len(x_missing)} missing MHW/MLW values using IDW...")
        mhw_interp = idw_interpolate(x_known, y_known, mhw_known, x_missing, y_missing)
        mlw_interp = idw_interpolate(x_known, y_known, mlw_known, x_missing, y_missing)
        mhw_array[~valid_mask] = mhw_interp
        mlw_array[~valid_mask] = mlw_interp

    # --- Check shape consistency ---
    assert len(z) == len(mhw_array), f"Length mismatch: bathymetry={len(z)}, MHW={len(mhw_array)}"

    # --- Compute HydroClass and assign labels ---
    hydroclass_index = np.full_like(z, 1, dtype=int)  # default to intertidal
    hydroclass_index[z > mhw_array] = 0  # dry
    hydroclass_index[z <= mlw_array] = 2  # submerged

    hydroclass_label = np.array(['intertidal'] * len(z), dtype=object)
    hydroclass_label[hydroclass_index == 0] = 'land'
    hydroclass_label[hydroclass_index == 2] = 'subtidal'

    # --- Fill manning and other values ---
    mannings_n = np.full_like(z, 0.035, dtype=float)

    df_out = pd.DataFrame({
        'node': np.arange(1, len(z) + 1),
        'x': x,
        'y': y,
        'z': z,
        'mann': mannings_n,
        'MLW_IDW': mlw_array,
        'MHW_IDW': mhw_array,
        'HydroClass': hydroclass_label
    })

    # --- Save output ---
    df_out.to_csv(domainIOFile, index=False)
    print(f"\n Saved domain input to {domainIOFile}")
    print(f"Done preprocessing. Time: {time.time() - start_time:.2f} seconds\n")
