# File: src_point/preprocessing_delft3d4.py
# Author: Shabnam and Jin Ikeda
# Modified: June 2025

import time
import numpy as np
import pandas as pd
from .basics import fileexists
import xarray as xr
from scipy.spatial import cKDTree
from .general_functions_delft3d4 import read_grid_file, read_dep_file, read_rgh_file
from .tidaldatums_delft3d4 import calculate_tidal_metrics_from_csv


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
    inputGrdFile,
    inputDepthFile,
    inputRghFile,
    inputWaterLevelCSV,
    inputShapeFile,         # not used yet
    domainIOFile,
    inEPSG,                 # not used yet
    outEPSG,                # not used yet
):
    start_time = time.time()
    print("\nLAUNCH: Running Delft3D Preprocessing Script with WEADS-style HydroClass + IDW Interpolation\n")

    # --- Check required files ---
    fileexists(inputGrdFile)
    fileexists(inputDepthFile)
    fileexists(inputWaterLevelCSV)
    # fileexists(inputRghFile)  # intentionally not checked now when user does not provide rgh file

    # --- Load grid file into x, y arrays ---
    df, nx, ny = read_grid_file(inputGrdFile, output_Flag=True)
    print(f"✔ Read  {df.shape} grid file...")

    # --- Load .dep file into flat z array ---
    dep = read_dep_file(inputDepthFile)
    z = -dep  # negative depth values (elevation is a positive direction)
    print(f"✔ Read {len(dep)} depth values from {inputDepthFile}")

    df['z'] = z
    ouptutGrdFile = inputGrdFile.replace(".grd", "_xyz.grd")
    df.to_csv(ouptutGrdFile, index=False)
    print(f"✔ Saved grid coordinates to {ouptutGrdFile}")


    # --- Load .rgh file into flat mannings n array ---  #rgh file is optional in delft3d4
    if inputRghFile is not None:
        # Read the .rgh file
        mannings_n = read_rgh_file(inputRghFile)
        print(f"✔ Read {inputRghFile} for roughness values...")
    else:
        # If no .rgh file is provided, create a default array of 0.03 for all nodes
        ########################################################################################################
        # future implementation: read rgh file and assign mannings n from mdf file if rgh file is not provided #
        ########################################################################################################
        mannings_n = np.full(len(z), 0.03)
        print("No roughness file provided. Using default value of 0.03 for all nodes...")

    # --- Load tidal metrics from NetCDF or CSV ---
    # Step 1: Extract tidal datums from CSV
    print("Calculating tidal datums from extracted CSV...")
    tidal_csv = "tidal_metrics.csv"
    coords_csv = ouptutGrdFile
    calculate_tidal_metrics_from_csv(inputWaterLevelCSV,coords_csv, output_csv=tidal_csv)


    if tidal_csv.endswith('.nc'):
        print("✔ Detected NetCDF file input for water levels...")
        ds = xr.open_dataset(tidal_csv)
        mhw_array = ds['MHW'].values
        mlw_array = ds['MLW'].values
        hp_array = ds['hydroperiod'].values

    else:
        print("Detected CSV file input for water levels...")
        df_wl = pd.read_csv(tidal_csv)
        mhw_array = df_wl['MHW'].values
        mlw_array = df_wl['MLW'].values
        hp_array = df_wl['hydroperiod'].values


    # --- Load shape file if provided ---
    # future implementation

    # --- convert to coordinate system if needed ---
    # future implementation


    # --- interpolate MHW/MLW in subtidal and land regions ---
    # Here if hp is less than 1 (not fully submerged), we need to interpolate mhw and mlw values
    valid_mask = hp_array >= 1  # Boolean mask
    # print("Valid indices:", np.where(valid_mask)[0])  # Optional to see the index numbers
    
    np.savetxt("valid_mask.txt", valid_mask)
    print(df.shape)
    x = df['x'].values
    y = df['y'].values
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
        mhw_IDW_array = np.full_like(mhw_array, -9999.0)
        mlw_IDW_array = np.full_like(mlw_array, -9999.0)
        mhw_IDW_array[~valid_mask] = mhw_interp
        mlw_IDW_array[~valid_mask] = mlw_interp

    # --- Check shape consistency ---
    assert len(dep) == len(mhw_array), f"Length mismatch: bathymetry={len(dep)}, MHW={len(mhw_array)}"

    # --- Compute HydroClass and assign labels ---
    hydroclass_index = np.full_like(z, 1, dtype=int)  # default to intertidal
    hydroclass_index[z > mhw_array] = 0  # dry
    hydroclass_index[z <= mlw_array] = 2  # submerged

    hydroclass_label = np.array(['intertidal'] * len(dep), dtype=object)
    hydroclass_label[hydroclass_index == 0] = 'land'
    hydroclass_label[hydroclass_index == 2] = 'subtidal'
    ############################################################################################################


    df_out = pd.DataFrame({
        'node_id': df['node_id'].values,
        'x': x,
        'y': y,
        'z': -dep,
        'mann': mannings_n,
        'hp': hp_array,
        'MLW': mlw_array,
        'MHW': mhw_array,
        'MLW_IDW': mlw_IDW_array,
        'MHW_IDW': mhw_IDW_array,
        'HydroClass': hydroclass_index,
        'HydroLabel': hydroclass_label
    })

    # --- Save output ---
    df_out.to_csv(domainIOFile, index=False)
    print(f"\n Saved domain input to {domainIOFile}")
    print(f"Done preprocessing. Time: {time.time() - start_time:.2f} seconds\n")
