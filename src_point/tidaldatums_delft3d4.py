# File: src_point/tidaldatums_delft3d4.py
# Author: Shabnam
# Purpose: Compute MHW, MLW, and percent inundation from wide-format CSV using Delft3D-MATLAB 
# Updated: May 2025 (Fix: properly match pt IDs to grid coords)

import pandas as pd
import numpy as np

def calculate_tidal_metrics_from_csv(water_level_csv, coords_csv, output_csv='tidal_metrics.csv', threshold=0.0):
    """
    Compute MHW, MLW, and percent inundation per grid cell.

    Parameters:
    - water_level_csv: CSV with time in first column, water level at pt1, pt2,... in subsequent columns.
    - coords_csv: grid coordinates with columns: pt,x,y
    - output_csv: where to save final tidal metrics
    - threshold: elevation threshold for percent inundation (default 0.0 m)

    Returns:
    - pandas DataFrame with MHW, MLW, percent_inundation per grid point
    """
    print(f"Reading wide-format water level CSV: {water_level_csv}")
    df = pd.read_csv(water_level_csv)
    time = df.iloc[:, 0]
    water_data = df.iloc[:, 1:]  # pt1, pt2, ...
    print(water_data.columns)

    print(f"Reading grid node coordinates: {coords_csv}")
    coords_df = pd.read_csv(coords_csv)
    z = coords_df['z'].values.T
    print
    # coords_df = coords_df.set_index('pt')  # 'pt' is assumed to be an integer

    results = []
    for i, col in enumerate(water_data.columns):
        wl = water_data[col].values
        wl = wl[np.isfinite(wl)]
        z_local = z[i]  # z is a 1D array, so we can index it directly
        if len(wl) < 3:
            mhw = np.nan
            mlw = np.nan
            percent_inundation = np.nan
        else:
            mask_max = (wl[1:-1] > wl[:-2]) & (wl[1:-1] > wl[2:])
            mask_min = (wl[1:-1] < wl[:-2]) & (wl[1:-1] < wl[2:])
            maxima = wl[1:-1][mask_max]
            minima = wl[1:-1][mask_min]

            mhw = np.mean(maxima) if len(maxima) >= 3 else np.nanmax(wl)
            mlw = np.mean(minima) if len(minima) >= 3 else np.nanmin(wl)
            
            #################################
            # wrong. # percent_inundation = np.sum(wl > z) / len(wl)
            # percent_inundation = np.sum(wl > threshold) / len(wl)
            percent_inundation = np.sum(wl > z_local) / len(wl)
            #################################

        # # Fix: Convert column name 'pt1' -> 1 (int)
        # try:
        #     pt_number = int(col.replace('pt', ''))
        #     x = coords_df.at[pt_number, 'x']
        #     y = coords_df.at[pt_number, 'y']
        # except (KeyError, ValueError):
        #     x = np.nan
        #     y = np.nan

        results.append({
            # 'grid_id': col,  # keep 'pt1', 'pt2', ...
            # 'x': x,
            # 'y': y,
            'MHW': mhw,
            'MLW': mlw,
            'percent_inundation': percent_inundation
        })

    df_out = coords_df.copy()
    df_result  = pd.DataFrame(results)
    df_out = pd.concat([df_out, df_result], axis=1)
    df_out.to_csv(output_csv, index=False)
    print(f"✓ Tidal metrics saved to: {output_csv}")
    return df_out


if __name__ == '__main__':
    calculate_tidal_metrics_from_csv('water_level_extracted.csv')
