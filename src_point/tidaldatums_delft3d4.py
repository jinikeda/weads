#!/usr/bin/env python3
# File: src_point/tidaldatums_delft3d4.py
# Purpose: Compute MHW, MLW, and percent inundation from wide-format CSV using /Delft3D-MATLAB 
# Updated: May 2025

import pandas as pd
import numpy as np

def calculate_tidal_metrics_from_csv(water_level_csv, coords_csv='grid_coords_extracted.csv', output_csv='tidal_metrics.csv', threshold=0.0):
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

    print(f"Reading grid node coordinates: {coords_csv}")
    coords_df = pd.read_csv(coords_csv)

    results = []
    for i, col in enumerate(water_data.columns):
        grid_id = i + 1  # pt1 → grid_id 1
        wl = water_data[col].values

        # Clean NaNs
        wl = wl[np.isfinite(wl)]

        if len(wl) < 3:
            mhw = np.nan
            mlw = np.nan
            percent_inundation = np.nan
        else:
            # Detect peaks and troughs
            mask_max = (wl[1:-1] > wl[:-2]) & (wl[1:-1] > wl[2:])
            mask_min = (wl[1:-1] < wl[:-2]) & (wl[1:-1] < wl[2:])
            maxima = wl[1:-1][mask_max]
            minima = wl[1:-1][mask_min]

            if len(maxima) >= 3 and len(minima) >= 3:
                mhw = np.mean(maxima)
                mlw = np.mean(minima)
            else:
                # Fallback to max/min
                mhw = np.nanmax(wl)
                mlw = np.nanmin(wl)

            percent_inundation = np.sum(wl > threshold) / len(wl)

        x = coords_df.iloc[i]['x'] if i < len(coords_df) else np.nan
        y = coords_df.iloc[i]['y'] if i < len(coords_df) else np.nan

        results.append({
            'grid_id': grid_id,
            'x': x,
            'y': y,
            'MHW': mhw,
            'MLW': mlw,
            'percent_inundation': percent_inundation
        })

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_csv, index=False)
    print(f"✓ Tidal metrics saved to: {output_csv}")
    return df_out


if __name__ == '__main__':
    calculate_tidal_metrics_from_csv('water_level_extracted.csv')
