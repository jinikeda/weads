# File: plots_delft3d4.py
# Author: Shabnam
# Modified: April 2025

import os
import pandas as pd
import matplotlib.pyplot as plt

# Setup base_dir relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.join(SCRIPT_DIR, "../yearly_outputs")
years = sorted([d for d in os.listdir(base_dir) if d.startswith("year")])

obs_points = {
    "intertidal": (24, 24),
    "subtidal": (5, 24),
    "supratidal": (36, 24)
}

# Process each observation point
for label, (row_idx, col_idx) in obs_points.items():
    depths, wl_means, mhw_list, mlw_list, biomass, accretion = [], [], [], [], [], []

    for year in years:
        path = os.path.join(base_dir, year, "ecology.csv")
        if not os.path.exists(path):
            print(f"⚠ Missing: {path}")
            continue
        df = pd.read_csv(path)
        center_index = row_idx * 47 + col_idx
        if center_index >= len(df):
            print(f"⚠ Invalid index {center_index} in {path}")
            continue
        row = df.loc[center_index]
        depths.append(row["z"])
        wl_means.append((row["MHW_IDW"] + row["MLW_IDW"]) / 2)
        mhw_list.append(row["MHW_IDW"])
        mlw_list.append(row["MLW_IDW"])
        biomass.append(row["B"])
        accretion.append(row["A"])

    # Plotting
    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    years_axis = list(range(1, len(depths)+1))

    # Depth
    axs[0].plot(years_axis, depths, marker='o', label="Depth")
    axs[0].set_ylabel("Depth (m)")
    axs[0].set_title(f"{label.title()} Cell Elevation Change Over Years")
    axs[0].grid(True)

    # Water level with tidal datums
    axs[1].plot(years_axis, wl_means, marker='o', label="Mean WL")
    axs[1].plot(years_axis, mhw_list, '--', color='red', label="MHW")
    axs[1].plot(years_axis, mlw_list, '--', color='green', label="MLW")
    axs[1].set_ylabel("Water Level (m)")
    axs[1].set_title("Tidal Datums and Mean Water Level")
    axs[1].legend()
    axs[1].grid(True)

    # Biomass
    axs[2].plot(years_axis, biomass, marker='o', color='green', label="Biomass")
    axs[2].set_ylabel("Biomass (g/m²)")
    axs[2].set_title("Biomass Over Time")
    axs[2].grid(True)

    # Accretion
    axs[3].plot(years_axis, accretion, marker='o', color='purple', label="Accretion Rate")
    axs[3].set_ylabel("Accretion (m/year)")
    axs[3].set_xlabel("Year")
    axs[3].set_title("Accretion Rate Over Time")
    axs[3].grid(True)

    plt.tight_layout()
    output_plot = os.path.join(base_dir, f"summary_{label}_plots.png")
    plt.savefig(output_plot)
    print("✅ Saved plot to:", output_plot)
