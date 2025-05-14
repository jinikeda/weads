#!/bin/bash

# File: clean_run.sh
# Purpose: Delete all intermediate files to ensure a clean rerun

echo "🧹 Cleaning up files for a fresh run..."

# Remove Delft3D output and restart files
rm -f trim-*.dat trim-*.def trih-*.dat trih-*.def
rm -f tri-diag.* tri-rst.* tri-ini.* tri-var.*

# Remove NetCDF and CSV outputs
rm -f trim_*.nc water_level_extracted.csv grid_coords_extracted.csv tidal_metrics.csv

# Remove MEM and coupling outputs
rm -f ecology.csv summary_center_plots.png

# Remove updated depth and roughness files
rm -f updated*.dep updated*.rgh

# Remove yearly output folder
rm -rf yearly_outputs

# Restore original model files
if [ -d "example_delft3d4" ]; then
    echo "🔄 Restoring initial model files from 'example_delft3d4/'..."
    cp -f example_delft3d4/* .
else
    echo "⚠️  Warning: 'example_delft3d4/' folder not found!"
fi

echo "✅ Cleanup complete. Ready for a fresh run!"
