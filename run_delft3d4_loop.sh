#!/bin/bash

# File: run_delft3d4_loop.sh
# Author: Shabnam
# Modified: April 2025

# Activate Delft3D environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /usr/local/delft3d

# Add Delft3D binaries to PATH
export PATH=/usr/local/delft3d4/bin:$PATH

# Set working directory
BASE_DIR=$(pwd)

# Number of years to simulate
NUM_YEARS=2

# Create output directory
mkdir -p "$BASE_DIR/yearly_outputs"

# Initial .dep file
DEP_FILE="45x45.dep"

for year in $(seq 1 $NUM_YEARS); do
    echo "================ YEAR ${year} ================"

    # Run Delft3D
    echo "Running Delft3D..."
    d_hydro "$BASE_DIR/config_dflow2d3d.xml"

    # Archive Delft3D outputs
    mkdir -p "$BASE_DIR/yearly_outputs/year${year}"
    cp "$BASE_DIR"/tri*-45x45* "$BASE_DIR/yearly_outputs/year${year}/" 2>/dev/null || echo "No tri-* files found to copy"
    
    # Run MATLAB water level extraction
    echo "Running MATLAB script..."
    export PATH=/usr/local/MATLAB/R2024b/bin:$PATH
    matlab -batch "run('${BASE_DIR}/src_delft3d4/waterlevel_delft3d4.m')"

    # Run WEADS Point Coupling for updated topo
    echo "Running WEADS_delft3d4.py..."
    conda activate WEADS_delft3d4_env
    
    python WEADS_delft3d4.py \
      --inputGrdFile 45x45.grd \
      --inputEncFile 45x45.enc \
      --inputDepFile "$DEP_FILE" \
      --updatedRghFile updated_45x45.rgh \
      --inputMdfFile 45x45.mdf \
      --wl_csv water_level_extracted.csv \
      --inEPSG 26914 \
      --outEPSG 26914 \
      --deltaT 1 \
      --outputMEMFile ecology.csv \
      --updatedDepFile updated_45x45.dep \
      --all

    # Replace .dep and .rgh file for next year
    cp updated_45x45.dep 45x45.dep
    cp updated_45x45.rgh 45x45.rgh
    DEP_FILE="45x45.dep"

    # Archive outputs
    cp ecology.csv "$BASE_DIR/yearly_outputs/year${year}/"
    cp updated_45x45.dep "$BASE_DIR/yearly_outputs/year${year}/"
    cp updated_45x45.rgh "$BASE_DIR/yearly_outputs/year${year}/"
    
    # Move MATLAB-generated NetCDF to correct folder if exists
    if [ -f "$BASE_DIR/trim_45x45.nc" ]; then
        cp "$BASE_DIR/trim_45x45.nc" "$BASE_DIR/yearly_outputs/year${year}/"
    fi
      
done

echo "All years complete. Outputs saved in yearly_outputs/"

echo "plotting...saved in yearly_outputs/"
conda activate WEADS_delft3d4_env
python src_delft3d4/plots_delft3d4.py

echo "Done"
