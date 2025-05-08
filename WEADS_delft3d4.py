#!/usr/bin/env python
# File: WEADS_delft3d4.py
# Author: Shabnam
# Modified: June 2025

import os
import sys
import time 
import argparse

# Allow importing src_point if not already on PYTHONPATH
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_POINT_PATH = os.path.join(CURRENT_DIR, "src_point")
if SRC_POINT_PATH not in sys.path:
    sys.path.insert(0, SRC_POINT_PATH)

from src_point import basics_delft3d4
from src_point import preprocessing_delft3d4, mem_delft3d4, postprocessing_delft3d4
from src_point.tidaldatums_delft3d4 import calculate_tidal_metrics_from_csv



startTime = time.time()

parser = argparse.ArgumentParser(description="Run WEADS-Delft3D grid-based coupling.")

# Main processing flags
parser.add_argument("--all", action="store_true", help="Run all steps")
parser.add_argument("--preprocessing", action="store_true", help="Run preprocessing")
parser.add_argument("--mem", action="store_true", help="Run MEM step")
parser.add_argument("--postprocessing", action="store_true", help="Run postprocessing step")

# Required input files
parser.add_argument("--inputGrdFile", type=str, help="Delft3D grid file (*.grd)")
parser.add_argument("--inputEncFile", type=str, help="Delft3D enclosure file (*.enc)")
parser.add_argument("--inputDepFile", type=str, help="Delft3D depth file (*.dep)")
parser.add_argument("--updatedRghFile", type=str, default="updated.rgh", help="Updated .rgh file name")
parser.add_argument("--inputMdfFile", type=str, help="Delft3D mdf config file (*.mdf)")

# Optional input
parser.add_argument("--inputvegetationFile", type=str, default=None, help="Vegetation raster file (*.tif)")

# Output
parser.add_argument("--outputMEMFile", type=str, default="ecology.csv", help="Output MEM CSV file")
parser.add_argument("--updatedDepFile", type=str, default="updated.dep", help="Updated .dep file name")

# Projection info
parser.add_argument("--inEPSG", type=int, required=True, help="Input EPSG")
parser.add_argument("--outEPSG", type=int, required=True, help="Output EPSG")

# Simulation setup
parser.add_argument("--deltaT", type=float, default=5, help="Delta time in years")

args = parser.parse_args()

# Step control
if args.all:
    args.preprocessing = True
    args.mem = True
    args.postprocessing = True

print("\n#################################################")
print("Running WEADS-Delft3D Coupling")
print("#################################################\n")

# Step 1: Extract tidal datums from CSV
print("Calculating tidal datums from extracted CSV...")
tidal_csv = "tidal_metrics.csv"
calculate_tidal_metrics_from_csv('water_level_extracted.csv', output_csv=tidal_csv)

# Step 2: Preprocessing step
if args.preprocessing:
    print("Running Preprocessing...")
    preprocessing_delft3d4.preprocessing_Delft3D(
        inputBathymetryFile=args.inputDepFile,
        inputWaterLevelCSV=tidal_csv,
        inputShapeFile=None,
        domainIOFile=args.outputMEMFile,
        inEPSG=args.inEPSG,
        outEPSG=args.outEPSG
    )

# Step 3: MEM step
if args.mem:
    print("Running MEM Step...")
    mem_delft3d4.mem_delft(
        domainIOFile=args.outputMEMFile,
        vegetationFile=args.inputvegetationFile,
        outputMEMFile=args.outputMEMFile,
        inEPSG=args.inEPSG,
        outEPSG=args.outEPSG,
        deltaT=args.deltaT
    )

# Step 4: Postprocessing step
if args.postprocessing:
    print("Running Postprocessing...")
    postprocessing_delft3d4.postprocessing_delft(
        inputDepFile=args.inputDepFile,
        inputMdfFile=args.inputMdfFile,
        outputMEMFile=args.outputMEMFile,
        dep_out=args.updatedDepFile,
        rgh_out=args.updatedRghFile
    )

print("\n#################################################")
print("WEADS-Delft3D Coupling Complete!")
print("Elapsed Time: %.2f seconds" % (time.time() - startTime))
print("#################################################\n")
