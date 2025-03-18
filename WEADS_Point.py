#!/usr/bin/env python
# File: WEADS_Point.py
# Modified: March 18, 2025 by Jin Ikeda

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------

import argparse
import src_point
from src_point.general_functions import *

# ----------------------------------------------------------
# F U N C T I O N    M A I N
# ----------------------------------------------------------
# Example:
# python WEADS_Point.py --inputMeshFile fort.14 --inputAttrFile fort.13
# --inputHarmonicsFile fort.53 --inputShapeFile Cut_domain.shp --inEPSG 4269 --outEPSG 26914
# --deltaT 25 --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14
# --inputInundationtimeFile inundationtime.63 --inputvegetationFile NWI_TX_wetlands.tif --all
# ----------------------------------------------------------


startTime = time.time()

# ----------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------

# Create an ArgumentParser instance
parser = argparse.ArgumentParser(description="Running point-based WEADS simulation.")

tstep = 1800.0  # Calculate time steps for tidal datum [seconds]
numIDWNeighbors = 12  # For idw interpolations

# Add arguments
parser.add_argument(
    "--all",
    action="store_true",
    help="Run all processing steps")
parser.add_argument(
    "--preprocessing",
    action="store_true",
    help="Run preprocessing step")
parser.add_argument(
    "--td",
    action="store_true",
    help="Run tidal datum step")
parser.add_argument("--mem", action="store_true", help="Run mem step")
parser.add_argument(
    "--postprocessing",
    action="store_true",
    help="Run postprocessing step")
parser.add_argument(
    "--inputMeshFile",
    type=str,
    help="Input mesh <fort.14 file>")
parser.add_argument(
    "--inputAttrFile",
    type=str,
    help="Input attribute <fort.13 file>")
parser.add_argument(
    "--inputHarmonicsFile",
    type=str,
    help="Input harmonics <fort.53 file>")
parser.add_argument(
    "--inputShapeFile",
    type=str,
    help="Input domain shape file <*.shp>")
parser.add_argument(
    "--outputMEMFile",
    type=str,
    default="ecology",
    help="Output MEM file:ecology.csv")
parser.add_argument(
    "--outputMeshFile",
    type=str,
    default="fort_new.14",
    help="Output mesh file <fort_new.14>")
parser.add_argument(
    "--outputAttrFile",
    type=str,
    default="fort_new.13",
    help="Output attribute file")
parser.add_argument(
    "--inEPSG",
    type=str,
    help="Input EPSG code <inEPSGCode>")
parser.add_argument(
    "--outEPSG",
    type=str,
    help="Output EPSG code <outEPSGCode>")
parser.add_argument(
    "--deltaT",
    type=float,
    default=5,
    help="time step for WEADS simulation in years")
parser.add_argument("--slr", type=float, help="Sea level rise")
parser.add_argument(
    "--inputInundationtimeFile",
    type=str,
    help="Use inundationtime file for running inunT <inundationtime.63 & .63.nc>")
parser.add_argument(
    "--inputMaxinundationdepthFile",
    type=str,
    default=None,
    required=False,
    help="Use ADCIRC MAx inundation depth file for plotting (e.g., maxinundepth.63 or maxinundepth.63.nc)"
)
parser.add_argument(
    "--inputvegetationFile",
    type=str,
    default=None,
    help="Path to vegetation file <*.tif>")
parser.add_argument(
    "--skip_extracting_raster",
    action="store_true",
    help="Skip extract original raster values")
parser.add_argument(
    "--no_spread_flag",
    action="store_false",
    help="Turn off spread flag for vegetation mapping")

# Parse the command line arguments
args = parser.parse_args()

# Retrieve the arguments
all_flag = args.all
preprocessing_flag = args.preprocessing
td_flag = args.td
mem_flag = args.mem
postprocessing_flag = args.postprocessing

inputMeshFile = args.inputMeshFile
inputAttrFile = args.inputAttrFile
inputHarmonicsFile = args.inputHarmonicsFile
inputInundationTFile = args.inputInundationtimeFile
inputInundationMaxDFile = args.inputMaxinundationdepthFile
inputShapeFile = args.inputShapeFile
outputMEMFile = args.outputMEMFile
outputMeshFile = args.outputMeshFile
outputAttrFile = args.outputAttrFile
inEPSG = args.inEPSG
outEPSG = args.outEPSG
deltaT = args.deltaT
slr = args.slr
inputvegetationFile = args.inputvegetationFile
skip_extracting_raster_Flag = args.skip_extracting_raster
spread_flag = args.no_spread_flag

inEPSG = int(inEPSG)
outEPSG = int(outEPSG)
deltaT = float(deltaT)
slr = float(slr)

# ----------------------------------------------------------

print('\n')
print('#################################################')
print('Starting Point-based WEADS Simulation')
print('#################################################')
print('\n')

# ----------------------------------------------------------
# Function calls
# ----------------------------------------------------------

domainIOFile = "domain_inputs.csv"
inputHarmonicFreqFile = "harmonics_Freq.txt"
outputHarmonicsFile = "domain_tide.csv"
interpolateHarmonicsFile = "tidal_prj.csv"
outputvegetationFile = 'domain_nwi.csv'

if all_flag:
    preprocessing_flag = True
    td_flag = True
    mem_flag = True
    postprocessing_flag = True

if preprocessing_flag:  # Read hydrodynamics
    print('\nReading input hydrodynamics...')

    src_point.basics.fileexists(inputMeshFile)
    src_point.basics.fileexists(inputShapeFile)
    src_point.basics.fileexists(inputHarmonicsFile)
    src_point.basics.fileexists(inputInundationTFile)

    src_point.preprocessing_ADCIRC(
        inputMeshFile,
        inputAttrFile,
        inputInundationTFile,
        inputHarmonicsFile,
        inputShapeFile,
        domainIOFile,
        inEPSG,
        outEPSG)

if td_flag:  # Compute tidal datums

    print('\n' + '\tComputing tidal datums...')
    src_point.basics.fileexists(inputHarmonicFreqFile)

    src_point.tidaldatums(
        domainIOFile,
        inputHarmonicFreqFile,
        outputHarmonicsFile,
        tstep)
    src_point.tidaldatumsidw(
        outputHarmonicsFile,
        interpolateHarmonicsFile,
        inEPSG,
        outEPSG,
        numIDWNeighbors)

if mem_flag:  # Run MEM
    print('\n' + '\tRunning MEM...')

    src_point.basics.fileexists(interpolateHarmonicsFile)

    if inputvegetationFile:  # Organize vegetation file
        print('\n' + '\tOrganizing vegetation file...')
        src_point.basics.fileexists(inputvegetationFile)
        # First simulation # --skip_extracting_raster_Flag = False
        # After first simulation # skip_extracting_raster_Flag = True
        #spread_flag = True

        src_point.vegetation.process_vegetation_file(
            inputvegetationFile,
            skip_extracting_raster_Flag,
            spread_flag,
            outputvegetationFile,
            inEPSG,
            outEPSG,
            deltaT=deltaT)

        print('\n' + '\tUse vegetation mapping...')
        src_point.mem(
            interpolateHarmonicsFile,
            outputvegetationFile,
            outputMEMFile + '.csv',
            inEPSG,
            outEPSG,
            deltaT=deltaT)

    else:  # Run MEM without vegetation
        print('\n' + '\tNo vegetation mapping references...')
        src_point.mem(
            interpolateHarmonicsFile,
            inputvegetationFile,
            outputMEMFile + '.csv',
            inEPSG,
            outEPSG,
            deltaT=deltaT)

if postprocessing_flag:
    src_point.postprocessing_ADCIRC(
        inputMeshFile,
        inputAttrFile,
        outputMeshFile,
        outputAttrFile,
        outputMEMFile + '.csv',
        slr,
        inputShapeFile,
        inEPSG,
        outEPSG,
        raster_resolution=100,
        inputInundationMaxDFile=inputInundationMaxDFile)

    print('Finished new fort.14, new fort.13 and rasterization')

print('\n' + '#################################################')
print('Point-based WEADS Complete!')
print("--- %s seconds ---" % (time.time() - startTime))
print('#################################################\n')



