#!/usr/bin/env python
# File: WEADS_Point.py
# Modified: April 22, 2025 by Jin Ikeda

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------

import argparse
import time
import sys
# import src_point
from src_point import basics, preprocessing_ADCIRC, tidaldatums, tidaldatumsidw, mem, vegetation, postprocessing_ADCIRC, mesh_edit

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
    help="Input domain shape file in a projected coordinate <*.shp>")
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
    help="Use ADCIRC Max inundation depth file for plots (e.g., maxinundepth.63 or maxinundepth.63.nc)"
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
parser.add_argument(
    "--z_adjustFile",
    type=str,
    default=None,
    help="Path to vector file for z adjustment <*.shp> or <*.geojson> Future tiff support")
parser.add_argument(
    "--zv",
    type=float,
    default=None,
    required=False,
    help="Value of z adjustment [meters]. If not provided, the script will not perform z adjustment.")

# Parse the command line arguments
args = parser.parse_args()

# Retrieve the arguments
all_flag = args.all
preprocessing_flag = args.preprocessing
td_flag = args.td
mem_flag = args.mem
postprocessing_flag = args.postprocessing

inputadjustFile = args.z_adjustFile
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
if inputadjustFile:
    pass
else:
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

if inputadjustFile:
    print('\n' + '\tAdjusting z values in the mesh file...')
    try:
        if args.zv is None:
            print("Error: Please provide a value for z adjustment using the --zv_adjust argument.")
            sys.exit(1)
        basics.fileexists(inputadjustFile)
        basics.fileexists(inputMeshFile)
        # basics.fileexists(inputAttrFile)  # not need now

        mesh_edit(inputMeshFile, inputadjustFile, inEPSG, z_adjust=args.zv)

        print(f"""
            Z values adjusted successfully in the mesh file.

        #################################################
        Z-adjustment Complete!
        --- {time.time() - startTime:.2f} seconds ---
        #################################################
        """)

        # Exit early after performing only the z adjustment
        sys.exit(0)

    except Exception as e:
        print(f" Error during z adjustment: {e}")
        sys.exit(1)  # Exit with error code


if all_flag:
    preprocessing_flag = True
    td_flag = True
    mem_flag = True
    postprocessing_flag = True

if preprocessing_flag:
    print('\nReading input hydrodynamics...')

    # Validate file existence
    for f in [inputMeshFile, inputShapeFile, inputHarmonicsFile, inputInundationTFile]:
        basics.fileexists(f)

    kwargs = dict(
        inputMeshFile=inputMeshFile,
        inputAttrFile=inputAttrFile,
        inputInundationtimeFile=inputInundationTFile,
        inputHarmonicsFile=inputHarmonicsFile,
        inputShapeFile=inputShapeFile,
        domainIOFile=domainIOFile,
        inEPSG=inEPSG,
        outEPSG=outEPSG
    )

    if inputInundationMaxDFile:
        basics.fileexists(inputInundationMaxDFile)
        kwargs["inputInundationMaxDFile"] = inputInundationMaxDFile

    preprocessing_ADCIRC(**kwargs)

if td_flag:  # Compute tidal datums

    print('\n' + '\tComputing tidal datums...')
    basics.fileexists(inputHarmonicFreqFile)

    tidaldatums(
        domainIOFile,
        inputHarmonicFreqFile,
        outputHarmonicsFile,
        tstep)
    tidaldatumsidw(
        outputHarmonicsFile,
        interpolateHarmonicsFile,
        inEPSG,
        outEPSG,
        numIDWNeighbors)

if mem_flag:  # Run MEM
    print('\n' + '\tRunning MEM...')

    basics.fileexists(interpolateHarmonicsFile)

    if inputvegetationFile:  # Organize vegetation file
        print('\n' + '\tOrganizing vegetation file...')
        basics.fileexists(inputvegetationFile)
        # First simulation # --skip_extracting_raster_Flag = False
        # After first simulation # skip_extracting_raster_Flag = True
        #spread_flag = True

        vegetation.process_vegetation_file(
            inputvegetationFile,
            skip_extracting_raster_Flag,
            spread_flag,
            outputvegetationFile,
            inEPSG,
            outEPSG,
            deltaT=deltaT)

        print('\n' + '\tUse vegetation mapping...')
        mem(
            interpolateHarmonicsFile,
            outputvegetationFile,
            outputMEMFile + '.csv',
            inEPSG,
            outEPSG,
            deltaT=deltaT)

    else:  # Run MEM without vegetation
        print('\n' + '\tNo vegetation mapping references...')
        mem(
            interpolateHarmonicsFile,
            inputvegetationFile,
            outputMEMFile + '.csv',
            inEPSG,
            outEPSG,
            deltaT=deltaT)

if postprocessing_flag:
    print("\n\tPostprocessing" + (" with inundation max depth..." if inputInundationMaxDFile else " without inundation max depth..."))

    kwargs = dict(
        inputMeshFile=inputMeshFile,
        inputAttrFile=inputAttrFile,
        outputMeshFile=outputMeshFile,
        outputAttrFile=outputAttrFile,
        outputMEMFile=outputMEMFile + '.csv',
        slr=slr,
        inputShapeFile=inputShapeFile,
        inEPSG=inEPSG,
        outEPSG=outEPSG,
        raster_resolution=100
    )

    if inputInundationMaxDFile:
        kwargs["inputInundationMaxDFile"] = inputInundationMaxDFile

    postprocessing_ADCIRC(**kwargs)

    print('Finished new fort.14, new fort.13 and rasterization')

print('\n' + '#################################################')
print('Point-based WEADS Complete!')
print("--- %s seconds ---" % (time.time() - startTime))
print('#################################################\n')



