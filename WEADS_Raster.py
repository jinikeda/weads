#!/usr/bin/python3
# File: WEADS_Raster.py
# Date: Aug 22, 2024
# Modified: Aug 22, 2024 by Jin Ikeda

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------

import os
import argparse
import src
import sys
import time

# ----------------------------------------------------------
# F U N C T I O N    M A I N
# ----------------------------------------------------------
# Example:
# python WEAD_Raster.py --inputMeshFile fort.14
# --inputAttrFile fort.13 --inputHarmonicsFile fort.53
# --inputEverdriedFile everdried.63 --inputShapeFile GB_bbox.shp
# --inEPSG 4269 --outEPSG 26916 --gridSize 20 --outputMEMRasterFile mem
# --slr 1.1 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --all
# ----------------------------------------------------------


def main(argv):

    startTime = time.time()

    # ----------------------------------------------------------
    # Command line arguments
    # ----------------------------------------------------------

    # Create an ArgumentParser instance
    parser = argparse.ArgumentParser(description="Running WEAD simulation.")

    tstep = 1800.0 # run time steps for tidal datum [seconds]
    numIDWNeighbors = 12 # For idw interpolations

    # Add arguments
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all processing steps")
    parser.add_argument(
        "--rasterize",
        action="store_true",
        help="Run rasterization step")
    parser.add_argument(
        "--hyconn",
        action="store_true",
        help="Run hyconn step")
    parser.add_argument(
        "--inunT",
        action="store_true",
        help="Run hyconn (subtidal, intertidal and land classification) step")
    parser.add_argument(
        "--td",
        action="store_true",
        help="Run tidal datum step")
    parser.add_argument("--mem", action="store_true", help="Run mem step")
    parser.add_argument(
        "--adc2rast",
        action="store_true",
        help="Run adc2rast step")
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
        "--inputEverdriedFile",
        type=str,
        help="Input everdried.63 file")
    parser.add_argument(
        "--inputShapeFile",
        type=str,
        help="Input domain shape file <*.shp>")
    parser.add_argument(
        "--outputMEMRasterFile",
        type=str,
        default="ecology",
        help="Output MEM raster file:ecology.tif")
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
        "--gridSize",
        type=float,
        help="Output raster resolution")
    parser.add_argument(
        "--deltaT",
        type=float,
        default=5,
        help="time step for WEADS simuiation in years")
    parser.add_argument("--slr", type=float, help="Sea level rise")
    parser.add_argument(
        "--inputInundationtimeFile",
        type=str,
        help="Use inundationtime file for running inunT <inundationtime.63>")
    parser.add_argument(
        "--inundationdepthFile",
        type=str,
        help="Use inundation depth file <maxinundepth.63>")
    parser.add_argument(
        "--InputvegetationFile",
        type=str,
        default=None,
        help="Path to vegetation file <*.tif>")
    parser.add_argument(
        "--skipresample",
        action="store_true",
        help="Skip reprojection and resample to raster domain")

    # Parse the command line arguments
    args = parser.parse_args()

    # Retrieve the arguments
    all_flag = args.all
    rasterize_flag = args.rasterize
    hyconn_flag = args.hyconn
    inunT_flag = args.inunT
    td_flag = args.td
    mem_flag = args.mem
    adc2rast_flag = args.adc2rast

    inputMeshFile = args.inputMeshFile
    inputAttrFile = args.inputAttrFile
    inputHarmonicsFile = args.inputHarmonicsFile
    inputEverdriedFile = args.inputEverdriedFile
    inputInundationTFile = args.inputInundationtimeFile
    inputShapeFile = args.inputShapeFile
    outputMEMRasterFile = args.outputMEMRasterFile
    outputMeshFile = args.outputMeshFile
    outputAttrFile = args.outputAttrFile
    inEPSG = args.inEPSG
    outEPSG = args.outEPSG
    gridSize = args.gridSize
    deltaT = args.deltaT
    slr = args.slr
    inundationdepthFile = args.inundationdepthFile
    vegetationFile = args.InputvegetationFile
    skipresample_flag = args.skipresample

    inEPSG = int(inEPSG)
    outEPSG = int(outEPSG)
    gridSize = float(gridSize)
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

    if all_flag:
        rasterize_flag = True
        hyconn_flag = False  # replaced by inunT_flag
        inunT_flag = True
        td_flag = True
        mem_flag = True
        adc2rast_flag = True

    if rasterize_flag:  # Create TIF images
        print('\n' + '\tCreating TIF images...')

        src.basics.fileexists(inputMeshFile)
        src.basics.fileexists(inputShapeFile)
        src.basics.fileexists(inputEverdriedFile)
        src.basics.fileexists(inputHarmonicsFile)
        src.basics.fileexists(inputInundationTFile)

        src.grd2dem(inputMeshFile, inputMeshFile, inputShapeFile,
                    'tbathy', inEPSG, outEPSG, gridSize, -1)
        # src.grd2dem(inputMeshFile,inputEverdriedFile,inputShapeFile,'everdried',inEPSG,outEPSG,gridSize,1,1,True)
        src.grd2dem(
            inputMeshFile,
            inputAttrFile,
            inputShapeFile,
            'manning',
            inEPSG,
            outEPSG,
            gridSize,
            1)
        src.grd2dem(
            inputMeshFile,
            inputHarmonicsFile,
            inputShapeFile,
            'harmonics',
            inEPSG,
            outEPSG,
            gridSize,
            1)
        src.grd2dem(
            inputMeshFile,
            inputInundationTFile,
            inputShapeFile,
            'inundationtime',
            inEPSG,
            outEPSG,
            gridSize,
            1,
            1,
            False)

        if inundationdepthFile:  # Run Inundation calculation
            print('\n' + '\tMake a raster of maximum inundation depth...')
            src.basics.fileexists(inundationdepthFile)
            src.grd2dem(
                inputMeshFile,
                inundationdepthFile,
                inputShapeFile,
                'maxinundation',
                inEPSG,
                outEPSG,
                gridSize,
                1,
                1,
                False)

    if hyconn_flag:  # Create TIF of hydraulically connected area
        print('\n' + '\tComputing hydraulic connectivity...')

        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('everdried.tif')

        src.hyconn('tbathy.tif', 'everdried.tif', 'hyconn.tif', False)

    if inunT_flag:  # Create TIF of hydraulically connected area
        print(
            '\n' +
            '\tComputing land classifcation and thier hydraulic connectivity...')

        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('inundationtime.tif')

        src.hydro_classify(
            'tbathy.tif',
            'inundationtime.tif',
            'hydro_class.tif',
            False)

    if td_flag:  # Compute tidal datums

        print('\n' + '\tComputing tidal datums...')

        src.basics.fileexists('harmonics.tif')
        src.basics.fileexists('hydro_class.tif')

        src.tidaldatums(
            'harmonics.README',
            'harmonics.tif',
            'hydro_class.tif',
            'TidalDatums.tif',
            tstep)

        src.basics.fileexists('TidalDatums.tif')

        src.tidaldatumsidw(
            'hydro_class.tif',
            'TidalDatums.tif',
            'TidalDatums_IDW.tif',
            numIDWNeighbors)

    if mem_flag:  # Run MEM
        print('\n' + '\tRunning MEM...')

        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('hydro_class.tif')
        src.basics.fileexists('TidalDatums_IDW.tif')

        if vegetationFile:  # Organize vegetation file
            print('\n' + '\tOrganizing vegetation file...')
            src.basics.fileexists(vegetationFile)
            src.nwi.read_tif(
                vegetationFile,
                outEPSG,
                gridSize,
                skipresample_flag,
                deltaT=deltaT)

            print('\n' + '\tUse vegetation mapping...')
            # Domain raster be careful
            Domain_raster = f'Domain_classification_distribution_resample{int(gridSize)}.tif'
            if not os.path.isfile(Domain_raster):
                print("Could not find " + Domain_raster)
                sys.exit(1)
            else:
                src.mem('hydro_class.tif', 'tbathy.tif', 'TidalDatums_IDW.tif', Domain_raster,
                        outputMEMRasterFile + '.tif', deltaT=deltaT)

        else:  # Run MEM without vegetation
            print('\n' + '\tNo vegetation mapping references...')
            src.mem(
                'hydro_class.tif',
                'tbathy.tif',
                'TidalDatums_IDW.tif',
                vegetationFile,
                outputMEMRasterFile +
                '.tif',
                deltaT=deltaT)

    if adc2rast_flag:
        # 3 is the annual accreation rate on MEM (m/yr). Caution: -1 is the
        # multiplier (for ADCIRC file)
        src.rast2adc(inputMeshFile, outputMeshFile,
                     outputMEMRasterFile + '.tif', inEPSG, 3, -deltaT)
        print('Finished new fort.14')
        src.update_nodal_attributes(inputMeshFile, outputMEMRasterFile + '.tif',
                                    inputAttrFile, outputAttrFile,
                                    inEPSG, slr, 6)
        print('Finished new fort.13')

    print('\n' + '#################################################')
    print('WEADSM Complete!')
    print("--- %s seconds ---" % (time.time() - startTime))
    print('#################################################\n')

# ----------------------------------------------------------


if __name__ == "__main__":
    main((sys.argv[1:]))
