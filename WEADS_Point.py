#!/usr/bin/python3
# File: WEADS_point.py
# Modified: July 15, 2024 by Jin Ikeda

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------

import argparse
import src_point
from src_point.general_functions import *
#----------------------------------------------------------
# F U N C T I O N    M A I N
#----------------------------------------------------------
# Example:
#python WEADS_point.py --inputMeshFile fort.14
#--inputAttrFile fort.13 --inputHarmonicsFile fort.53
#--inputShapeFile GB_bbox.shp
#--inEPSG 4269 --outEPSG 26916 --gridSize 20 --outputMEMRasterFile mem
#--slr 1.1 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --all
#----------------------------------------------------------
def main(argv):
    
    startTime = time.time()
    
    #----------------------------------------------------------
    # Command line arguments
    #----------------------------------------------------------

    # Create an ArgumentParser instance
    parser = argparse.ArgumentParser(description="Running point-based WEADS simulation.")

    rts = 3600.0
    numIDWNeighbors = 12

    # Add arguments
    parser.add_argument("--all", action="store_true", help="Run all processing steps")
    parser.add_argument("--preprocessing", action="store_true", help="Run preprocessing step")
    parser.add_argument("--td", action="store_true", help="Run tidal datum step")
    parser.add_argument("--mem", action="store_true", help="Run mem step")
    parser.add_argument("--postprocessing", action="store_true", help="Run postprocessing step")
    parser.add_argument("--inputMeshFile", type=str, help="Input mesh <fort.14 file>")
    parser.add_argument("--inputAttrFile", type=str, help="Input attribute <fort.13 file>")
    parser.add_argument("--inputHarmonicsFile", type=str, help="Input harmonics <fort.53 file>")
    parser.add_argument("--inputShapeFile", type=str, help="Input domain shape file <*.shp>")
    parser.add_argument("--outputMEMFile", type=str, help="Output MEM file:mem.csv")
    parser.add_argument("--outputMeshFile", type=str, default="fort_new.14", help="Output mesh file <fort_new.14>")
    parser.add_argument("--outputAttrFile", type=str, default="fort_new.13", help="Output attribute file")
    parser.add_argument("--inEPSG", type=str, help="Input EPSG code <inEPSGCode>")
    parser.add_argument("--outEPSG", type=str, help="Output EPSG code <outEPSGCode>")
    parser.add_argument("--deltaT", type=float,default=5, help="time step for WEADS simuiation in years")
    parser.add_argument("--slr", type=float, help="Sea level rise")
    parser.add_argument("--inputInundationtimeFile", type=str, help="Use inundationtime file for running inunT <inundationtime.63>")
    parser.add_argument("--InputvegetationFile", type=str,default=None, help="Path to vegetation file <*.tif>")
    parser.add_argument("--outputVegetationFile", type=str, help="Output MEM file:mem.csv")
    parser.add_argument("--skipresample", action="store_true", help="Skip reprojection and resample to raster domain")

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
    inputShapeFile = args.inputShapeFile
    outputMEMFile = args.outputMEMFile
    outputMeshFile = args.outputMeshFile
    outputAttrFile = args.outputAttrFile
    inEPSG = args.inEPSG
    outEPSG = args.outEPSG
    deltaT = args.deltaT
    slr = args.slr
    InputvegetationFile = args.InputvegetationFile


    inEPSG = int(inEPSG)
    outEPSG = int(outEPSG)
    deltaT = float(deltaT)
    slr = float(slr)

    #----------------------------------------------------------
    
    print('\n' + '#################################################')
    print('#################################################')
    print('Starting Point-based WEADS Simulation')
    print('#################################################')
    
    #----------------------------------------------------------
    # Function calls
    #----------------------------------------------------------

    domainIOFile = "domain_inputs.csv"
    inputHarmonicFreqFile = "harmonics_Freq.txt"
    outputHarmonicsFile = "domain_tide.csv"
    interpolateHarmonicsFile = "tidal_prj.csv"

    if all_flag:
        preprocessing_flag = True
        td_flag = True
        mem_flag = True
        postprocessing_flag = True

    if preprocessing_flag: # Read hydrodynamics
        print('\n' + '\tReading input hydrodynamics...')
        
        src_point.basics.fileexists(inputMeshFile)
        src_point.basics.fileexists(inputShapeFile)
        src_point.basics.fileexists(inputHarmonicsFile)
        src_point.basics.fileexists(inputInundationTFile)

        src_point.preprocessing_ADCIRC(inputMeshFile, inputAttrFile, inputInundationTFile, inputHarmonicsFile, inputShapeFile, domainIOFile, inEPSG, outEPSG)

    if td_flag: # Compute tidal datums

        print('\n' + '\tComputing tidal datums...')
        src_point.basics.fileexists(inputHarmonicFreqFile)

        src_point.tidaldatums(domainIOFile,inputHarmonicFreqFile,outputHarmonicsFile,rts)
        src_point.tidaldatumsidw(outputHarmonicsFile,interpolateHarmonicsFile,inEPSG,outEPSG, numIDWNeighbors)

    if mem_flag: # Run MEM
        print('\n' + '\tRunning MEM...')

        src_point.basics.fileexists(interpolateHarmonicsFile)

        if InputvegetationFile:  # Organize vegetation file
            print('\n' + '\tOrganizing vegetation file...')
            src_point.basics.fileexists(InputvegetationFile)
            # First simulation
            skip_raster_extracting_Flag = False
            # After first simulation
            #skip_raster_extracting_Flag = True
            spread_flag = True
            OutputvegetationFile = 'domain_nwi.csv'

            src_point.process_vegetation_file(InputvegetationFile,skip_raster_extracting_Flag,spread_flag,OutputvegetationFile,inEPSG,outEPSG,deltaT=deltaT)

            print('\n' + '\tUse vegetation mapping...')
            src_point.mem(interpolateHarmonicsFile, OutputvegetationFile, outputMEMFile + '.csv',inEPSG,outEPSG, deltaT=deltaT)

        else:  # Run MEM without vegetation
            print('\n' + '\tNo vegetation mapping references...')
            src_point.mem(interpolateHarmonicsFile, InputvegetationFile, outputMEMFile + '.csv',inEPSG,outEPSG, deltaT=deltaT)

    if postprocessing_flag:
        src_point.postprocessing_ADCIRC(inputMeshFile,inputAttrFile,outputMeshFile,outputAttrFile,outputMEMFile + '.csv',slr)  # 3 is the annual accreation rate on MEM (m/yr). Caution: -1 is the multiplier (for ADCIRC file)
        print('Finished new fort.14 and new fort.13')

    print('\n' + '#################################################')
    print('Point-based WEADS Complete!')
    print("--- %s seconds ---" % (time.time() - startTime))
    print('#################################################\n')

#----------------------------------------------------------

if __name__ == "__main__":
    main((sys.argv[1:]))