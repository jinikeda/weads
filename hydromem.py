#!/usr/bin/python3
# File: hydromem.py
# Date: Aug 15, 2023

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import os
import argparse
import src
import sys
import time

#----------------------------------------------------------
# F U N C T I O N    M A I N
#----------------------------------------------------------
# Example:
#python hydromem.py --inputMeshFile fort.14
#--inputAttrFile fort.13 --inputHarmonicsFile fort.53
#--inputEverdriedFile everdried.63 --inputShapeFile GB_bbox.shp
#--inEPSG 4269 --outEPSG 26916 --gridSize 20 --outputMEMRasterFile mem
#--slr 1.1 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --all
#----------------------------------------------------------
def main(argv):
    
    startTime = time.time()
    
    #----------------------------------------------------------
    # Command line arguments
    #----------------------------------------------------------

    # Create an ArgumentParser instance
    parser = argparse.ArgumentParser(description="Running WEAD")

    rts = 900 #rts = 3600
    numIDWNeighbors = 12

    # Add arguments
    parser.add_argument("--all", action="store_true", help="Run all processing steps")
    parser.add_argument("--rasterize", action="store_true", help="Run rasterization step")
    parser.add_argument("--hyconn", action="store_true", help="Run hyconn step")
    parser.add_argument("--td", action="store_true", help="Run tidal datum step")
    parser.add_argument("--mem", action="store_true", help="Run mem step")
    parser.add_argument("--adc2rast", action="store_true", help="Run adc2rast step")
    parser.add_argument("--inputMeshFile", type=str, help="Input mesh <fort.14 file>")
    parser.add_argument("--inputAttrFile", type=str, help="Input attribute <fort.13 file>")
    parser.add_argument("--inputHarmonicsFile", type=str, help="Input harmonics <fort.53 file>")
    parser.add_argument("--inputEverdriedFile", type=str, help="Input everdried.63 file")
    parser.add_argument("--inputShapeFile", type=str, help="Input domain shape file <*.shp>")
    parser.add_argument("--outputMEMRasterFile", type=str, help="Output MEM raster file:HydroMEM.tif")
    parser.add_argument("--outputMeshFile", type=str, help="Output mesh file <fort_new.14>")
    parser.add_argument("--outputAttrFile", type=str, help="Output attribute file")
    parser.add_argument("--inEPSG", type=str, help="Input EPSG code <inEPSGCode>")
    parser.add_argument("--outEPSG", type=str, help="Output EPSG code <outEPSGCode>")
    parser.add_argument("--gridSize", type=float, help="Output raster resolution")
    parser.add_argument("--slr", type=float, help="Sea level rise")
    parser.add_argument("--inundationFile", type=str, help="Use inundation file <maxinundepth.63>")
    parser.add_argument("--vegetationFile", type=str, help="Path to vegetation file <*.tif>")
    parser.add_argument("--skipresample", action="store_true", help="Skip reprojection and resample to raster domain")

    # Parse the command line arguments
    args = parser.parse_args()

    # Retrieve the arguments
    all_flag = args.all
    rasterize_flag = args.rasterize
    hyconn_flag = args.hyconn
    td_flag = args.td
    mem_flag = args.mem
    adc2rast_flag = args.adc2rast

    inputMeshFile = args.inputMeshFile
    inputAttrFile = args.inputAttrFile
    inputHarmonicsFile = args.inputHarmonicsFile
    inputEverdriedFile = args.inputEverdriedFile
    inputShapeFile = args.inputShapeFile
    outputMEMRasterFile = args.outputMEMRasterFile
    outputMeshFile = args.outputMeshFile
    outputAttrFile = args.outputAttrFile
    inEPSG = args.inEPSG
    outEPSG = args.outEPSG
    gridSize = args.gridSize
    slr = args.slr
    inundationFile = args.inundationFile
    vegetationFile = args.vegetationFile
    skipresample_flag = args.skipresample

    inEPSG = int(inEPSG)
    outEPSG = int(outEPSG)
    gridSize = float(gridSize)
    slr = float(slr)

    #----------------------------------------------------------
    
    print('\n' + '#################################################')
    print('#################################################')
    print('Starting pyHydroMEM')
    print('#################################################')
    
    #----------------------------------------------------------
    # Function calls
    #----------------------------------------------------------

    if all_flag:
        rasterize_flag = True
        hyconn_flag = True
        td_flag = True
        mem_flag = True
        adc2rast_flag = True

    if rasterize_flag: # Create TIF images
        print('\n' + '\tCreating TIF images...')
        
        # src.basics.fileexists(inputMeshFile)
        # src.basics.fileexists(inputShapeFile)
        # src.basics.fileexists(inputEverdriedFile)
        # src.basics.fileexists(inputHarmonicsFile)
        #
        # src.grd2dem(inputMeshFile,inputMeshFile,inputShapeFile,'tbathy',inEPSG,outEPSG,gridSize,-1)
        # src.grd2dem(inputMeshFile,inputEverdriedFile,inputShapeFile,'everdried',inEPSG,outEPSG,gridSize,1,1,True)
        # src.grd2dem(inputMeshFile,inputAttrFile,inputShapeFile,'manning',inEPSG,outEPSG,gridSize,1)
        # src.grd2dem(inputMeshFile,inputHarmonicsFile,inputShapeFile,'harmonics',inEPSG,outEPSG,gridSize,1)

        if inundationFile:  # Run Inundation calculation
            print('\n' + '\tMake a raster of maximum inundation depth...')
            src.basics.fileexists(inundationFile)
            src.grd2dem(inputMeshFile,inundationFile,inputShapeFile, 'maxinundation', inEPSG, outEPSG, gridSize,1,1,False)
    
    if hyconn_flag: # Create TIF of hydraulically connected area
        print('\n' + '\tComputing hydraulic connectivity...')
        
        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('everdried.tif')
        
        src.hyconn('tbathy.tif','everdried.tif','hyconn.tif',False)
    
    if td_flag: # Compute tidal datums
        
        print('\n' + '\tComputing tidal datums...')
        
        src.basics.fileexists('harmonics.tif')
        src.basics.fileexists('hyconn.tif')
        
        src.tidaldatums('harmonics.README','harmonics.tif','hyconn.tif','TidalDatums.tif',rts)
        
        src.basics.fileexists('TidalDatums.tif')
        
        src.tidaldatumsidw('hyconn.tif','TidalDatums.tif','TidalDatums_IDW.tif',numIDWNeighbors)
    
    if vegetationFile: # Organize vegetation file
        print('\n' + '\tOrganizing vegetation file...')
        src.basics.fileexists(vegetationFile)
        src.nwi.read_tif(vegetationFile, outEPSG, gridSize, skipresample_flag)

    if mem_flag: # Run MEM
        print('\n' + '\tRunning MEM...')

        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('hyconn.tif')
        src.basics.fileexists('TidalDatums_IDW.tif')

        if vegetationFile is None: # Run MEM without vegetation
            print('\n' + '\tNo vegetation mapping references...')
            src.mem('hyconn.tif', 'tbathy.tif', 'TidalDatums_IDW.tif', vegetationFile, outputMEMRasterFile + '.tif')
        else:
            print('\n' + '\tUse vegetation mapping...')
            Domain_raster = 'Domain_classification_distribution_resample100.tif'
            if not os.path.isfile(Domain_raster):
                print("Could not find " + Domain_raster)
                sys.exit(1)
            else:
                src.mem('hyconn.tif', 'tbathy.tif', 'TidalDatums_IDW.tif',Domain_raster, outputMEMRasterFile + '.tif')

    if adc2rast_flag:
        src.rast2adc(inputMeshFile,outputMeshFile,outputMEMRasterFile+'.tif',inEPSG,4,1)
        print('Finished new fort.14')
        src.update_nodal_attributes(inputMeshFile,outputMEMRasterFile+'.tif',
                           inputAttrFile,outputAttrFile,
                           inEPSG,slr,6)
        print('Finished new fort.13')

    print('\n' + '#################################################')
    print('pyHydro-MEM Complete!')
    print("--- %s seconds ---" % (time.time() - startTime))
    print('#################################################\n')

#----------------------------------------------------------

if __name__ == "__main__":
    main((sys.argv[1:]))
