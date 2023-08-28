#!/usr/bin/python3
# File: hydromem.py
# Date: Aug 15, 2023

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import src
import getopt, sys
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
    inputMeshFile = ''
    inputAttrFile = ''
    inputHarmonicsFile = ''
    inputEverdriedFile = ''
    inputShapeFile = ''
    outputMEMRasterFile = ''
    outputMeshFile = ''
    outputAttrFile = ''
    inEPSG = ''
    outEPSG = ''
    gridSize = ''
    slr = ''
    dorasterize = False
    dohyconn = False
    dotd = False
    domem = False
    dofinal = False
    inundationFile = False
    vegetationFile = None
    
    rts = 900
    #rts = 3600
    numIDWNeighbors = 12
    
    if len(argv) < 1:
        print('USAGE: hydromem.py --inputMeshFile <fort.14>'\
              '--inputAttrFile <fort.13>',
              '--inputHarmonicsFile <fort.53>',
              '--inputEverdriedFile <everdried.63>',
              '--inputShapeFile <*.shp>',
              '--outputMEMRasterFile <HydroMEM.tif>',
              '--outputMeshFile <outputMeshFile.14>',
              '--outputAttrFile <outputAttrFile.13>',
              '--inEPSG <inEPSGCode>',
              '--outEPSG <outEPSGCode>',
              '--gridSize <outputRasterResolution>',
              '--slr <sea level rise>','\n')
    try:
        opts, args = getopt.getopt(
            argv,'hxi:a:h:e:s:r:m:n:p:q:g:l:',
                ['all',
                 'rasterize',
                 'hyconn',
                 'td',
                 'mem',
                 'adc2rast',
                 'inputMeshFile=',
                 'inputAttrFile=',
                 'inputHarmonicsFile=',
                 'inputEverdriedFile=',
                 'inputShapeFile=',
                 'outputMEMRasterFile=',
                 'outputMeshFile=',
                 'outputAttrFile=',
                 'inEPSG=',
                 'outEPSG=',
                 'gridSize=',
                 'slr=',
                 'inundationFile=',
                 'vegetationFile='])
    except getopt.GetoptError as e:
        print(e)
        quit()
        print('Incorrect command line arguments.\n')
        print('grd2dem.py --inputMeshFile <fort.14>'\
              '-- inputAttrFile <fort.13>',
              '--inputHarmonicsFile <fort.53>',
              '--inputEverdriedFile <everdried.63>',
              '--inputShapeFile <*.shp>',
              '--outputMEMRasterFile <HydroMEM.tif>',
              '--outputMeshFile <fort_new.14>',
              '--outputAttrFile <outputAttrFile.13>',
              '--inEPSG <inEPSGCode>',
              '--outEPSG <outEPSGCode>',
              '--gridSize <outputRasterResolution>',
              '--slr <sea level rise>',
              '--inundationFile <maxele.63>',
              '--vegetationFile <x.tif>',
              '--<all; rasterize; hyconn; td; mem>\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('grd2dem.py --inputMeshFile <inputMeshFile>',
                    '--- inputAttrFile <fort.13>',
                    '--inputHarmonicsFile <fort.53>',
                    '--inputEverdriedFile <everdried.63>',
                    '--inputShapeFile <*.shp>',
                    '--outputMEMRasterFile <HydroMEM.tif>',
                    '--outputMeshFile <fort_new.14>',
                    '--outputAttrFile <outputAttrFile.13>',
                    '--inEPSG <inEPSGCode>',
                    '--outEPSG <outEPSGCode>',
                    '--gridSize <outputRasterResolution>',
                    '--slr',
                    '--inundationFile <maxele.63>.',
                    '--vegetationFile <x.tif>',
                    '--<all; rasterize; hyconn; td; mem>\n')
            sys.exit(2)
        elif opt in ('--all'):
            dorasterize = True
            dohyconn = True
            dotd = True
            domem = True
            dofinal = True
        elif opt in ('--rasterize'):
            dorasterize = True
        elif opt in ('--hyconn'):
            dohyconn = True
        elif opt in ('--td'):
            dotd = True
        elif opt in ('--mem'):
            domem = True
        elif opt in ('--adc2rast'):
            dofinal = True
        elif opt in ('--inputMeshFile'):
            inputMeshFile = arg
        elif opt in ('--inputAttrFile'):
            inputAttrFile= arg
        elif opt in ('--inputHarmonicsFile'):
            inputHarmonicsFile = arg
        elif opt in ('--inputEverdriedFile'):
            inputEverdriedFile = arg
        elif opt in ('--inputShapeFile'):
            inputShapeFile = arg
        elif opt in ('--outputMEMRasterFile'):
            outputMEMRasterFile = arg
        elif opt in ('--outputMeshFile'):
            outputMeshFile = arg
        elif opt in ('--outputAttrFile'):
            outputAttrFile = arg
        elif opt in ('--inEPSG'):
            inEPSG = arg
        elif opt in ('--outEPSG'):
            outEPSG = arg
        elif opt in ('--gridSize'):
            gridSize = arg
        elif opt in ('--slr'):
            slr = arg
        elif opt in ('--inundationFile'):
            inundationDepth = arg
        elif opt in ('--vegetationFile'):
            vegetationFile = arg
        
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
   
    if dorasterize: # Create TIF images
        print('\n' + '\tCreating TIF images...')
        
        src.basics.fileexists(inputMeshFile)
        src.basics.fileexists(inputShapeFile)
        src.basics.fileexists(inputEverdriedFile)
        src.basics.fileexists(inputHarmonicsFile)
        
        src.grd2dem(inputMeshFile,inputMeshFile,inputShapeFile,'tbathy',inEPSG,outEPSG,gridSize,-1)
        src.grd2dem(inputMeshFile,inputEverdriedFile,inputShapeFile,'everdried',inEPSG,outEPSG,gridSize,1,1,True)
        src.grd2dem(inputMeshFile,inputAttrFile,inputShapeFile,'manning',inEPSG,outEPSG,gridSize,1)
        src.grd2dem(inputMeshFile,inputHarmonicsFile,inputShapeFile,'harmonics',inEPSG,outEPSG,gridSize,1)
    
    if dohyconn: # Create TIF of hydraulically connected area
        print('\n' + '\tComputing hydraulic connectivity...')
        
        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('everdried.tif')
        
        src.hyconn('tbathy.tif','everdried.tif','hyconn.tif',False)
    
    if dotd: # Compute tidal datums  
        
        print('\n' + '\tComputing tidal datums...')
        
        src.basics.fileexists('harmonics.tif')
        src.basics.fileexists('hyconn.tif')
        
        src.tidaldatums('harmonics.README','harmonics.tif','hyconn.tif','TidalDatums.tif',rts)
        
        src.basics.fileexists('TidalDatums.tif')
        
        src.tidaldatumsidw('hyconn.tif','TidalDatums.tif','TidalDatums_IDW.tif',numIDWNeighbors)
    
    if inundationDepth: # Run Inundation calculation
        print('\n' + '\tCalculating maximum inundation depth...')

    if domem: # Run MEM
        print('\n' + '\tRunning MEM...')

        src.basics.fileexists('tbathy.tif')
        src.basics.fileexists('hyconn.tif')
        src.basics.fileexists('TidalDatums_IDW.tif')

        if vegetationFile is None: # Run MEM without vegetation
            print('\n' + '\tNo vegetation mapping references...')
        else:
            print('\n' + '\tUse vegetation mapping...')

        src.mem('hyconn.tif', 'tbathy.tif', 'TidalDatums_IDW.tif',vegetationFile, outputMEMRasterFile + '.tif')

    if dofinal:
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
