#!/usr/bin/python3
# File: grd2dem.py
# Name: Matthew V. Bilskie
# Date: June 11, 2021

#----------------------------------------------------------
# M O D U L E S
#----------------------------------------------------------
#----------------------------------------------------------
import os, sys
import shapefile
import pyadcircmodules
from pyproj import Transformer
import numpy as np
from osgeo import gdal
from osgeo import osr
import src.basics
import rasterio

#----------------------------------------------------------
# F U N C T I O N    G R D 2 D E M
#----------------------------------------------------------
#
# Converts an ADCIRC unstructured mesh file to a structured
# raster file.
# result = grd2dem((inputMeshFile,inputDataFile,inputShapefile,
                #outputRasterFile,inEPSG,outEPSG,
                #gridSize,multFac,
                #snap (Optional, Default=0),
                #everdried (Optiona, Default = False)):
#----------------------------------------------------------
def grd2dem(inputMeshFile,inputDataFile,inputShapefile,\
                outputRasterFile,inEPSG,outEPSG,\
                gridSize,multFac,\
                snap=0,everdried=False):

    #----------------------------------------------------------
    # Check that input files exist
    #----------------------------------------------------------
    src.basics.fileexists(inputMeshFile)
    src.basics.fileexists(inputShapefile)
    #----------------------------------------------------------

    #----------------------------------------------------------
    # Set up in/out projection
    #----------------------------------------------------------
    transformer = Transformer.from_crs(inEPSG,outEPSG)
    #----------------------------------------------------------

    #----------------------------------------------------------
    # Read in shapefile, compute bounding box, and compute
    # the number of rows and columns for output raster given
    # the output raster grid cell size.
    #----------------------------------------------------------
    sf = shapefile.Reader(inputShapefile)
    bboxGeo = sf.shape(0).bbox

    # Convert coordinates from input projection to output projection
    bbox = [None] * 4
    bbox[0], bbox[1] = transformer.transform(bboxGeo[1],bboxGeo[0])
    bbox[2], bbox[3] = transformer.transform(bboxGeo[3], bboxGeo[2])

    #----------------------------------------------------------

    #----------------------------------------------------------
    # Read in ADCIRC mesh file
    #----------------------------------------------------------
    myMesh = pyadcircmodules.Mesh(inputMeshFile)
    ierr = myMesh.read()
    if ierr == 0:
        exit(ierr)
    
    for i in range(myMesh.numNodes()):
        nx = myMesh.node(i).x()
        ny = myMesh.node(i).y()
        
        nx2, ny2 = transformer.transform(ny,nx)
        myMesh.node(i).setX(nx2)
        myMesh.node(i).setY(ny2)
        
    myMesh.defineProjection(outEPSG,False)
        
    # myMesh.reproject(26919) This does not seem to work...
    
    # Convert coordinates
    xt,yt = transformer.transform(myMesh.y(),myMesh.x())

    #----------------------------------------------------------

    #----------------------------------------------------------
    # Find type and prepare ADCIRC data files for processing
    #----------------------------------------------------------
    useMeshZ = False
    useWaterLevelOutput = False
    useNodalAttribute = False
    useHarmonics = False
    isNetCDFFile = False
    isASCIIFile = False

    # Find the file extension
    fileExtension = os.path.splitext(inputDataFile)[1]
    if fileExtension == '.53':
        
        #print('Processing ' + inputDataFile)
        he=pyadcircmodules.HarmonicsOutput(inputDataFile); he.read();
        amp = np.zeros(shape=(he.numNodes(),he.numConstituents()))
        pha = np.zeros(shape=(he.numNodes(),he.numConstituents()))
        constitNames = []
        
            
        for i in range(he.numNodes()):
            for j in range(he.numConstituents()):
                amp[i][j] = he.amplitude(j).value(i)
                pha[i][j] = he.phase(j).value(i)
                
        # Create a metadata file
        metaDataFile = open('harmonics.README','w')
        print('raster.file\t:\t' + outputRasterFile + '.tif', file = metaDataFile)
        print('number.constituents\t:\t' + str(he.numConstituents()), file = metaDataFile)
        for i in range(he.numConstituents()):
            constitNames.append(he.name(i))
            # Print to metadatafile
            print(he.amplitude(i).name() + '.frequency\t:\t' + str(he.amplitude(i).frequency()), file = metaDataFile)
                
        useHarmonics = True

        #print('Success! \n')

    elif fileExtension == '.63':

        g = pyadcircmodules.ReadOutput(inputDataFile)
        g.open()

        if snap == 0  and g.numSnaps() > 1:
            print('Please select a snap from the *.63 file.')
            print('grd2dem.py --snap <snap> ...\n')
            sys.exit(2)
        else :
            snap = 1

        g.read()
        useWaterLevelOutput = True
        isASCIIFile = True

    elif fileExtension == '.13':

        n13 = pyadcircmodules.NodalAttributes(inputDataFile)
        n13.read()

        idx = n13.locateAttribute('mannings_n_at_sea_floor')

        useNodalAttribute = True

    elif fileExtension == '.nc' :

        fileName = os.path.splitext(inputDataFile)[0]
        fileExtension2 = os.path.splitext(fileName)[1]
        isNetCDFFile = True

        if fileExtension2 == '.63' :

            g = pyadcircmodules.ReadOutput(inputDataFile)
            g.open()

            if snap == 0  and g.numSnaps() > 1:
                print('Please select a snap from the *.63 file.')
                print('grd2dem.py --snap <snap> ...\n')
                sys.exit(2)
            else :
                snap = 1

            g.read(snap)

            useWaterLevelOutput = True

        else :
            print('Only .63.nc files are currently handeled by this code. Sorry.\n')
            sys.exit(2)

    elif inputDataFile == inputMeshFile:
        print('Creating raster from ' + inputMeshFile)
        useMeshZ = True
    #----------------------------------------------------------

        #----------------------------------------------------------
    # Main portion of the program.
    #
    # Much of the below needs to be inserted into a function/sub-routine.
    # Namely to compute the basis function
    #----------------------------------------------------------
    #print('Converting ADCIRC mesh to raster using ' + inputDataFile + '...')

    #if useMeshZ or useWaterLevelOutput or useNodalAttribute:
    if useMeshZ:
        
        values = np.empty((myMesh.numNodes(),1))
        values.fill(-99999)
        for n in range(myMesh.numNodes()):
            values[n] = myMesh.node(n).z()
        values[:,0] = values[:,0] * multFac
        values = np.squeeze(values)
        
        myMesh.toRaster(outputRasterFile+'.tif',
                        values,
                        bbox,
                        gridSize,
                        -99999,
                        "none",
                        "none",
                        False)
        
        tempRast = gdal.Open(outputRasterFile+'.tif', gdal.GA_Update)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(outEPSG)
        tempRast.SetProjection(srs.ExportToWkt())
        del tempRast 

    elif useWaterLevelOutput:

        values = np.empty((myMesh.numNodes(),1))
        values.fill(-99999)

        if isASCIIFile and everdried:
                        
            for n in range(myMesh.numNodes()):
                values[n] = g.data(snap-1).z(n)
                if abs(values[n]) > 0.95 and abs(values[n]) < 1.05:
                    values[n] = 1.0
                else:
                    values[n] = -99999
                    
            values[:,0] = values[:,0] * multFac
            values = np.squeeze(values)

        elif isNetCDFFile:

            for n in range(myMesh.numNodes()):
                values[n] = g.data(snap).z(n)
                
        elif isASCIIFile:
            
            for n in range(myMesh.numNodes()):
                values[n] = g.data(snap-1).z(n)
            
            values[:,0] = values[:,0] * multFac
            values = np.squeeze(values)
                   
        myMesh.toRaster(outputRasterFile+'.tif',
                        values,
                        bbox,
                        gridSize,
                        -99999,
                        "none",
                        "none",
                        False)
            
        tempRast = gdal.Open(outputRasterFile+'.tif', gdal.GA_Update)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(outEPSG)
        tempRast.SetProjection(srs.ExportToWkt())
        del tempRast
        
    elif useNodalAttribute:

        values = np.empty((myMesh.numNodes(),1))
        values.fill(-99999)

        for n in range(myMesh.numNodes()):
            values[n] = n13.attribute(idx,n).value(0)

        values[:,0] = values[:,0] * multFac
        values = np.squeeze(values)

        myMesh.toRaster(outputRasterFile+'.tif',
                        values,
                        bbox,
                        gridSize,
                        -99999,
                        "none",
                        "none",
                        False)
            
        tempRast = gdal.Open(outputRasterFile+'.tif', gdal.GA_Update)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(outEPSG)
        tempRast.SetProjection(srs.ExportToWkt())
        del tempRast 
        
    elif useHarmonics:

        for c in range(he.numConstituents()):
            
            values = np.empty((myMesh.numNodes(),1))
            valuesp = np.empty((myMesh.numNodes(),1))

            values = amp[:,c]
            valuesp = pha[:,c]

            values = values * multFac
            valuesp = valuesp * multFac
            
            myMesh.toRaster(
                outputRasterFile+'_'+constitNames[c]+'_amplitude.tif',
                values,
                bbox,
                gridSize,
                -99999,
                "none",
                "none",
                False)
            
            tempRast = gdal.Open(outputRasterFile+'_'+constitNames[c]+'_amplitude.tif', gdal.GA_Update)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(outEPSG)
            tempRast.SetProjection(srs.ExportToWkt())
            del tempRast 
            
            myMesh.toRaster(
                outputRasterFile+'_'+constitNames[c]+'_phase.tif',
                valuesp,
                bbox,
                gridSize,
                -99999,
                "none",
                "none",
                False)
            
            tempRast = gdal.Open(outputRasterFile+'_'+constitNames[c]+'_phase.tif', gdal.GA_Update)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(outEPSG)
            tempRast.SetProjection(srs.ExportToWkt())
            del tempRast
           
        # Combine into a single multi-banded raster            
        file_list = []
        for c in range(he.numConstituents()):
            file_list.append('harmonics_'+constitNames[c]+'_amplitude.tif')
            file_list.append('harmonics_'+constitNames[c]+'_phase.tif')
                       
        # Read metadata of first file
        with rasterio.open(file_list[0]) as src0:
            meta = src0.meta
            
            
        # Update meta to reflect the number of layers
        meta.update(count = len(file_list))
            
        # Read each layer and write it to stack
        with rasterio.open(outputRasterFile+'.tif', 'w', **meta) as dst:
            for id, layer in enumerate(file_list, start=1):
                with rasterio.open(layer) as src1:
                    dst.write_band(id, src1.read(1))
                    
        # Clean up
        for f in file_list:
            os.remove(f)
            
    #----------------------------------------------------------

#----------------------------------------------------------
# Clean up and close files.
#----------------------------------------------------------
    #dst_ds.FlushCache()
    if fileExtension == '.53': metaDataFile.close()

    #print('Success! \n')
#----------------------------------------------------------
