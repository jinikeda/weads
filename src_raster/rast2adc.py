#!/usr/bin/python3
# File: rast2adc.py
# Name: Matthew V. Bilskie
# Date: March 8, 2022

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------
# ----------------------------------------------------------
import pyadcircmodules
import numpy as np
from osgeo import gdal
import src_raster.basics
from src_raster.raster import get_numrowcol, get_boundingbox,\
    get_nodatavalue, isinraster, coord2pixel,\
    get_rastvalue

# ----------------------------------------------------------
#
# ----------------------------------------------------------
# F U N C T I O N    R A S T 2 A D C
# ----------------------------------------------------------
#
# ----------------------------------------------------------


def rast2adc(inputMeshFile, outputMeshFile, inputRasterFile,
             outEPSG, band, multFac):

    # ----------------------------------------------------------
    # Check that input files exist
    # ----------------------------------------------------------
    src_raster.basics.fileexists(inputMeshFile)
    src_raster.basics.fileexists(inputRasterFile)
    # ----------------------------------------------------------
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Set up in/out projection
    # ----------------------------------------------------------
    #transformer = Transformer.from_crs(inEPSG,outEPSG)
    #rast = rast.reproject(outEPSG)
    # ----------------------------------------------------------

    # Compute number of rows and columns for output raster
    #numRows = math.ceil( (bbox[3] - bbox[1]) / gridSize )
    #numCols = math.ceil( (bbox[2] - bbox[0]) / gridSize )
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Read in ADCIRC mesh file
    # ----------------------------------------------------------
    mesh = pyadcircmodules.Mesh(inputMeshFile)
    print('Reading ADCIRC mesh file...')
    ierr = mesh.read()
    if ierr == 0:
        exit(ierr)
    print('Success! \n')
    # ----------------------------------------------------------

    rast = gdal.Open(inputRasterFile, gdal.GA_ReadOnly)
    rast = gdal.Warp("", rast, format="vrt", dstSRS="EPSG:" + str(outEPSG))
    numcols, numrows = get_numrowcol(rast)

    b = rast.GetRasterBand(band)
    vals = b.ReadAsArray()

    #bbox = get_boundingbox(rast)

    newZ = np.zeros(mesh.numNodes())

    noDataValue = get_nodatavalue(rast, band)

    for i in range(mesh.numNodes()):
        newZ[i] = mesh.node(i).z()

    for i in range(mesh.numNodes()):

        # if not isinraster(mesh.node(i).x(),mesh.node(i).y(),rast):
        # print(i,mesh.numNodes(),'Out')
        # continue

        col, row = coord2pixel(mesh.node(i).x(), mesh.node(i).y(), rast)

        if col <= 0:
            continue
        if row <= 0:
            continue
        if col >= numcols:
            continue
        if row >= numrows:
            continue

        # if get_rastvalue(col,row,rast,band) != noDataValue:
        if (vals[row][col] != noDataValue) and (vals[row][col] >
                                                1e-6):  # Threshold for biomass accumulation is 1e-6
            # print(i,mesh.numNodes(),'New')
            newZ[i] = multFac * vals[row][col] + mesh.node(i).z()
            # newZ[i] = multFac * vals[row][col] # coure grid polute adcirc mesh
            # print(mesh.node(i).z(),newZ[i])
        # else:
            # print(i,mesh.numNodes(),'Old')
            #newZ[i] = mesh.node(i).z()

    #print (max(vals.flatten()),min(vals.flatten()))
    mesh.setZ(newZ)
    mesh.write(outputMeshFile)
