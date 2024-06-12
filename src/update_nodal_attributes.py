#!/usr/bin/python3
# File: rast2adc.py
# Name: Matthew V. Bilskie, Jin Ikeda
# Date: August 23, 2023

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------
# ----------------------------------------------------------
import pyadcircmodules
import numpy as np
from osgeo import gdal
import src.basics
from src.raster import get_numrowcol, \
    get_nodatavalue, coord2pixel


# ----------------------------------------------------------
#
# ----------------------------------------------------------

# ----------------------------------------------------------
# U P D A T E _ N O D A L _ A T T R I B U T E S
# ----------------------------------------------------------
def update_nodal_attributes(inputMeshFile, inputRasterFile, \
                            inputAttrFile, outputAttrFile, \
                            outEPSG, slr, band):
    # ----------------------------------------------------------
    # Check that input files exist
    # ----------------------------------------------------------
    src.basics.fileexists(inputMeshFile)
    src.basics.fileexists(inputRasterFile)
    src.basics.fileexists(inputAttrFile)
    # ----------------------------------------------------------
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

    # ----------------------------------------------------------
    # Read in nodal attribute file
    # ----------------------------------------------------------
    n13 = pyadcircmodules.NodalAttributes(inputAttrFile)
    n13.read()

    # ----------------------------------------------------------
    # Update sea_surface_heigh_above_geoid attribute
    # ----------------------------------------------------------
    idx_sshag = n13.locateAttribute('sea_surface_height_above_geoid')
    idx_manning = n13.locateAttribute('mannings_n_at_sea_floor')
    print ('slr \t =', slr)
    for n in range(mesh.numNodes()):
#        sshag = n13.attribute(idx_sshag, n).value(0) + slr # this is probably not good idea for sequential runs
        sshag = slr # directly add the diff from the base year: Modified by Jin Ikeda 06/08/2024
        n13.attribute(idx_sshag, n).setValue(sshag)

    print(n13.attribute(idx_sshag, n).value(0),sshag)

    # ----------------------------------------------------------
    # Update manning's n attribute based on biomass
    # ----------------------------------------------------------
    values = np.empty((mesh.numNodes(), 1))
    values.fill(-99999)

    rast = gdal.Open(inputRasterFile, gdal.GA_ReadOnly);
    rast = gdal.Warp("", rast, format="vrt", dstSRS="EPSG:" + str(outEPSG))
    numcols, numrows = get_numrowcol(rast)
    b = rast.GetRasterBand(band)
    vals = b.ReadAsArray()
    noDataValue = get_nodatavalue(rast, band)

    for i in range(mesh.numNodes()):
       col, row = coord2pixel(mesh.node(i).x(), mesh.node(i).y(), rast)

       if 0 < col < numcols and 0 < row < numrows and vals[row][col] != noDataValue:
           if 10 <= vals[row][col] < 20: # 16: low productivity
               manning = 0.035
           elif 20 <= vals[row][col] < 28: # 23 medium productivity
               manning = 0.05
           elif 28 <= vals[row][col] < 38: # 32 high productivitiy
               manning = 0.07
           else:
               continue  # Skip if value doesn't match any condition
            
           n13.attribute(idx_manning, i).setValue(manning)

    n13.write(outputAttrFile)



