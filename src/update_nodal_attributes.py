#!/usr/bin/python3
# File: rast2adc.py
# Name: Matthew V. Bilskie
# Date: March 8, 2022

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import pyadcircmodules
import numpy as np
from osgeo import gdal
import src.basics
from src.raster import get_numrowcol,\
        get_nodatavalue, coord2pixel

#----------------------------------------------------------
#
#----------------------------------------------------------

#----------------------------------------------------------
# U P D A T E _ N O D A L _ A T T R I B U T E S
#----------------------------------------------------------
def update_nodal_attributes(inputMeshFile,inputRasterFile,\
             inputAttrFile,outputAttrFile,\
             outEPSG,slr,band):

    #----------------------------------------------------------
    # Check that input files exist
    #----------------------------------------------------------
    src.basics.fileexists(inputMeshFile)
    src.basics.fileexists(inputRasterFile)
    src.basics.fileexists(inputAttrFile)
    #----------------------------------------------------------
    #----------------------------------------------------------
    
    #----------------------------------------------------------
    # Read in ADCIRC mesh file
    #----------------------------------------------------------
    mesh = pyadcircmodules.Mesh(inputMeshFile)
    print('Reading ADCIRC mesh file...')
    ierr = mesh.read()
    if ierr == 0:
        exit(ierr)
    print('Success! \n')
    #----------------------------------------------------------
    
    #----------------------------------------------------------
    # Read in nodal attribute file
    #----------------------------------------------------------
    n13 = pyadcircmodules.NodalAttributes(inputAttrFile)
    n13.read()
    
    #----------------------------------------------------------
    # Update sea_surface_heigh_above_geoid attribute
    #----------------------------------------------------------
    idx_sshag = n13.locateAttribute('sea_surface_height_above_geoid')
    idx_manning = n13.locateAttribute('mannings_n_at_sea_floor')
    for n in range(mesh.numNodes()):
        sshag = n13.attribute(idx_sshag,n).value(0) + slr
        n13.attribute(idx_sshag,n).setValue(sshag)

    #----------------------------------------------------------
    # Update manning's n attribute based on biomass
    #----------------------------------------------------------
    values = np.empty((mesh.numNodes(),1))
    values.fill(-99999)

    rast = gdal.Open(inputRasterFile, gdal.GA_ReadOnly);
    rast = gdal.Warp("", rast, format="vrt", dstSRS="EPSG:"+str(outEPSG))
    numcols, numrows = get_numrowcol(rast)
    b = rast.GetRasterBand(band)
    vals = b.ReadAsArray()
    noDataValue = get_nodatavalue(rast,band)
    
    for i in range(mesh.numNodes()):
          
        col, row = coord2pixel(mesh.node(i).x(), mesh.node(i).y(), rast)
        
        if col <= 0:
            continue
        if row <= 0:
            continue
        if col >= numcols:
            continue
        if row >= numrows:
            continue
 
        if vals[row][col] != noDataValue:
            
            if vals[row][col] > 10 and vals[row][col] < 20:
                manning = 0.035
                continue
            if vals[row][col] > 20 and vals[row][col] < 28:
                manning = 0.05
                continue
            if vals[row][col] > 28 and vals[row][col] < 38:
                manning = 0.07
                continue
            
            n13.attribute(idx_manning,n).setValue(manning)
                
        
    n13.write(outputAttrFile)
