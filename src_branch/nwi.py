#!/usr/bin/python3
# File: nwi.py
# coding: utf-8
# Gdal_reprojecion,resampling, etc.
# Developed by the Center for Computation & Technology and Center for Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# Developer: Jin Ikeda
# Last modified July 11, 2024

#######
# Caution: Segmentation fault might happen when read a large tiff file 

### Step 1 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 1.: Import modules")
print("----------------------------------------------------------------------\n")
######################################################################

#--- Load internal modules ---
from general_functions import *
import os, sys
from scipy import stats
import numpy as np
import pandas as pd
import time
import rasterio
from scipy import ndimage



# import spatial analysis
from osgeo import gdal,ogr,osr
from osgeo.gdalconst import *
gdal.UseExceptions()


# --- GLOBAL PARAMETERS ---
ndv = -99999.0  # No data value (ndv) using ADCIRC convention
ndv_byte = 128

# # Print the current working directory
# print("Current working directory: {0}".format(os.getcwd()))
#
# # The target Working directory
# # Workspace = "Z:/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO/inputs/Wetlands_NWI_forMeshRefinement/"
# Workspace = os.getcwd()
#
# # Change the current working directory
# os.chdir(Workspace)
# print("Current working directory":, Workspace)

### Function #####################################################################################################
def gdal_reading(file): # 0 GA_ReadOnly or 1
    Rasterdata = gdal.Open(file, GA_ReadOnly)
    if Rasterdata is None:
        print("Could not open " + Rasterdata)
        sys.exit(1)
    print("Reading raster file (georeference, etc)")

    # Coordinate system
    prj = Rasterdata.GetProjection()  # Read projection
    print("Projection:", prj)

    # Get raster size and band
    rows = Rasterdata.RasterYSize  # number of rows
    cols = Rasterdata.RasterXSize  # number of columns
    bandnum = Rasterdata.RasterCount  # band number
    print("rows=", rows, "cols=", cols)
    # print("band=", bandnum)

    # Get georeference info
    transform = Rasterdata.GetGeoTransform()
    xOrigin = transform[0]  # Upperleft x
    yOrigin = transform[3]  # Upperleft y
    pixelWidth = transform[1]  # cell size x
    pixelHeight = transform[5]  # cell size y (value is negative)
    print("xOrigin=", xOrigin, "m", "yOrigin=", yOrigin, "m")
    print("pixelWidth=", pixelWidth, "m", "pixelHeight=", -pixelHeight, "m")  # pixelHeight is always negative

    # Read the raster band
    band = Rasterdata.GetRasterBand(1)
    # Data type of the values
    print('data type is', gdal.GetDataTypeName(band.DataType))  # Each raster file has a different data type
    # Get band value info
    RV = band.ReadAsArray()          # raster values in the band
    del Rasterdata

    return prj,rows,cols,transform,RV

def create_raster_dilate(file, prj,rows,cols,transform,zarray, nodata):
    # Create the output raster dataset
    gtiff_driver = gdal.GetDriverByName('GTiff')
    #output_raster_path = os.path.join(Workspace, file)
    out_ds = gtiff_driver.Create(file, cols, rows, gdal.GDT_Byte)
    out_ds.SetProjection(prj)
    geotransform = list(transform)

    out_ds.SetGeoTransform(geotransform)
    # Write the resampled data to the output band
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(zarray)
    out_band.SetNoDataValue(ndv_byte)  # 0xff7fffee
    # out_band.ComputeStatistics(0)
    out_ds = None
    print('Each dilation has done')
    return file

def extract_point_values(raster_path, points_gdf, points_path):
    # Load the shapefile of points
    points = points_gdf
    # Load the DEM raster
    dem = rasterio.open(raster_path)

    # extract xy from point geometry
    raster_values = np.full(len(points), ndv_byte)

    array = dem.read(1)
    for i, point in enumerate(points.geometry):
        # print(point.xy[0][0],point.xy[1][0])
        x = point.xy[0][0]
        y = point.xy[1][0]
        row, col = dem.index(x, y)

        # Append the z value to the list of z values
        if 0 <= row < array.shape[0] and 0 <= col < array.shape[1]:
            raster_values[i] = array[row, col]

        # print("Point correspond to row, col: %d, %d"%(row,col))
        # print(array[row, col])
        # print("Raster value on point %.2f \n"%dem.read(1)[row,col])

    points['NWI'] = raster_values
    # points.to_file(points_path, driver='ESRI Shapefile')
    points.to_csv(points_path, index=False)  # Save the DataFrame to a CSV file with headers
    del points

    return raster_values

def read_tif(file,outEPSG,gridSize,skipresample_flag,deltaT=5):
    start = time.time()
    ### Step 2 ###########################################################
    print("\n----------------------------------------------------------------------")
    print("Step 2.: Reading files")
    print("----------------------------------------------------------------------\n")
    ######################################################################

    # Input raster data
    Raster_file = file
    print(Raster_file)

    prj,rows,cols,transform,band=gdal_reading(Raster_file) # Reading a raster file

    skipresample_flag = False

    if skipresample_flag == False:
        ### Step 3 ###########################################################
        print("\n----------------------------------------------------------------------")
        print("Step 3.: Reprojection for input raster")
        print("----------------------------------------------------------------------\n")
        ######################################################################
        # NEED to make a function later

        df = pd.read_csv("domain_inputs.csv")
        print(df.shape, df.columns, df.dtypes)

        df = df.loc [:, ['node','x','y',]]
        df ['NWI'] = ndv_byte

        #################################################################################
        inEPSG = 4269  # GCS_North_American_1983
        outEPSG = 26914
        #################################################################################

        xy_list = ['x', 'y']
        drop_list = []
        raster_prj = prj.split('],AUTHORITY["EPSG",')[-1].split(']]')[0] # Get EPSG code
        PRJ = 'EPSG:'+ raster_prj.replace('"', '')
        print ('Raster projection is', PRJ)
        crs_points = 'EPSG:'+str(inEPSG)

        # Create a GeoDataFrame for filtering
        gdf_ADCIRC = create_df2gdf(df, xy_list, drop_list, crs_points, None)
        points_prj = convert_gcs2coordinates(gdf_ADCIRC, PRJ, None)
        print(points_prj.head(20))

        points_path = 'domain_nwi.csv'

        extract_point_values(Raster_file, points_prj, points_path)

        # add dilation later Jin July 11, 2024

# internal command
# read_tif('NWI_TX_wetlands4m.tif', 26914, 40, False,deltaT=5.0)
read_tif('NWI_TX_wetlands.tif', 26914, 40, False,deltaT=5.0)