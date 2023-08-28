#!/usr/bin/env python
# coding: utf-8
# Gdal_reprojecion,resampling
# Developed by the Center for Computation & Technology and Center for Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# Developer: Jin Ikeda, Peter Bacopoulos, and Christopher E. Kees
# Last modified August 26, 2023

### Step 1 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 1.: Import modules")
print("----------------------------------------------------------------------\n")
######################################################################

# import os module
import os, sys
from scipy import stats
import numpy as np
import pandas as pd
import warnings
import geopandas as gpd
import glob
import time
import rasterio
from scipy import ndimage

# import spatial analysis
from osgeo import gdal, ogr,osr
from osgeo.gdalconst import *

start=time.time()

print("Done")

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

### Step 2 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 2.: Reading files")
print("----------------------------------------------------------------------\n")
######################################################################

# --- GLOBAL PARAMETERS ---
gridSize=100
ndv=-99999.0                      # No data value (ndv) using ADCIRC convention
ndv_byte=128
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# The target Working directory
#Workspace = "Z:/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO/inputs/Wetlands_NWI_forMeshRefinement/"
Workspace =os.getcwd()

# Change the current working directory
os.chdir(Workspace)

# Input raster data
Raster_file = os.path.join(Workspace, "Region3_NWI_LC_Reclassify_wetlandsOnly.tif")
prj,rows,cols,transform,band=gdal_reading(Raster_file) # Reading a raster file

### Step 3 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 3.: Reprojection")
print("----------------------------------------------------------------------\n")
######################################################################
# NEED to make a function later
reprj_out="reproject_NN_gdal.tif"

# Desired spatial reference
output_srs = osr.SpatialReference()
output_srs.ImportFromEPSG(26919)  # Set the desired coordinate system using EPSG code
print(output_srs.ExportToWkt())  # Print the WKT representation of the spatial reference

# Parameters for output raster
output_pixel_width = transform[1]  # Keep the same pixel width
output_pixel_height = -transform[5]  # Keep the same pixel height (make it positive)
output_cols = cols
output_rows = rows

# Warp options
warp_options = gdal.WarpOptions(
    format="GTiff",
    xRes=output_pixel_width,
    yRes=output_pixel_height,
    srcSRS=prj,
    dstSRS=output_srs.ExportToWkt(),
    resampleAlg=gdal.GRA_NearestNeighbour  # Nearest neighbor resampling
)

# Perform the warp operation
output_raster_path = os.path.join(Workspace, reprj_out)
Rasterdata = gdal.Open(Raster_file, GA_ReadOnly)
output_ds = gdal.Warp(output_raster_path, Rasterdata, options=warp_options)
output_ds = None
Rasterdata = None

output_raster_path = os.path.join(Workspace, reprj_out)
print(output_raster_path)

# Open the output raster dataset to add statistics
out_ds = gdal.Open(output_raster_path, gdal.GA_Update)

# Check if the dataset was opened successfully
if out_ds is None:
    print("Failed to open output raster dataset.")
else:
    # Get the first band of the raster dataset
    out_band = out_ds.GetRasterBand(1) # 0 means read-only. 1 means writeable.
    out_band.ComputeStatistics(0)
    out_ds = None

### Step 4 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 4.: Resample raster size")
print("----------------------------------------------------------------------\n")
######################################################################
# NEED to make a function later

resample_factor = 50 # Define the resampling factor

output_raster_path = os.path.join(Workspace,reprj_out)
in_ds = gdal.Open(output_raster_path, 0)
in_band = in_ds.GetRasterBand(1)
original_data=in_band.ReadAsArray()
#original_data = np.random.randint(0, 3, size=(1000, 1000))  # Example random data
print('Original row and cols: ',np.shape(original_data))


# Calculate the dimensions of the resampled array
resampled_rows = original_data.shape[0] // resample_factor
resampled_cols = original_data.shape[1] // resample_factor

print('Resampled row and cols: ',resampled_rows,resampled_cols)

# Calculate the new dimensions to ensure divisibility by the resample factor
end_rows = resampled_rows * resample_factor
end_cols = resampled_cols * resample_factor

# Crop the original data to match the new dimensions
original_data = original_data[:end_rows, :end_cols]

# Initialize an empty array to store the resampled blocks
resampled_blocks = []

# Loop through the rows and columns to create blocks of size (resample_factor, resample_factor)
for i in range(0, end_rows, resample_factor):
    for j in range(0, end_cols, resample_factor):
        block = original_data[i:i+resample_factor, j:j+resample_factor]
        resampled_blocks.append(block)

# Initialize an empty array to store the mode values for each block
mode_values = []

# Loop through each block and calculate the mode value
for block in resampled_blocks:
    mode_value = stats.mode(block.flatten(), keepdims=True)[0]
    mode_values.append(mode_value)

mode_values = np.array(mode_values).reshape(resampled_rows,resampled_cols)
#mode_values[mode_values > 127] = ndv

#print(mode_values)

# Create the output raster dataset
gtiff_driver = gdal.GetDriverByName('GTiff')
output_raster_path = os.path.join(Workspace, 'Region3_NWI_LC_UTM19N_resample100.tif')
out_ds = gtiff_driver.Create(output_raster_path, resampled_cols, resampled_rows,gdal.GDT_Byte)
out_ds.SetProjection(in_ds.GetProjection())
geotransform = list(in_ds.GetGeoTransform())
geotransform[1] *= resample_factor  # Adjust pixel width
geotransform[5] *= resample_factor  # Adjust pixel height

print('resample size: ',geotransform[1])
out_ds.SetGeoTransform(geotransform)

# Write the resampled data to the output band
out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(mode_values)
out_band.SetNoDataValue(ndv_byte)
out_band.ComputeStatistics(0)
out_ds = None

### Step 5 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 5.: Create a new raster for study area domain")
print("----------------------------------------------------------------------\n")
######################################################################

# domain raster data
Domain_raster_file = os.path.join(Workspace,"../hyconn.tif")
prj2,rows2,cols2,transform2,hc=gdal_reading(Domain_raster_file) # Reading a raster file on domain
prj3,rows3,cols3,transform3,band3=gdal_reading(output_raster_path) # Reading a original marsh classification file

xOrigin = transform2[0]  # Upperleft x
yOrigin = transform2[3]  # Upperleft y
pixelWidth = transform2[1]  # cell size x
pixelHeight = transform2[5]  # cell size y (value is negative)

# Create points
p=0
nodes = np.zeros((rows2*cols2,2),'d')
for j in range(rows2):
    for i in range(cols2):
        nodes[p,0]=xOrigin + (i+0.5) * pixelWidth # Get center values
        nodes[p,1]=yOrigin + (j+0.5) * pixelHeight # Get center values
        p+=1
headers = ["x", "y"]
df = pd.DataFrame(nodes, columns=headers)

# Open raster file
with rasterio.open(output_raster_path) as marsh_Raster:
    z = []
    for idx, point in df.iterrows():
        x, y = point['x'], point['y']
        row, col = marsh_Raster.index(x, y)
        # Check if the row and column indices are within valid range
        #print(row,col)
        if 0 <= row < rows3 and 0 <= col < cols3:
            z.append(marsh_Raster.read(1)[row, col])
        else:
            z.append(ndv)  # Assign a placeholder value for out-of-range indices

# Reshape the array to match the shape of the original grid
z = np.array(z)  # Convert the list to a NumPy array
df['z']=z
df.to_csv('domain_xyz.csv', index=False) # Save the DataFrame to a CSV file with headers
z = z.reshape(rows2, cols2)
z[z < -9999.] = ndv_byte  # Replace values less than -9999 with 128
#print(z)

#
# # gdf = gpd.GeoDataFrame(nodes, columns=headers, geometry=gpd.points_from_xy(nodes[:, 0], nodes[:, 1]), crs="EPSG:26919")
# # gdf.to_file('domain_point.shp', driver='ESRI Shapefile')
# #
# # #open point shapefile
# # pointData = gpd.read_file('domain_point.shp')
#
# Create the output raster dataset
gtiff_driver = gdal.GetDriverByName('GTiff')
output_raster_path = os.path.join(Workspace, 'Resampled_raster_domain.tif')
out_ds = gtiff_driver.Create(output_raster_path, cols2, rows2, gdal.GDT_Byte)
out_ds.SetProjection(prj2)
geotransform = list(transform2)

out_ds.SetGeoTransform(geotransform)
# Write the resampled data to the output band
out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(z)
out_band.SetNoDataValue(ndv_byte) #0xff7fffee
out_band.ComputeStatistics(0)
out_ds = None

### Step 6 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 6.: Determine the range of  extend area")
print("----------------------------------------------------------------------\n")
######################################################################
# dilate the
speed = 200 #[m/year]
range = speed*5 #[m per 5 years calculation]
print('The allowable range of expansion speed per simulation is',range, '[m-5yrs]\n')
#
# Input raster data
prj,rows,cols,transform,RV=gdal_reading(output_raster_path) # Reading a raster file

# Create the kernel
num=range//gridSize
size_karnel=(int(2*num+1)) # should be odd
print ('dilation pixel is', size_karnel)
kernel = np.ones((size_karnel, size_karnel), np.uint8)

# Apply binary dilation using scipy.ndimage
SRF = ndimage.binary_dilation(RV == 8, structure=kernel, iterations=1) # salt marsh (regularly flooded) 8
MG  = ndimage.binary_dilation(RV == 9, structure=kernel, iterations=1) # mangrove 9
SIRF= ndimage.binary_dilation(RV == 20, structure=kernel, iterations=1) # irregularly flooded marsh 20

def create_raster_dilate(file, zarray, nodata):
    # Create the output raster dataset
    gtiff_driver = gdal.GetDriverByName('GTiff')
    output_raster_path = os.path.join(Workspace, file)
    out_ds = gtiff_driver.Create(output_raster_path, cols, rows, gdal.GDT_Byte)
    out_ds.SetProjection(prj)
    geotransform = list(transform)

    out_ds.SetGeoTransform(geotransform)
    # Write the resampled data to the output band
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(zarray)
    out_band.SetNoDataValue(ndv_byte) #0xff7fffee
    #out_band.ComputeStatistics(0)
    out_ds = None
    print('Each dilation has done')
    return output_raster_path

SRF_raster_path = create_raster_dilate('SRF_resample100_dilation.tif',SRF,ndv_byte)
MG_raster_path = create_raster_dilate('MG_resample100_dilation.tif',MG,ndv_byte)
SIRF_raster_path = create_raster_dilate('SIRF_resample100_dilation.tif',SIRF,ndv_byte)

tolerance = 1e-6  # Set a tolerance threshold for dimension comparison

if abs(rows2 - rows) > tolerance or abs(cols2 - cols) > tolerance:
    print('something wrong\n')
    print('hyconn row and cols:', rows2, ':', cols2)
    print('marsh row and cols:', rows, ':', cols)
else:
    Domain = np.full((rows2,cols2), ndv_byte) # Create an array of default values (ndv)

    # We merge five classification into a single values
    water_mask = ((-0.5 < hc) & (hc < 0.5)) | (1.5 < hc) # water regions (water and pond in hyconn.tiff)
    land_mask = (hc >= 0.5) & (hc <= 1.5) # land regions in hy conn.tiff

    # Create a mask for different conditions
    mask1 = (MG == 1) # priority order
    mask2 = (MG != 1) & (SRF==1)
    mask3 = (MG != 1) & (SRF!=1) & (SIRF==1)

    # Assign values based on conditions
    # background raster
    Domain[land_mask] = 55

    # Marsh productivity
    Domain[mask1] = 9 # mangrove 9
    Domain[mask2] = 8 # salt marsh (regularly flooded) 8
    Domain[mask3] = 20 # irregularly flooded marsh 20

    Domain[water_mask] = 40 # Need to consider later Aug 26th

    SRF_raster_path = create_raster_dilate('Domain_classification_distribution_resample100.tif', Domain,ndv_byte) # We are still using Byte
    print('Potential vegetation classification map has done')

#--- Exit script ---
print("\n"); print("EXIT: Existing script!\n");
end=time.time(); print ("Time elapsed (seconds):",end-start);
print("\n"); sys.exit();
# ######################################################################################################################
