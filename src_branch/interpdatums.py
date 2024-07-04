#!/usr/bin/python3
# File: interpdatums.py
# Name: Jin Ikeda and Peter Bacopoulos
# Date: July 3, 2024
# Command line: python tidalDatumsIDW.py -i HyControl.img -s tidalDatums.img -o tidalDatumsIDW.img

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------
# ----------------------------------------------------------

import numpy as np
import time
import pandas as pd
import geopandas as gpd
from KDTree_idw import Invdisttree # Need KDTree_idw.py

def convert_gcs2coordinates(gdf, PRJ,Type_str):
    gdf_proj = gdf.to_crs(PRJ)
    if Type_str == "Point":
        gdf_proj["x"] = gdf_proj.geometry.apply(lambda point: point.x)
        gdf_proj["y"] = gdf_proj.geometry.apply(lambda point: point.y)
    elif Type_str == "Polygon":
        gdf_proj["x"] = gdf_proj.geometry.apply(lambda polygon: polygon.centroid.x)
        gdf_proj["y"] = gdf_proj.geometry.apply(lambda polygon: polygon.centroid.y)
    else:
        print("Type_str is not defined")
    return gdf_proj


# ----------------------------------------------------------
# F U N C T I O N	T I D A L D A T U M S I D W
# ----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums fully dried areas.
# result = function(inputRasterHyControl,inputRasterTidalDatums,\
#                   outputRaster,numNeighbors = 12)
# ----------------------------------------------------------
# def tidaldatumsidw(inputRasterHyControl, inputRasterTidalDatums, \
#                    outputRaster, numNeighbors=12):
'''
print ("")
print ("")
print ("--------------------------------------------------------")
print ("--- LAUNCHING HYDRO-MEM IN PYTHON - TIDAL DATUMS IDW ---")
print ("--------------------------------------------------------")
print ("")
'''

# --- GLOBAL PARAMETERS ---
start=time.time()
ndv = -99999.0  # No data value (ndv) using ADCIRC convention

inEPSG = 4269 # GCS_North_American_1983
outEPSG = 26914
#################################################################################
PRJ = "EPSG:" + str(outEPSG)

# --- READ INPUTS ---
# Read the CSV file
df = pd.read_csv("domain_tide.csv")
print ("  Read HyControl (HC) and tidal datums (TD) successfully")

print(df.shape,df.columns)

tidal_gdf = gpd.GeoDataFrame(df[['x','y','HydroClass','msl','mlw','mhw']], geometry=gpd.points_from_xy(df['x'], df['y'], crs='EPSG:'+str(inEPSG)))
tidal_prj = convert_gcs2coordinates(tidal_gdf, PRJ, None)
tidal_prj['x_prj'] = tidal_prj.geometry.apply(lambda point: point.x)
tidal_prj['y_prj'] = tidal_prj.geometry.apply(lambda point: point.y)
print(tidal_prj.head())

HC = tidal_prj['HydroClass'].values
# interpolate tidal datums (MSL, MLW, MHW) using IDW
tidal_prj['MSL_IDW'] = ndv
tidal_prj['MLW_IDW'] = ndv
tidal_prj['MHW_IDW'] = ndv

indices = np.where(HC < 2) # The indices of interpolated areas (not fully wetted areas)
inverse_indices = np.where(~(HC < 2)) # The indices of fully wetted areas (base points)
print('Nodes of interpolation:\t',len(indices[0]),', Nodes of references:\t', len(inverse_indices[0]))

# Use a projected system
scale_factor =  1000  # scale factor for KDTree unit: [m to 1km]
xy_base = tidal_prj[['x_prj', 'y_prj']].iloc[inverse_indices].to_numpy() / scale_factor
xy_interp = tidal_prj[['x_prj', 'y_prj']].iloc[indices].to_numpy() / scale_factor

print (len(xy_base))

# eps = .1  # approximate nearest, dist <= (1 + eps) * true nearest
# weights ~ 1 / distance**p
leafsize = 10  # leafsize for KDTree
power = 2  # the power of the inverse distance
knn = 12  # the number of nearest neighbors

# Call KDTree with the base points
for i in ['msl','mlw','mhw']:
    z = tidal_prj[i].values
    invdisttree = Invdisttree(xy_base, z[inverse_indices], leafsize=leafsize, stat=1)  # create i
    # interpolate tidal datums (MSL, MLW, MHW) using IDW
    interpol, weight_factors = invdisttree(xy_interp, nnear=knn, eps=0.0, p=power)  # interpolated using nnear numbers
    tidal_prj.loc[indices[0], i.upper()+'_IDW'] = interpol

tidal_prj.to_csv("tidal_prj.csv",index=False)

#--- Exit script ---
print("EXIT: Existing script!\n")
end=time.time(); print ("Time elapsed (seconds):",end-start);


