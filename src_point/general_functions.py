#!/usr/bin/python
# coding: utf-8
# Provide general functions for point-based WEADS
# Developed by the Center for Computation & Technology and Center for Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# Developer: Jin Ikeda, Peter Bacopoulos and Christopher E. Kees
# Last modified Jul 13, 2024

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import sys
import zipfile
import shutil
from scipy.spatial import KDTree
from datetime import datetime
import time
import glob
from pathlib import Path

# import spatial analysis
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import rasterio
gdal.UseExceptions()


def read_text_file(fileName):
    with open(fileName, "r") as f:
        lines = f.readlines()
    return lines


def read_fort14(inputMeshFile, output_Flag = False):
    print("Processing mesh...")
    with open(inputMeshFile, "r") as f:
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and elements  ######
    skip_index = 1  # skip the first line
    # eN:number of elemetns, #nN: number of nodes
    eN, nN = lines[skip_index].split()
    eN = int(eN)  # string to integer
    nN = int(nN)  # string to integer
    print("number of elements:eN={0},number of nodes:nN={1}".format(eN, nN))

    ##### Step.3 Output nodes ######
    nodes = []
    for i in range(nN):
        # skip eN and nN information
        nodeNum, x, y, z = lines[(skip_index + 1) + i].split()
        nodes.append([int(nodeNum), float(x), float(y),
                      float((-1) * float(z))])  # output longitude [degree],latitude [degree], and topobathy h [m, NAVD88] # Caution: ADCIRC z value (water depth direction is positive)

    ADCIRC_nodes = np.array(nodes)
    if output_Flag:
        np.savetxt("mesh_original.txt", ADCIRC_nodes, fmt='%d\t%.8f\t%.8f\t%.8f')
    return ADCIRC_nodes, nN, eN


def read_fort13(inputMeshFile):
    with open(inputMeshFile, "r") as f:
        lines = f.readlines()
    f.close()

    skip_index = 1  # skip the first line
    nN = int(lines[skip_index].split()[0])  # nN: number of nodes

    # skip_index2 = 2  # read for the number of attributes
    # numAttributes = int(lines[skip_index2].split()[0])  # number of attributes

    attribute = "mannings_n_at_sea_floor"
    mann_indices = []
    local_mann_indices = []
    global_mann = None  # initialize global_mann

    for i, line in enumerate(lines):
        if attribute in line:
            mann_indices.append(i)
            if not global_mann:  # if global_mann is not yet set
                global_mann = float(lines[i + 3].split()[0])
                print("global_manning\t", global_mann)
                mann = np.full((nN, 1), global_mann, dtype=float)  # initialize using a global value
            else:
                local_mann_Num = int(lines[i + 1].split()[0])
                print("local_manning_Num\t", local_mann_Num)
                for j in range(local_mann_Num):
                    NN, value = map(float, lines[i + 2 + j].split())
                    local_mann_indices.append((int(NN) - 1))
                    mann[int(NN) - 1][0] = value  # overwrite the value

    return mann, mann_indices, local_mann_indices, global_mann


# Jin's comments. We also need to change the function to read netcdf file.
def read_inundationtime63(inputInundationTFile):
    with open(inputInundationTFile, "r") as f:
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and maximum simulation time ######
    skip_index = 1  # skip the first line
    nN = int(lines[skip_index].split()[1])  # nN: number of nodes
    skip_index2 = 2  # skip the first line
    # max_time: maximum simulation time [s]
    max_time = float(lines[skip_index2].split()[1])
    print("node number\t", nN, "max_time\t", max_time)

    ##### Step.3 Output inundationtime ######
    # inundationtime = np.zeros(nN, dtype=[('nodeNum', int), ('time', float)])
    # node number, time of inundation (0: dry, 1: wet)
    inundationtime = np.zeros(nN, dtype=[('nodeNum', int), ('time', float)])
    for i in range(nN):
        # skip before inundation number
        nodeNum, time = lines[(skip_index2 + 1) + i].split()
        inundationtime[i][0] = int(nodeNum)
        if float(time) >= 0:
            # float(float(time) / max_time)  # time of inundation (0: dry, 1:
            # wet)
            inundationtime[i][1] = float(time) / max_time
        else:
            pass

    return inundationtime, nN, max_time


def read_fort53(inputHarmonicsFile):
    print("   Processing harmonics...\n")
    # Jin's note: We may need to change the function to read netcdf file.
    with open(inputHarmonicsFile, "r") as f:
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and maximum simulation time ######
    skip_index = 0  # skip the first line
    numHarm = int(lines[skip_index].split()[0])  # numHarm: number of Harmonics

    tidal_frequncies = []
    tidal_constitunents = []

    skip_index2 = 1  # skip the first line
    for i in range(numHarm):
        # tidal frequencies (rad/s)
        tidal_frequncies.append(float(lines[skip_index2 + i].split()[0]))
        tidal_constitunents.append(
            lines[skip_index2 + i].split()[3])  # tidal constituents

    skip_index3 = 1 + numHarm  # skip the line until the number of nodes
    nN = int(lines[skip_index3].split()[0])  # nN: number of nodes
    print("number of Harmonics\t", numHarm, ", node number\t", nN)

    # Tidal harmonics of each node
    # 3D array [AMP, PHASE] cannot save as a text file
    Harmonics_nodes = np.zeros((nN, numHarm, 2), dtype=float)

    for i in range(nN):
        # print('check',lines[skip_index3 + 1 + i * (numHarm + 1)].split())
        for ii in range(numHarm):
            # skip before EMAGT(k,j), PHASEDE(k,j) here number of node is
            # skipped.
            AMP, PHASE = lines[skip_index3 + 2 +
                               i * (numHarm + 1) + ii].split()
            Harmonics_nodes[i][ii][0] = float(AMP)
            Harmonics_nodes[i][ii][1] = float(PHASE)
            # print(ii, Harmonics_nodes[i][ii][0], Harmonics_nodes[i][ii][1])

    return Harmonics_nodes, nN, numHarm, tidal_frequncies, tidal_constitunents


def read_domain(inputShapeFile):
    gdf = gpd.read_file(inputShapeFile)
    return gdf


def create_gdf(ADCIRC_nodes, xy_list, drop_list, crs, output_file=None):
    # Convert the numpy array to a pandas DataFrame
    df = pd.DataFrame(ADCIRC_nodes, columns=["node_id", "x", "y", "z"])
    # Or Read the text file into a pandas DataFrame
    # df = pd.read_csv(ADCIRC_Mesh_File, sep="\t", header=None, names=["node_id", "x", "y", "z"])
    gdf = gpd.GeoDataFrame(df.drop(drop_list, axis=1), geometry=gpd.points_from_xy(
        df[xy_list[0]], df[xy_list[1]], crs=crs))

    if output_file is not None:
        gdf.to_file(output_file, driver='ESRI Shapefile')

    return gdf


def create_df2gdf(df, xy_list, drop_list, crs, output_file=None):
    gdf = gpd.GeoDataFrame(df.drop(drop_list, axis=1), geometry=gpd.points_from_xy(
        df[xy_list[0]], df[xy_list[1]], crs=crs))

    if output_file is not None:
        gdf.to_file(output_file, driver='ESRI Shapefile')

    return gdf


def convert_gcs2coordinates(gdf, PRJ, Type_str): # Type_str means convert the type of vector data.
    gdf_proj = gdf.to_crs(PRJ)
    if Type_str == "Point":
        gdf_proj["x_prj"] = gdf_proj.geometry.apply(lambda point: point.x)
        gdf_proj["y_prj"] = gdf_proj.geometry.apply(lambda point: point.y)
    elif Type_str == "Polygon":
        gdf_proj["x"] = gdf_proj.geometry.apply(
            lambda polygon: polygon.centroid.x)
        gdf_proj["y"] = gdf_proj.geometry.apply(
            lambda polygon: polygon.centroid.y)
    else:
        print("Type_str is not defined")
    return gdf_proj


def convert_coordinates2gcs(gdf_proj, PRJ):
    gdf = gdf_proj.to_crs(PRJ)
    gdf["x"] = gdf.geometry.apply(lambda point: point.x)
    gdf["y"] = gdf.geometry.apply(lambda point: point.y)
    return gdf


def filter_points_within_domain(points_gdf, polygon_gdf):
    # Perform a spatial join between the points and the polygons
    points_within_domain = points_gdf[points_gdf.geometry.within(
        polygon_gdf.unary_union)]

    # Get the indices of the rows in points_within_domain
    true_indices = points_within_domain.index

    return points_within_domain, true_indices


def gdal_reading(file):  # 0 GA_ReadOnly or 1
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
    print("pixelWidth=", pixelWidth, "m", "pixelHeight=", -
          pixelHeight, "m")  # pixelHeight is always negative

    # Read the raster band
    band = Rasterdata.GetRasterBand(1)
    # Data type of the values
    # Each raster file has a different data type
    print('data type is', gdal.GetDataTypeName(band.DataType))
    # Get band value info
    RV = band.ReadAsArray()          # raster values in the band

    return prj, rows, cols, transform, RV, Rasterdata


def create_raster(file, rasterdata, zarray, dtype,
                  no_data_value, stats_flag=False):
    # Create the output raster dataset
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(
        file,
        rasterdata.RasterXSize,
        rasterdata.RasterYSize,
        rasterdata.RasterCount,
        dtype)  # dtype is e.g. gdal.GDT_Int32 and gdal.GDT_Float32
    out_ds.SetProjection(rasterdata.GetProjection())
    out_ds.SetGeoTransform(rasterdata.GetGeoTransform())
    dst_band = out_ds.GetRasterBand(1)
    dst_band.WriteArray(zarray)
    dst_band.SetNoDataValue(no_data_value)  # Exclude nodata value
    stats = dst_band.ComputeStatistics(0)
    min_val, max_val, mean_val, std_dev_val = stats
    if stats_flag:
        print(
            f'Made a raster file. Statistics:\n\tMinimum: {min_val}, Maximum: {max_val}, Mean: {mean_val}, Standard Deviation: {std_dev_val}')
    else:
        print('Made a raster file')
    out_ds = None

    return

def extract_point_values(raster_path, points_gdf, points_path, ndv,ndv_byte):
    # Load the shapefile of points
    points = points_gdf
    # Load the DEM raster
    dem = rasterio.open(raster_path)

    # extract xy from point geometry
    raster_id = np.zeros(len(points))
    raster_values = np.full(len(points), ndv_byte)

    array = dem.read(1)
    for i, point in enumerate(points.geometry):
        # print(point.xy[0][0],point.xy[1][0])
        x = point.xy[0][0]
        y = point.xy[1][0]
        row, col = dem.index(x, y)

        # Append the z value to the list of z values
        if 0 <= row < array.shape[0] and 0 <= col < array.shape[1]:
            raster_id[i] = row * array.shape[1] + col
            raster_values[i] = array[row, col]
        else:
            raster_id[i] = ndv

        # print("Point correspond to row, col: %d, %d"%(row,col))
        # print(array[row, col])
        # print("Raster value on point %.2f \n"%dem.read(1)[row,col])

    points['Raster_id'] = raster_id
    points['NWI'] = raster_values
    # points.to_file(points_path, driver='ESRI Shapefile')
    # Save the DataFrame to a CSV file with headers
    points.to_csv(points_path, index=False)

    return raster_values, points


def expand_nodes(nodes_positions, node_states,
                 target_value, infection_distance):
    # nodes_positions: the positions of the nodes (x, y coordinates) as a numpy array (N,2)
    # node_states: the states of the nodes (0: not infected, target_value:
    # infected) as a numpy array (N,1)

    # Create a KDTree for efficient neighbor search
    tree = KDTree(nodes_positions)

    # Find all nodes within the infection distance of the infected node(s) and
    # update their states
    # Get the indices of infected nodes here node_state a tuple
    potential_expansion_nodes = np.where(node_states == target_value)[0]
    for node in potential_expansion_nodes:
        neighbors = tree.query_ball_point(
            nodes_positions[node],
            infection_distance,
            p=2.0)  # p=2 (circle) for Euclidean distance
        for neighbor in neighbors:
            if neighbor != node:  # Exclude the node itself
                node_states[neighbor] = target_value

    return node_states
