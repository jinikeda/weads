#!/usr/bin/python3
# File: preprocessing.py
# Name: Jin Ikeda & Peter Bacopoulos
# Date: June 28, 2024

# Development Notes:
""" This script is used to preprocess the input files for the WEAD project. The input files are the ADCIRC mesh file (fort.14), the attributes file (fort.13), the inundation time file (inundationtime.63), and the harmonics file (fort.53). The script reads the input files and filters the nodes within the domain shapefile (Cut_domain.shp). The script outputs the filtered nodes and attributes in the domain as a text file (domain_inputs.csv). The script is used in the WEAD project."""

########################################################################################################################
# module and function should be merged once the codes are completed
#--- Load modules ---
import numpy as np
import time
import geopandas as gpd
import pandas as pd

# import spatial analysis using Gdal
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import rasterio
gdal.AllRegister()  # Register all of drivers
gdal.UseExceptions()  # Enable exceptions

#--- Initialize code ---
start_time = time.time()
print("\n"); print("LAUNCH: Launching script!\n")

########################################################################################################################
def read_fort14(inputMeshFile):
    print("Processing mesh...")
    with open(inputMeshFile, "r") as f:
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and elements  ######
    skip_index = 1  # skip the first line
    eN, nN = lines[skip_index].split()  # eN:number of elemetns, #nN: number of nodes
    eN = int(eN)  # string to integer
    nN = int(nN)  # string to integer
    print("number of elements:eN={0},number of nodes:nN={1}".format(eN, nN))

    ##### Step.3 Output nodes ######
    nodes = []
    for i in range(nN):
        nodeNum, x, y, z = lines[(skip_index + 1) + i].split()  # skip eN and nN information
        nodes.append([int(nodeNum), float(x), float(y),
                      float((-1)*float(z))])  # output longitude [degree],latitude [degree], and topobathy h [m, NAVD88] # Caution: ADCIRC z value (water depth direction is positive)

    ADCIRC_nodes = np.array(nodes)
    np.savetxt("mesh_original.txt", ADCIRC_nodes, fmt='%d\t%.8f\t%.8f\t%.8f')
    return ADCIRC_nodes, nN, eN

def read_fort13(inputMeshFile):
    with open(inputMeshFile, "r") as f:
        lines = f.readlines()
    f.close()

    skip_index = 1  # skip the first line
    nN = int(lines[skip_index].split()[0])  # nN: number of nodes

    skip_index2 = 2  # read for the number of attributes
    numAttributes = int(lines[skip_index2].split()[0])  # nN: number of nodes

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
                mann = np.full((nN, 1), global_mann, dtype=float)
            else:
                local_mann_Num = int(lines[i + 1].split()[0])
                print("local_manning_Num\t", local_mann_Num)
                for j in range(local_mann_Num):
                    NN, value = map(float, lines[i + 2 + j].split())
                    local_mann_indices.append((int(NN) - 1))
                    mann[int(NN) - 1][0] = value

    return mann, mann_indices, local_mann_indices, global_mann

def read_inundationtime63(inputInundationTFile): # Jin's comments. We may need to change the function to read netcdf file.
    with open(inputInundationTFile, "r") as f:
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and maximum simulation time ######
    skip_index = 1  # skip the first line
    nN = int(lines[skip_index].split()[1])  # nN: number of nodes
    skip_index2 = 2  # skip the first line
    max_time = float(lines[skip_index2].split()[1])  # max_time: maximum simulation time [s]
    print("node number\t", nN, "max_time\t", max_time)

    ##### Step.3 Output inundationtime ######
    inundationtime = np.zeros(nN, dtype=[('nodeNum', int), ('time', float)])  # node number, time of inundation (0: dry, 1: wet)
    for i in range(nN):
        nodeNum, time = lines[(skip_index2 + 1) + i].split()  # skip before inundation number
        inundationtime[i][0] = int(nodeNum)
        if float(time) >= 0:
            inundationtime[i][1] = float(time) / max_time  #float(float(time) / max_time)  # time of inundation (0: dry, 1: wet)
        else:
            pass

    return inundationtime, nN, max_time

def read_fort53(inputHarmonicsFile):
    print("   Processing harmonics...\n")
    with open(inputHarmonicsFile, "r") as f: # Jin's note: We may need to change the function to read netcdf file.
        lines = f.readlines()
    f.close()

    ##### Step.2 Get the number of nodes and maximum simulation time ######
    skip_index = 0  # skip the first line
    numHarm = int(lines[skip_index].split()[0])  # numHarm: number of Harmonics

    tidalconstitunents = []
    skip_index2 = 1  # skip the first line
    for i in range(numHarm):
        tidalconstitunents.append(lines[skip_index2 + i].split()[3])

    skip_index3 = 1+numHarm  # skip the first line
    nN = int(lines[skip_index3].split()[0])  # nN: number of nodes
    print("number of Harmonics\t",numHarm,"node number\t", nN)

    # Tidal harmonics of each node
    Harmonics_nodes = np.zeros((nN, numHarm, 2), dtype=float) # 3D array [AMP, PHASE] cannot save as a text file

    for i in range(nN):
        for ii in range(numHarm):
            AMP, PHASE = lines[skip_index3 + 2 + i*(numHarm+1) + ii].split() # skip before inundation number
            Harmonics_nodes[i][ii][0] = float(AMP)
            Harmonics_nodes[i][ii][1] = float(PHASE)

    return Harmonics_nodes, nN, numHarm,tidalconstitunents

def read_domain(inputShapeFile):
    gdf = gpd.read_file(inputShapeFile)
    return gdf

def create_gdf(ADCIRC_nodes,xy_list, drop_list, crs, output_file=None):
    # Convert the numpy array to a pandas DataFrame
    df = pd.DataFrame(ADCIRC_nodes, columns=["node_id", "x", "y", "z"])
    # Or Read the text file into a pandas DataFrame
    # df = pd.read_csv(ADCIRC_Mesh_File, sep="\t", header=None, names=["node_id", "x", "y", "z"])
    gdf = gpd.GeoDataFrame(df.drop(drop_list, axis=1), geometry=gpd.points_from_xy(df[xy_list[0]], df[xy_list[1]], crs=crs))

    if output_file != None:
        gdf.to_file(output_file, driver='ESRI Shapefile')

    return gdf

def create_df2gdf(df, xy_list, drop_list, crs, output_file=None):
    gdf = gpd.GeoDataFrame(df.drop(drop_list, axis=1), geometry=gpd.points_from_xy(df[xy_list[0]], df[xy_list[1]], crs=crs))

    if output_file != None:
        gdf.to_file(output_file, driver='ESRI Shapefile')

    return gdf

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
def convert_coordinates2gcs(gdf_proj, PRJ):
    gdf = gdf_proj.to_crs(PRJ)
    gdf["x"] = gdf.geometry.apply(lambda point: point.x)
    gdf["y"] = gdf.geometry.apply(lambda point: point.y)
    return gdf

def filter_points_within_domain(points_gdf, polygon_gdf):
    # Perform a spatial join between the points and the polygons
    points_within_domain = points_gdf[points_gdf.geometry.within(polygon_gdf.unary_union)]

    # Get the indices of the rows in points_within_domain
    true_indices = points_within_domain.index

    return points_within_domain, true_indices

########################################################################################################################
#--- Read a domain shapefile and coordinate conversions ---
########################################################################################################################

# --- GLOBAL PARAMETERS ---
ndv = -99999.0


print("   Reading a domain ...\n")
gdf = read_domain("Cut_domain.shp")

print("   Processing coordinate conversion...\n")
# Convert inEPSG(GCS) to outEPSG(target projection)
### This part should be provided from the parser file #### Jin June 27, 2024
inEPSG = 4269 # GCS_North_American_1983
outEPSG = 26914
#################################################################################
PRJ = "EPSG:" + str(outEPSG)
domain_prj = convert_gcs2coordinates(gdf, PRJ, None)

# domain_prj.plot()
# plt.show()

########################################################################################################################
#--- Read mesh (inputMeshFile) ---
########################################################################################################################
ADCIRC_nodes, numNodes, numElements = read_fort14("TX2008_T35H_Refined_052224.grd")

xy_list = ['x', 'y']
drop_list = []
crs_points = 'EPSG:'+str(inEPSG)

# Create a GeoDataFrame for filtering
gdf_ADCIRC = create_gdf(ADCIRC_nodes,xy_list, drop_list, crs_points, None)
points_prj = convert_gcs2coordinates(gdf_ADCIRC, PRJ, None)
print(points_prj.head(20))

# Filter points within the domain
points_within_polygon, true_indices = filter_points_within_domain(points_prj, domain_prj)
ADCIRC_nodes_domain = ADCIRC_nodes[true_indices]
np.savetxt("true_indices.txt", true_indices, fmt='%d')
print("number of nodes in the domain\t", len(ADCIRC_nodes_domain))
np.savetxt("tbathy.txt", ADCIRC_nodes_domain, fmt='%d\t%.8f\t%.8f\t%.8f')

########################################################################################################################
#--- Read attributes (inputMeshFile) ---
########################################################################################################################
mann, mann_indices, local_mann_indices, global_mann = read_fort13("TX2008_T35H_Refined_052224.13")
mann_domain = mann[true_indices]
np.savetxt("mannings_n_at_sea_floor.txt", mann_domain, fmt='%f')
print("manning_indices\t", mann_indices)

########################################################################################################################
#--- Read inundationtime (inputInundationtimeFile) ---
########################################################################################################################
inundationtime, nN, max_time = read_inundationtime63("inundationtime.63")

inundationtime_domain = inundationtime[true_indices]
np.savetxt("inundationtime.txt", inundationtime_domain, fmt='%d\t%.4f')
########################################################################################################################

### Pending to be implemented
# #--- Read everdried ---
# print("\n"); print("   Processing everdried...\n");
# myFile=open("everdried.63","r")
# myLine=myFile.readline()
# myLine=myFile.readline()
# myLine=myFile.readline()
# ed=np.ones((numNodes,1),dtype=float)
# for j in range(numNodes):
#     myLine=myFile.readline(); myRow=myLine.split(); ed[j][0]=float(myRow[1]);
# myFile.close()


########################################################################################################################
# Pending to be implemented

# #--- Read maxinundepth ---
# print("   Processing maxinundepth...\n")
# myFile=open("maxinundepth.63","r")
# myLine=myFile.readline()
# myLine=myFile.readline()
# myLine=myFile.readline()
# mi=np.ones((numNodes,1),dtype=float)
# for j in range(numNodes):
#     myLine=myFile.readline(); myRow=myLine.split(); mi[j][0]=float(myRow[1]);
# myFile.close()

########################################################################################################################
#--- Read harmonics (inputHarmonicsFile) ---
########################################################################################################################
# Harmonics_nodes, numHarm, nN,tidal_constituents = read_fort53("fort.53")
# print(tidal_constituents)
# Harmonics_nodes_domain = Harmonics_nodes[true_indices]
# np.savetxt("harmonics_AMP.txt", Harmonics_nodes_domain[:, :, 0], fmt='%.8f')
# np.savetxt("harmonics_PHASE.txt", Harmonics_nodes_domain[:, :, 1], fmt='%.4f')

########################################################################################################################
# merged_array = np.hstack((ADCIRC_nodes_domain, mann_domain, inundationtime_domain, Harmonics_nodes_domain[:, :, 0], Harmonics_nodes_domain[:, :, 1]))
# dummy_node ='node_2'
#
# header_list =['node','x','y','z', 'mann', dummy_node,'inundationtime']
# for i in ['_amp','_phase']:
#     for j in tidal_constituents:
#         header_list.append(j+i)
#
# df = pd.DataFrame(merged_array, columns=header_list)
# df.drop(dummy_node, axis=1, inplace=True)
# df.to_csv("domain_inputs.csv", index=None)

#######################################################################################################################
##


# --- READ INPUTS ---
# print ("Reading ADCIRC node")

print(" Read an normalized inundation time (inpuT)")
inputInundationTFile = "inundationtime.txt"
inunT = np.loadtxt(inputInundationTFile, usecols=1) # domain_inundationtime.txt
print (inunT)

print(" Read an tbathy (TB)")
inputElevation = "tbathy.txt"
TB = np.loadtxt(inputElevation)

# compare inundationtime and tbathy and Clean up values to only have 1 and -99999
accuracy = 1.0 * 10 ** -6

""" inunTBN = -99999.0: nodata,
            0: land, 1: intertidal, 2: subtidal """  # , 3: pond/lake"""

# Create the masks
mask_land = (0 - accuracy < inunT) & (inunT < 0 + accuracy)
mask_intertidal = (accuracy < inunT) & (inunT < 1 - accuracy)
mask_water = (1 - accuracy < inunT) & (inunT < 1 + accuracy)
mask_outdomain = (TB[:, 3] == ndv)  # Assuming the fourth column is at index 3

# Update the fourth column of TB based on the masks
TB[mask_land, 3] = int(0)  # fully dried (land) region
TB[mask_intertidal, 3] = int(1)  # intertidal region
TB[mask_water, 3] = int(2)  # temporary set to water region and will separate to ocean (subtidal zone) and pond/lake
TB[mask_outdomain, 3] = int(ndv)  # set nodata value to -99999.0 in the domain of inunT

# Save the modified TB array to a file
np.savetxt("hydro_class.txt", TB, fmt='%d\t%.8f\t%.8f\t%d')

########################################################################################################################
# Calculate the elapsed time
end_time = time.time()
elapsed_time = end_time - start_time

# Print the elapsed time
print("Done preprosessing")
print("Time to Compute: \t\t\t", elapsed_time, " seconds")
print("Job Finished ʕ •ᴥ•ʔ")


