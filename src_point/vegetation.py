#!/usr/bin/python3
# File: vegetation.py
# Gdal reading,extract raster values,potential expansion, etc.
# Developed by the Center for Computation & Technology and Center for Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# Developer: Jin Ikeda
# Last modified Aug 24, 2024

########################################################################################################################
### Caution: Segmentation fault might happen when read a large tiff file
########################################################################################################################

# ----------------------------------------------------------
# F U N C T I O N    P R O C E S S V E G E T A T I O N F I L E
# ----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW),
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = process_vegetation_file(inputvegetationFile, skip_raster_extracting_Flag,
#                             spread_flag, outputvegetationFile, inEPSG, outEPSG, deltaT=5)
# deltaT is user defined time on the argument file.

# --- Load internal modules ---
from .general_functions import *

### Step 1 ###########################################################
print("\n----------------------------------------------------------------------")
print("Step 1.: Import modules")
print("----------------------------------------------------------------------\n")
######################################################################

# --- GLOBAL PARAMETERS ---
ndv = -99999.0  # No data value (ndv) using ADCIRC convention
ndv_byte = 128

### Functions #############################################################


def extract_point_values(raster_path, points_gdf, points_path):
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


def calculate_potential_expansion(
        node_positions, NWI_values, target_value, radius):
    Vegetation_mask = np.where(
        NWI_values == target_value,
        target_value,
        0).flatten()  # salt marsh (regularly flooded) 8
    Potential_expansion = expand_nodes(
        node_positions,
        Vegetation_mask,
        target_value,
        radius)  # refer general function
    indices = np.where(Potential_expansion == target_value)[0]
    return Potential_expansion, indices


def apply_priority_order(df, indices_values, output_file):
    df.rename(columns={'NWI': 'NWI_org'},
              inplace=True)  # Rename the column name
    # Copy the original NWI values to potential NWI values
    df['NWI'] = df['NWI_org']
    for indices, value in indices_values:
        for i in indices:
            # Assign the potential NWI values to the original NWI values
            df.loc[i, 'NWI'] = value
    # Save the DataFrame to a CSV file with headers
    df.to_csv(output_file, index=False)
    return df


def process_vegetation_file(inputvegetationFile, skip_raster_extracting_Flag,
                            spread_flag, outputvegetationFile, inEPSG, outEPSG, deltaT=5):

    start_time = time.time()

    #outputvegetationFile = 'domain_nwi.csv'

    ### potential expansion speed
    speed = 40  # [m/year]   # User defined values
    radius = speed * deltaT  # [m per deltaT years calculation]
    print(
        'The allowable range of expansion speed per simulation is',
        radius,
        '[m-deltaT-yrs]\n')

    ### Step 2 ###########################################################
    print("\n----------------------------------------------------------------------")
    print("Step 2.: Reading files")
    print("----------------------------------------------------------------------\n")
    ######################################################################

    # Input raster data (.tiff)
    print(inputvegetationFile)

    prj, rows, cols, transform, RV, _ = gdal_reading(
        inputvegetationFile)  # Reading a raster file
    print('row and cols: ', np.shape(RV))

    xy_list = ['x', 'y']
    drop_list = []
    raster_prj = prj.split(
        '],AUTHORITY["EPSG",')[-1].split(']]')[0]  # Get EPSG code
    PRJ = 'EPSG:' + raster_prj.replace('"', '')
    print('Raster projection is', PRJ)
    crs_points = 'EPSG:' + str(inEPSG)

    if skip_raster_extracting_Flag == False:

        ### Step 3 ###########################################################
        print("\n----------------------------------------------------------------------")
        print("Step 3.: Reprojection for input points")
        print("----------------------------------------------------------------------\n")
        ######################################################################

        df = pd.read_csv("domain_inputs.csv")
        print(df.shape, df.columns, df.dtypes)

        df = df.loc[:, ['node', 'x', 'y', 'z']]
        df['NWI'] = ndv_byte # add no_values

        # Create a GeoDataFrame for filtering
        gdf_ADCIRC = create_df2gdf(df, xy_list, drop_list, crs_points, None)
        points_prj = convert_gcs2coordinates(gdf_ADCIRC, PRJ, "Point") # Convert point-based datasets
        print(points_prj.head(20))

        process_file = 'domain_nwi_original.csv'

        # extract values from raster file
        NWI_values, NWI_df = extract_point_values(
            inputvegetationFile, points_prj, process_file)

    else:

        # reading simulated mem file
        df = pd.read_csv('previous_ecology.csv') # for sequential running
        df.rename(columns={'new_NWI': 'NWI'}, inplace=True) # rename and use as original values
        NWI_values = df['NWI'].values

        ### Step 3 ###########################################################
        print("\n----------------------------------------------------------------------")
        print("Step 3.: Reprojection for input points")
        print("----------------------------------------------------------------------\n")
        ######################################################################

        NWI_df = df.loc[:, ['node', 'x', 'y', 'z', 'NWI']]
        print(NWI_df.columns)

        # Create a GeoDataFrame for filtering
        gdf_ADCIRC = create_df2gdf(
            NWI_df, xy_list, drop_list, crs_points, None)
        points_prj = convert_gcs2coordinates(gdf_ADCIRC, PRJ, "Point")
        print(points_prj.head(20))

    if spread_flag == True: # calculate vegetation expansions

        Point_x = points_prj[['x_prj']]
        Point_y = points_prj[['y_prj']]
        node_positions = np.column_stack((Point_x, Point_y))

        assert Point_x.shape[0] == len(
            NWI_values), 'The number of points and NWI values are not matched'

        print("\n----------------------------------------------------------------------")
        print("Step 4.: Spread vegetation nodes and make potential expansion")
        print("----------------------------------------------------------------------\n")

        # Apply binary dilation using scipy.ndimage
        Potential_SRF, SRF_indices = calculate_potential_expansion(
            node_positions, NWI_values, 8, radius)  # salt marsh (regularly flooded) 8
        Potential_MG, MG_indices = calculate_potential_expansion(
            node_positions, NWI_values, 9, radius)  # mangrove 9
        Potential_SIRF, SIRF_indices = calculate_potential_expansion(
            node_positions, NWI_values, 20, radius)  # irregularly flooded marsh 20

        # Save the potential expansion nodes

        # priority order #
        # mask1 = (MG == 1)  # priority order
        # mask2 = (MG != 1) & (SRF == 1)
        # mask3 = (MG != 1) & (SRF != 1) & (SIRF == 1)

        NWI_df = apply_priority_order(NWI_df, [(SIRF_indices, 20), (SRF_indices, 8), (
            MG_indices, 9)], outputvegetationFile)  # order from low to high priority and save as a single csv file

    if spread_flag == False:
        shutil.copy(process_file, outputvegetationFile)

    print('The NWI data has been saved')

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")

##########################################################################
# internal command
# process_vegetation_file('NWI_TX_wetlands4m.tif', False,True,deltaT=5.0) #process_vegetation_file(inputvegetationFile,skip_raster_extracting_Flag,spread_flag,deltaT=5):
# process_vegetation_file('NWI_TX_wetlands.tif', False, True, deltaT=5.0) #