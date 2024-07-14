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
from scipy import stats
import rasterio
from scipy import ndimage

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
    points.to_csv(points_path, index=False)  # Save the DataFrame to a CSV file with headers
    del points

    return raster_values

def read_tif(Raster_file,Resample_gridSize,skipresample_flag,diration_flag,deltaT=5):

    start_time = time.time()

    ### Step 2 ###########################################################
    print("\n----------------------------------------------------------------------")
    print("Step 2.: Reading files")
    print("----------------------------------------------------------------------\n")
    ######################################################################

    # Input raster data
    print(Raster_file)

    prj,rows,cols,transform,RV, _ = gdal_reading(Raster_file) # Reading a raster file

    # Create a resampled raster file
    if skipresample_flag == False:
        ### Step 3 ###########################################################
        print("\n----------------------------------------------------------------------")
        print("Step 3-1.: Resample raster size")
        print("----------------------------------------------------------------------\n")
        ######################################################################

        # Resample the raster data to speed up the dilation process
        original_data = RV.copy()

        print('Original row and cols: ', np.shape(original_data))

        # Calculate the dimensions of the resampled array
        resample_factor = int(Resample_gridSize / transform[1])  # Define the resampling factor only available interger
        print("Convert size factor is: ", resample_factor)
        resampled_rows = original_data.shape[0] // resample_factor
        resampled_cols = original_data.shape[1] // resample_factor

        print('Resampled row and cols: ', resampled_rows, resampled_cols)

        # Calculate the new dimensions to ensure divisibility by the resample factor
        end_rows = int(resampled_rows * resample_factor)
        end_cols = int(resampled_cols * resample_factor)

        # Crop the original data to match the new dimensions
        original_data = original_data[:end_rows, :end_cols]

        # Initialize an empty array to store the resampled blocks
        resampled_blocks = []

        # Loop through the rows and columns to create blocks of size (resample_factor, resample_factor)
        for i in range(0, end_rows, resample_factor):
            for j in range(0, end_cols, resample_factor):
                block = original_data[i:i + resample_factor, j:j + resample_factor]
                resampled_blocks.append(block)

        # Initialize an empty array to store the mode values for each block
        mode_values = []

        # Loop through each block and calculate the mode value
        for block in resampled_blocks:
            mode_value = stats.mode(block.flatten(), keepdims=True)[0]
            mode_values.append(mode_value)

        mode_values = np.array(mode_values).reshape(resampled_rows, resampled_cols)
        # mode_values[mode_values > 127] = ndv

        # print(mode_values)

        # Create the output raster dataset
        gtiff_driver = gdal.GetDriverByName('GTiff')
        Resample_raster_path = Raster_file.split('.')[0]+'_resample' + str(int(Resample_gridSize)) + 'm.tif'
        out_ds = gtiff_driver.Create(Resample_raster_path, resampled_cols, resampled_rows, gdal.GDT_Byte)
        out_ds.SetProjection(prj)
        geotransform = list(transform)
        geotransform[1] *= int(resample_factor)  # Adjust pixel width
        geotransform[5] *= int(resample_factor)  # Adjust pixel height

        print('Resample pixel size: ', geotransform[1])
        out_ds.SetGeoTransform(geotransform)

        # Write the resampled data to the output band
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(mode_values)
        out_band.SetNoDataValue(ndv_byte)
        out_band.ComputeStatistics(0)
        out_ds = None

    # Create a new raster file with dilation
    if diration_flag == True:

        print("\n----------------------------------------------------------------------")
        print("Step 4.: Dilate raster files")
        print("----------------------------------------------------------------------\n")

        # dilate the
        speed = 40  # [m/year]
        area = speed * deltaT  # [m per deltaT years calculation]
        print('The allowable range of expansion speed per simulation is', area, '[m-5yrs]\n')

        if skipresample_flag == True:
            Raster_path = Raster_file
        else:
            Raster_path = Resample_raster_path

        prj,rows,cols,transform,RV, _ = gdal_reading(Raster_path)  # Reading a raster file

        # Create the kernel
        gridSize = transform[1]  # [m]
        num = area // gridSize
        size_karnel = (int(2 * num + 1))  # should be odd
        print('Dilation pixel is', size_karnel)
        kernel = np.ones((size_karnel, size_karnel), np.uint8)

        # Apply binary dilation using scipy.ndimage
        SRF = ndimage.binary_dilation(RV == 8, structure=kernel, iterations=1)  # salt marsh (regularly flooded) 8
        MG = ndimage.binary_dilation(RV == 9, structure=kernel, iterations=1)  # mangrove 9
        SIRF = ndimage.binary_dilation(RV == 20, structure=kernel, iterations=1)  # irregularly flooded marsh 20

        SRF_raster_path = 'SRF_resample' + str(int(transform[1])) + '_dilation.tif'
        MG_raster_path = 'MG_resample' + str(int(transform[1])) + '_dilation.tif'
        SIRF_raster_path = 'SIRF_resample' + str(int(transform[1])) + '_dilation.tif'

        create_raster_dilate(SRF_raster_path, prj, rows, cols, transform, SRF, ndv_byte)
        create_raster_dilate(MG_raster_path, prj, rows, cols, transform, MG, ndv_byte)
        create_raster_dilate(SIRF_raster_path, prj, rows, cols, transform, SIRF, ndv_byte)

        tolerance = 1e-6  # Set a tolerance threshold for dimension comparison

        # Read the CSV file
        df = pd.read_csv("tidal_prj.csv")
        print("  Read HyControl (HC) and tidal datums (TD) successfully")
        print(df.shape, df.columns, df.dtypes)

        hc = df['HydroClass']  # 'HydroClass' 0: land, 1: intertidal, 2: subtidal(water)

        # Create a new raster file with dilation
        Domain = np.full((rows, cols), ndv_byte)  # Create an array of default values (ndv)

        # Create a mask for different conditions
        mask1 = (MG == 1)  # priority order
        mask2 = (MG != 1) & (SRF == 1)
        mask3 = (MG != 1) & (SRF != 1) & (SIRF == 1)

        # Assign values based on conditions

        # Marsh productivity
        Domain[mask1] = 9  # mangrove 9
        Domain[mask2] = 8  # salt marsh (regularly flooded) 8
        Domain[mask3] = 20  # irregularly flooded marsh 20

        # Point-based approach is not able to handle the background process

        # # We merge five classification into a single values
        # water_mask = ((1.5 < hc) & (hc < 2.5) | (0.5 < hc) & (hc < 1.5) | (
        #             2.5 < hc))  # submergence region (hc = 1.0, 2.0, 3.0 including lake and pond)
        # land_mask = (-0.5 < hc) & (hc <= 0.5)  # land region (hc = 0.0) # land regions in hydro_class.tif
        # Domain[land_mask] = 55
        # Domain[water_mask] = 40  # Need to consider later Aug 26th

        Output = 'Domain_classification_distribution_dilation' + str(int(transform[1])) + 'm.tif'
        SRF_raster_path = create_raster_dilate(Output, prj, rows, cols, transform, Domain,
                                               ndv_byte)  # We are still using Byte
        print('Potential vegetation classification map has done')

    ######################################################################

    ### Step 4 ###########################################################
    print("\n----------------------------------------------------------------------")
    print("Step 4.: Reprojection for input points")
    print("----------------------------------------------------------------------\n")
    ######################################################################
    # NEED to make a function later

    df = pd.read_csv("domain_inputs.csv")
    print(df.shape, df.columns, df.dtypes)

    df = df.loc [:, ['node','x','y','z']]
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

    ########################################################################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Done ecological response")
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")



# internal command
# read_tif('NWI_TX_wetlands4m.tif', 26914, 40,False, False,deltaT=5.0) #read_tif(Raster_file,Resample_gridSize,skipresample_flag,diration_flag,deltaT=5):
#read_tif('NWI_TX_wetlands_resample10m.tif',10, True, True, deltaT=5.0) #

# # Create a function to get the mode with priority (internal function for the mode_calculate2raster)
# def prioritized_mode(series, priority_order):
#     mode_values = series.mode()
#     #print(mode_values)
#     for value in priority_order:
#         if value in mode_values.values:
#             return value
#     return mode_values.iloc[0]
#
# # For renew raster value need to assign the cell id and value
# def mode_calculate2raster(file_path,target_list, drop_value, priority_order, output_file): # target_list = [target group, target value]= ['Raster_id', 'NWI']
#
#     df = pd.read_csv(file_path) # Read the CSV file
#     df = df[df[target_list[0]] != drop_value] # Remove rows with ndv values
#
#     # Get the unique values and their counts
#     value_counts = df[target_list[0]].value_counts()
#
#     # Get target_list[0] values with count greater than one
#     print(value_counts[value_counts > 1]) # Print values with count greater than one
#     raster_id_over_2 = value_counts[value_counts > 1].index
#
#     # Select rows in the DataFrame where target_list[0] is in 'raster_id' greater than one
#     df_id_deplicate = df[df[target_list[0]].isin(raster_id_over_2)]
#     df_process = df.loc[~df[target_list[0]].isin(raster_id_over_2), target_list] # Remove rows with target_list[0] in 'raster_id greater than one'
#
#     # Group the DataFrame by 'Raster_id'
#     grouped = df_id_deplicate.groupby(target_list[0])
#
#     # Initialize a list to store the results
#     results = []
#
#     for raster_id, group in grouped:
#         nwi_mode = prioritized_mode(group['NWI'], priority_order)
#         results.append([raster_id, nwi_mode])
#
#     # Create a DataFrame from the results
#     mode_df = pd.DataFrame(results, columns=target_list)
#     # mode_df.to_csv("check_mode.csv", index=False)
#
#     # Concatenate the mode_df and df
#     df_merge = pd.concat([df_process, mode_df], ignore_index=True)
#     df_merge.to_csv(output_file, index=False)
#
# mode_calculate2raster("domain_nwi.csv",['Raster_id','NWI'], ndv, [9, 8, 20, 128], "domain_nwi_mode.csv")