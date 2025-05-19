#!/usr/bin/env python
# File: tidaldatumsidw.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: May 9, 2025

# ----------------------------------------------------------
# F U N C T I O N	T I D A L D A T U M S I D W
# ----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums not fully wet areas.
# result = function(outputHarmonicsFile,interpolateHarmonicsFile,inEPSG,outEPSG, numNeighbors=12):
# ----------------------------------------------------------


##########################################################################
# --- Load internal modules ---
from .general_functions import *
from .KDTree_idw import Invdisttree  # Need KDTree_idw.py


def tidaldatumsidw(outputHarmonicsFile, interpolateHarmonicsFile,
                   inEPSG, outEPSG, numNeighbors=12):

    start_time = time.time()

    '''
    print ("")
    print ("")
    print ("--------------------------------------------------------")
    print ("-------------- LAUNCHING TIDAL DATUMS IDW --------------")
    print ("--------------------------------------------------------")
    print ("")
    '''

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0  # No data value (ndv) using ADCIRC convention

    PRJ = int(outEPSG)
    print("PRJ\t", PRJ)

    # --- READ INPUTS ---
    # Read the CSV file
    df = pd.read_csv(outputHarmonicsFile)
    print("  Read HyControl (HC) and tidal datums (TD) successfully")

    print(df.shape, df.columns)

    tidal_gdf = gpd.GeoDataFrame(df[['node',
                                     'x',
                                     'y',
                                     'z',
                                     'HydroClass',
                                     'inundationtime',
                                     'msl',
                                     'mlw',
                                     'mhw']],
                                 geometry=gpd.points_from_xy(df['x'],
                                                             df['y'],
                                                             crs=int(inEPSG)))  # crs='EPSG:'+str(inEPSG)
    # Convert the projection
    tidal_prj = convert_gcs2coordinates(tidal_gdf, PRJ, None)
    tidal_prj['x_prj'] = tidal_prj.geometry.apply(lambda point: point.x) # add projected x values
    tidal_prj['y_prj'] = tidal_prj.geometry.apply(lambda point: point.y) # add projected y values
    print(tidal_prj.head())

    HC = tidal_prj['HydroClass'].values
    # interpolate tidal datums (MSL, MLW, MHW) using IDW
    tidal_prj['MSL_IDW'] = ndv
    tidal_prj['MLW_IDW'] = ndv
    tidal_prj['MHW_IDW'] = ndv

    # The indices of interpolated areas (not fully wetted areas)
    indices = np.where(HC < 2)
    # The indices of fully wetted areas (base points)
    inverse_indices = np.where(~(HC < 2))
    print(
        'Nodes of interpolation:\t', len(
            indices[0]), ', Nodes of references:\t', len(
            inverse_indices[0]))

    # Use a projected system
    scale_factor = 1000  # scale factor for KDTree unit: [m to 1km]
    xy_base = tidal_prj[['x_prj', 'y_prj']
                        ].iloc[inverse_indices].to_numpy() / scale_factor
    xy_interp = tidal_prj[['x_prj', 'y_prj']
                          ].iloc[indices].to_numpy() / scale_factor

    print(len(xy_base))

    # Use KDTree for IDW interpolation
    leafsize = 10  # leafsize for KDTree
    power = 2  # the power of the inverse distance
    knn = numNeighbors  # the number of nearest neighbors

    # Call KDTree with the base points
    for i in ['msl', 'mlw', 'mhw']:
        z = tidal_prj[i].values
        invdisttree = Invdisttree(
            xy_base,
            z[inverse_indices],
            leafsize=leafsize,
            stat=0)  # create i

        # interpolate tidal datums (MSL, MLW, MHW) using IDW
        interpol, weight_factors = invdisttree(
            xy_interp, nnear=knn, eps=0.0, p=power)  # interpolated using k nnear numbers
        tidal_prj.loc[indices[0], i.upper() + '_IDW'] = interpol # add interpolated values in the dataframe

    tidal_prj.to_csv(interpolateHarmonicsFile, index=False)

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # --- PRINT TO SCREEN ---
    print("\n----------------------------------------------------")
    print("----------- COMPLETED - TIDAL DATUMS IDW ------------")
    print("----------------------------------------------------\n")
    print(f"Time to Compute: {elapsed_time} seconds")
    print("Job Finished ʕ •ᴥ•ʔ")
