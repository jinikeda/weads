#!/usr/bin/python3
# File: preprocessing.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: Jul 15, 2024

# Development Notes:
""" This script is used to preprocess the input files for the WEAD project. The input files are the ADCIRC mesh file (fort.14), the attributes file (fort.13), the inundation time file (inundationtime.63), and the harmonics file (fort.53). The script reads the input files and filters the nodes within the domain shapefile (Cut_domain.shp). The script outputs the filtered nodes and attributes in the domain as a text file (domain_inputs.csv). The script is used in the WEAD project."""

##########################################################################
# --- Load internal modules ---

from .general_functions import *


def preprocessing_ADCIRC(inputMeshFile, inputAttrFile, inputInundationtimeFile,
                         inputHarmonicsFile, inputShapeFile, domainIOFile, inEPSG, outEPSG):

    #from src_point.general_functions import *
    # --- Initialize code ---
    start_time = time.time()
    print("\n")
    print("LAUNCH: Launching script!\n")

    ##########################################################################
    # --- Read a domain shapefile and coordinate conversions ---
    ##########################################################################

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0

    print("   Reading a domain ...\n")
    gdf = read_domain(inputShapeFile)
    gdf = gdf.set_crs(int(inEPSG))
    print(gdf.head(5))

    print("   Processing coordinate conversion...\n")
    PRJ = int(outEPSG)
    print("PRJ\t", PRJ)
    domain_prj = convert_gcs2coordinates(gdf, PRJ, None)

    # User's option for plot the domain
    # domain_prj.plot()
    # plt.show()

    ##########################################################################
    # --- Read mesh (inputMeshFile) ---
    ##########################################################################
    ADCIRC_nodes, numNodes, numElements = read_fort14(inputMeshFile)

    xy_list = ['x', 'y']
    drop_list = []
    crs_points = 'EPSG:' + str(inEPSG)

    # Create a GeoDataFrame for filtering
    gdf_ADCIRC = create_gdf(ADCIRC_nodes, xy_list, drop_list, crs_points, None)
    points_prj = convert_gcs2coordinates(gdf_ADCIRC, PRJ, None)
    print(points_prj.head(20))

    # Filter points within the domain
    points_within_polygon, true_indices = filter_points_within_domain(
        points_prj, domain_prj)
    ADCIRC_nodes_domain = ADCIRC_nodes[true_indices]
    # np.savetxt("true_indices.txt", true_indices, fmt='%d') # optional to turn on the output file
    print("number of nodes in the domain\t", len(ADCIRC_nodes_domain))
    inputElevation = "tbathy.txt"
    np.savetxt(inputElevation, ADCIRC_nodes_domain, fmt='%d\t%.8f\t%.8f\t%.8f')

    ##########################################################################
    # --- Read attributes (inputMeshFile) ---
    ##########################################################################
    mann, mann_indices, local_mann_indices, global_mann = read_fort13(
        inputAttrFile)
    mann_domain = mann[true_indices]
    # np.savetxt("mannings_n_at_sea_floor.txt", mann_domain, fmt='%f') # optional to turn on the output file
    print("manning_indices\t", mann_indices)

    ##########################################################################
    # --- Read inundationtime (inputInundationtimeFile) ---
    ##########################################################################
    inundationtime, nN, max_time = read_inundationtime63(
        inputInundationtimeFile)

    inundationtime_domain = inundationtime[true_indices]

    # Convert structured array to 2D array
    inundationtime_domain_2d = np.column_stack(
        [inundationtime_domain[name] for name in inundationtime_domain.dtype.names])
    np.savetxt("inundationtime.txt", inundationtime_domain_2d, fmt='%d\t%.4f')
    ##########################################################################
    ##########################################################################
    # --- Read harmonics (inputHarmonicsFile) ---
    ##########################################################################
    Harmonics_nodes, numHarm, nN, tidal_frequencies, tidal_constituents = read_fort53(
        inputHarmonicsFile)
    #np.savetxt("harmonics_AMP_pre.txt", Harmonics_nodes[:, :, 0], fmt='%.8f')
    #print("Harmonics_nodes dimensions: ", Harmonics_nodes.shape)
    print(tidal_constituents)
    Harmonics_nodes_domain = Harmonics_nodes[true_indices]
    np.savetxt("harmonics_Freq.txt", tidal_frequencies, fmt='%.12f')
    # np.savetxt("harmonics_AMP.txt",
    #            Harmonics_nodes_domain[:, :, 0], fmt='%.8f') # optional to turn on the output file
    # np.savetxt("harmonics_PHASE.txt",
    #            Harmonics_nodes_domain[:, :, 1], fmt='%.4f') # optional to turn on the output file

    ##########################################################################
    # --- Save the filtered nodes and attributes in the domain as a text file (domain_inputs.csv) ---
    ##########################################################################
    print("ADCIRC_nodes_domain dimensions: ", ADCIRC_nodes_domain.shape)
    print("mann_domain dimensions: ", mann_domain.shape)
    print("inundationtime_domain dimensions: ", inundationtime_domain_2d.shape)
    print("Harmonics_nodes_domain AMP dimensions: ",
          Harmonics_nodes_domain[:, :, 0].shape)
    print("Harmonics_nodes_domain PHASE dimensions: ",
          Harmonics_nodes_domain[:, :, 1].shape)
    merged_array = np.hstack((ADCIRC_nodes_domain,
                              mann_domain,
                              inundationtime_domain_2d,
                              Harmonics_nodes_domain[:,
                                                     :,
                                                     0],
                              Harmonics_nodes_domain[:,
                                                     :,
                                                     1]))
    dummy_node = 'node_2'

    header_list = ['node', 'x', 'y', 'z', 'mann', dummy_node, 'inundationtime']
    for i in ['_amp', '_phase']:
        for j in tidal_constituents:
            header_list.append(j + i)

    df = pd.DataFrame(merged_array, columns=header_list)
    df.drop(dummy_node, axis=1, inplace=True)
    #df.to_csv("domain_inputs.csv", index=None)

    ##########################################################################
    # --- Processing hydro_class based on inputMeshFile and inputInundationtimeFile ---
    ##########################################################################
    print(" Read an normalized inundation time (inpuT)")
    inputInundationTFile = "inundationtime.txt"
    inunT = np.loadtxt(
        inputInundationTFile,
        usecols=1)  # domain_inundationtime.txt
    print(inunT)

    print(" Read an tbathy (TB)")
    TB = np.loadtxt(inputElevation)

    # compare inundationtime and tbathy and Clean up values to only have 1 and
    # -99999
    accuracy = 1.0 * 10 ** -6

    """ inunTBN = -99999.0: nodata,
                0: land, 1: intertidal, 2: subtidal """  # , 3: pond/lake"""

    # Create the masks
    mask_land = (0 - accuracy < inunT) & (inunT < 0 + accuracy)  # land mask
    mask_intertidal = (accuracy < inunT) & (inunT < 1 - accuracy)  # intertidal mask
    mask_water = (1 - accuracy < inunT) & (inunT < 1 + accuracy)  # water mask

    # Forth column ("HydroClass")  on TB is at index 3
    mask_outdomain = (TB[:, 3] == ndv)  # initialize the value

    # Update the fourth column of TB based on the masks
    TB[mask_land, 3] = int(0)  # fully dried (land) region
    TB[mask_intertidal, 3] = int(1)  # intertidal region

    # Temporary set to water region and may separate to ocean (subtidal zone)
    # and pond/lake
    TB[mask_water, 3] = int(2)
    # set nodata value to -99999.0 in the domain of inunT
    TB[mask_outdomain, 3] = int(ndv)

    # Save the modified TB array to a file
    # np.savetxt("hydro_class.txt", TB, fmt='%d\t%.8f\t%.8f\t%d') # optional to turn on the output file

    # Reoder the columns of the DataFrame later Jin July 3, 2024
    df["HydroClass"] = TB[:, 3]

    df.to_csv(domainIOFile, index=None)

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Done preprosessing")
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")
    print("Job Finished ʕ •ᴥ•ʔ")

#preprocessing_ADCIRC_file("fort.14", "fort.13", "inundationtime.63", "fort.53", "Cut_domain.shp", 4269, 26914)
