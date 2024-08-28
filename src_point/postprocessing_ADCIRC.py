#!/usr/bin/python3
# File: postprocessing.py
# Developer: Jin Ikeda
# Last modified: Aug 27, 2024

##########################################################################
# --- Load internal modules ---
from .general_functions import *

def update_ADCIRC_mesh(outputMeshFile, node_id, x, y, new_z):
    """
    Updates the z-values in an ADCIRC mesh file based on input arrays.

    Parameters:
    - outputMeshFile (str): Path to the ADCIRC mesh file to be updated.
    - node_id (array-like): Array of node indices.
    - x (array-like): X-coordinates of the nodes.
    - y (array-like): Y-coordinates of the nodes.
    - new_z (array-like): New z-values to update in the mesh file.

    """
    # Read the entire file into a list of lines
    lines = read_text_file(outputMeshFile)

    # Update the z value for the ADCIRC mesh
    for i in range(len(node_id)):
        idx = int(node_id[i])
        # Update the line corresponding to the node index
        # Node indices in the file are assumed to start from line 2 (index 1 in Python)
        lines[idx + 1] = "{0:>10}     {1:.6f}      {2:.6f} {3:.8E}\n".format(
            idx, x[i], y[i], -new_z[i])

    # Determine the output file path
    output_path = outputMeshFile

    # Write the updated lines back to the file
    with open(output_path, 'w') as file:
        file.writelines(lines)
    file.close()

    print(f"Updated ADCIRC file '{output_path}' with new node data from mem results")


def update_ADCIRC_attributes(outputAttrFile, slr, manning_node, manning):
    """
    Updates sea_surface_height_above_geoid and Manning's n values in an ADCIRC attribute file.

    Parameters:
    - outputAttrFile (str): Path to the ADCIRC attribute file to be updated.
    - slr (float): Sea level rise to be added to the sea surface height.
    - manning_node (array-like): Array of node indices for updating Manning's n.
    - manning (array-like): Array of new Manning's n values corresponding to manning_node.
    """
    print("   Processing attributes...\n")

    with open(outputAttrFile, 'r') as file:
        lines = file.readlines()
    file.close()

    skip_index = 1  # Skip the first line
    nN = int(lines[skip_index].split()[0])  # Number of nodes
    print("Number of nodes:", nN)

    idx_sshag = [i for i, line in enumerate(lines) if 'sea_surface_height_above_geoid' in line]
    idx_manning = [i for i, line in enumerate(lines) if 'mannings_n_at_sea_floor' in line]
    print('idx_sshag:\t', idx_sshag, '\nidx_manning:\t', idx_manning)

    # Sea_surface_height_above_geoid
    if idx_sshag:
        global_SSH = float(lines[idx_sshag[0] + 3].split()[0])
        new_global_SSH = global_SSH + slr
        # Update the global sea_surface_height_above_geoid # this part differs
        # from raster version. Directory change global and local values
        lines[idx_sshag[0] + 3] = "{0:.6f}\n".format(new_global_SSH)
        print("Updated global SSH:\t", new_global_SSH)

        # Read local sea_surface_height_above_geoid
        # number of local sea_surface_height_above_geoid
        num_local_sshag = int(lines[idx_sshag[1] + 1].split()[0])
        print(
            "num_local_sshag:\t",
            num_local_sshag)

        # Update local sea_surface_height_above_geoid values
        num_local_sshag = int(lines[idx_sshag[1] + 1].split()[0])
        if num_local_sshag != 0:
            for i in range(num_local_sshag):
                local_ssh_idx = int(lines[idx_sshag[1] + 2 + i].split()[0])
                local_ssh = float(lines[idx_sshag[1] + 2 + i].split()[1])
                new_local_ssh = local_ssh + slr
                # overwrite the local sea_surface_height_above_geoid
                lines[idx_sshag[1] + 2 +
                      i] = "{0:>10}     {1:.6f}\n".format(local_ssh_idx, new_local_ssh)
                lines[idx_sshag[1] + 2 + i] = "{0:>10}     {1:.6f}\n".format(local_ssh_idx, new_local_ssh)

    # Manning's n
    if idx_manning:
        # Grobal manning
        global_mann = float(lines[idx_manning[0] + 3].split()[0])
        # number of local manning's n
        num_local_manning = int(lines[idx_manning[1] + 1].split()[0])
        print(f"\nglobal manning value:\t{global_mann}\nnum_local_manning:\t{num_local_manning}")

        # Update Manning's n values
        node_value_updated = np.full(nN, global_mann, dtype=float)  # Initialize the value with global_mann
        if num_local_manning != 0:
            for i in range(num_local_manning):
                local_node, local_value = lines[idx_manning[1] + 2 + i].split()
                local_node = int(local_node)
                local_value = float(local_value)

                node_value_updated[local_node - 1] = local_value
        else:
            pass    # If there is no local manning's n, then pass

        print('size of node_value_updated:\t', len(node_value_updated))
        # Update the manning's n for the ADCIRC mesh based on the MEM results
        for i in range(len(manning_node)):
            node_value_updated[int(manning_node[i]) - 1] = manning[i]

        if num_local_manning != 0:
            # delete the local manning's n
            del lines[idx_manning[1] + 2: idx_manning[1] + 2 + num_local_manning]
        else:
            pass

        # Update the number of local manning's n
        lines[idx_manning[1] + 1] = "{0:>10}\n".format(nN)

        # Create all the lines to be inserted
        new_lines = ["{0:>10}      {1:.6f}\n".format(
            i + 1, node_value_updated[i]) for i in range(nN)]
        # Insert all the lines at once
        lines[idx_manning[1] + 2: idx_manning[1] + 2] = new_lines

    # Future add Average_horizontal_eddy_viscosity_in_sea_water_wrt_depth (Aug 28 2024)

    # Write the updated lines back to the ADCIRC attribute file
    with open(outputAttrFile, 'w') as file:
        file.writelines(lines)
    file.close()

    print(f"Updated ADCIRC file '{outputAttrFile}' with new node data from mem results")


def postprocessing_ADCIRC(inputMeshFile, inputAttrFile,
                          outputMeshFile, outputAttrFile, outputMEMFile, slr):

    # --- Initialize code ---
    start_time = time.time()
    print("\n")
    print("LAUNCH: Launching script!\n")
    ##########################################################################
    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0  # No data value (ndv) using ADCIRC conversion
    # slr = 0.0 # Sea level rise [m]

    # --- READ INPUTS ---
    # Read the mem file
    df = pd.read_csv(outputMEMFile)
    print("  Read MEM I/O file successfully")

    print(df.shape, df.columns, df.dtypes)
    new_z = df['tb_update'].values  # Update the z value for the ADCIRC mesh
    node_id = df['node'].values
    x = df['x'].values
    y = df['y'].values

    # Copy the original mesh file to the new file
    shutil.copy(inputMeshFile, outputMeshFile)
    # Copy the original attribute file to the new file
    shutil.copy(inputAttrFile, outputAttrFile)

    ####################################################################################################################
    ### Update MeshFile
    ####################################################################################################################

    update_ADCIRC_mesh(outputMeshFile, node_id, x, y, new_z)

    ####################################################################################################################
    ### Update Attribute file
    ####################################################################################################################

    # Filter out rows where 'manning' is equal to ndv
    manning_df = df[df['manning'] != ndv]

    # only keep the rows where 'manning' is not equal to ndv
    print('update manning:\t', manning_df.shape)
    manning_node = manning_df['node'].values  # update nodes
    manning = manning_df['manning'].values  # update values

    update_ADCIRC_attributes(outputAttrFile, slr, manning_node, manning)


    # # --- Read and write attributes ---
    # print("   Processing attributes...\n")
    #
    # lines = read_text_file(outputAttrFile)
    # skip_index = 1  # skip the first line
    # nN = int(lines[skip_index].split()[0])  # nN: number of nodes
    # print("number of nodes:", nN)
    #
    # # Select the line numbers that include 'sea_surface_height_above_geoid'
    # idx_sshag = [i for i, line in enumerate(
    #     lines) if 'sea_surface_height_above_geoid' in line]
    # idx_manning = [i for i, line in enumerate(
    #     lines) if 'mannings_n_at_sea_floor' in line]
    # print('idx_sshag:\t', idx_sshag, '\nidx_manning:\t', idx_manning)
    #
    # # Read the global sea_surface_height_above_geoid and manning's n
    # # global sea_surface_height_above_geoid
    # global_SSH = float(lines[idx_sshag[0] + 3].split()[0])
    # # global manning's n
    # global_mann = float(lines[idx_manning[0] + 3].split()[0])
    # print("global SSH:\t", global_SSH, "\nglobal manning's n:\t", global_mann)
    #
    # # Update the sea_surface_height_above_geoid and manning's n for the ADCIRC
    # # mesh
    # new_global_SSH = global_SSH + slr
    # # Update the global sea_surface_height_above_geoid # this part differs
    # # from raster version. Directory change global and local values
    # lines[idx_sshag[0] + 3] = "{0:.6f}\n".format(new_global_SSH)
    # print("Updated global SSH:\t", new_global_SSH)
    #
    # # Read local sea_surface_height_above_geoid and manning's n and update them
    # # number of local sea_surface_height_above_geoid
    # num_local_sshag = int(lines[idx_sshag[1] + 1].split()[0])
    # # number of local manning's n
    # num_local_manning = int(lines[idx_manning[1] + 1].split()[0])
    # print(
    #     "num_local_sshag:\t",
    #     num_local_sshag,
    #     "\nnum_local_manning:\t",
    #     num_local_manning)
    # if num_local_sshag != 0:
    #     for i in range(num_local_sshag):
    #         local_ssh_idx = int(lines[idx_sshag[1] + 2 + i].split()[0])
    #         local_ssh = float(lines[idx_sshag[1] + 2 + i].split()[1])
    #         new_local_ssh = local_ssh + slr
    #         # overwrite the local sea_surface_height_above_geoid
    #         lines[idx_sshag[1] + 2 +
    #               i] = "{0:>10}     {1:.6f}\n".format(local_ssh_idx, new_local_ssh)
    #
    # node_value_updated = np.full(nN, global_SSH, dtype=float)
    #
    # if num_local_manning != 0:
    #     for i in range(num_local_manning):
    #         local_node, local_value = lines[idx_manning[1] + 2 + i].split()
    #         local_node = int(local_node)
    #         local_value = float(local_value)
    #
    #         node_value_updated[local_node - 1] = local_value
    # else:
    #     pass    # If there is no local manning's n, then pass
    #
    # print('size of node_value_updated:\t', len(node_value_updated))
    # # Update the manning's n for the ADCIRC mesh based on the MEM results
    # for i in range(len(manning_node)):
    #     node_value_updated[int(manning_node[i]) - 1] = manning[i]
    #
    # if num_local_manning != 0:
    #     # delete the local manning's n
    #     del lines[idx_manning[1] + 2: idx_manning[1] + 2 + num_local_manning]
    # else:
    #     pass
    #
    # # Update the number of local manning's n
    # lines[idx_manning[1] + 1] = "{0:>10}\n".format(nN)
    #
    # # Create all the lines to be inserted
    # new_lines = ["{0:>10}      {1:.6f}\n".format(
    #     i + 1, node_value_updated[i]) for i in range(nN)]
    # # Insert all the lines at once
    # lines[idx_manning[1] + 2: idx_manning[1] + 2] = new_lines
    #
    # # Write the updated lines back to the ADCIRC attribute file
    # with open(outputAttrFile, 'w') as f:
    #     f.writelines(lines)
    #
    # print(
    #     f"Updated ADCIRC file '{outputAttrFile}' with new node data from mem results")

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Done interpolating tidal datums using IDW")
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")
    print("Job Finished ʕ •ᴥ•ʔ")
