#!/usr/bin/python3
# File: postprocessing.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: Jul 15, 2024

##########################################################################
# --- Load internal modules ---
from .general_functions import *


def postprocessing_ADCIRC(inputMeshFile, inputAttrFile,
                          outputMeshFile, outputAttrFile, outputMEMFile, slr):

    # --- Initialize code ---
    start_time = time.time()
    print("\n")
    print("LAUNCH: Launching script!\n")
    ##########################################################################
    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0  # No data value (ndv) using ADCIRC conversion
    # slr = 0.5 # Sea level rise [m]

    # --- READ INPUTS ---
    # Read the mem file
    df = pd.read_csv(outputMEMFile)
    print("  Read MEM I/O file successfully")

    print(df.shape, df.columns, df.dtypes)
    new_z = df['tb_update'].values  # Update the z value for the ADCIRC mesh
    node_id = df['node'].values
    x = df['x'].values
    y = df['y'].values

    # Filter out rows where 'manning' is equal to ndv
    manning_df = df[df['manning'] != ndv]

    # only keep the rows where 'manning' is not equal to ndv
    print('update manning:\t', manning_df.shape)
    manning_node = manning_df['node'].values
    manning = manning_df['manning'].values
    # Convert manning_id and manning to a dictionary
    manning_dict = dict(zip(manning_node, manning))

    # Copy the original mesh file to the new file
    shutil.copy(inputMeshFile, outputMeshFile)
    # Copy the original attribute file to the new file
    shutil.copy(inputAttrFile, outputAttrFile)

    lines = read_text_file(outputMeshFile)

    # Update the z value for the ADCIRC mesh
    for i in range(df.shape[0]):
        idx = int(node_id[i])
        # Note: z value is negated to match the ADCIRC convention. Also node
        # start from line 2 to nN-1
        lines[idx + 1] = "{0:>10}     {1:.6f}      {2:.6f} {3:.8E}\n".format(
            idx, x[i], y[i], - new_z[i])

    # Write the updated lines back to the ADCIRC file
    with open(outputMeshFile, 'w') as f:
        f.writelines(lines)

    print(
        f"Updated ADCIRC file '{outputMeshFile}' with new node data from mem results")

    # --- Read and write attributes ---
    print("   Processing attributes...\n")

    lines = read_text_file(outputAttrFile)
    skip_index = 1  # skip the first line
    nN = int(lines[skip_index].split()[0])  # nN: number of nodes
    print("number of nodes:", nN)

    # Select the line numbers that include 'sea_surface_height_above_geoid'
    idx_sshag = [i for i, line in enumerate(
        lines) if 'sea_surface_height_above_geoid' in line]
    idx_manning = [i for i, line in enumerate(
        lines) if 'mannings_n_at_sea_floor' in line]
    print('idx_sshag:\t', idx_sshag, '\nidx_manning:\t', idx_manning)

    # Read the global sea_surface_height_above_geoid and manning's n
    # global sea_surface_height_above_geoid
    global_SSH = float(lines[idx_sshag[0] + 3].split()[0])
    # global manning's n
    global_mann = float(lines[idx_manning[0] + 3].split()[0])
    print("global SSH:\t", global_SSH, "\nglobal manning's n:\t", global_mann)

    # Update the sea_surface_height_above_geoid and manning's n for the ADCIRC
    # mesh
    new_global_SSH = global_SSH + slr
    # Update the global sea_surface_height_above_geoid # this part differs
    # from raster version. Directory change global and local values
    lines[idx_sshag[0] + 3] = "{0:.6f}\n".format(new_global_SSH)
    print("Updated global SSH:\t", new_global_SSH)

    # Read local sea_surface_height_above_geoid and manning's n and update them
    # number of local sea_surface_height_above_geoid
    num_local_sshag = int(lines[idx_sshag[1] + 1].split()[0])
    # number of local manning's n
    num_local_manning = int(lines[idx_manning[1] + 1].split()[0])
    print(
        "num_local_sshag:\t",
        num_local_sshag,
        "\nnum_local_manning:\t",
        num_local_manning)
    if num_local_sshag != 0:
        for i in range(num_local_sshag):
            local_ssh_idx = int(lines[idx_sshag[1] + 2 + i].split()[0])
            local_ssh = float(lines[idx_sshag[1] + 2 + i].split()[1])
            new_local_ssh = local_ssh + slr
            # overwrite the local sea_surface_height_above_geoid
            lines[idx_sshag[1] + 2 +
                  i] = "{0:>10}     {1:.6f}\n".format(local_ssh_idx, new_local_ssh)

    node_value_updated = np.full(nN, global_SSH, dtype=float)

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

    # Write the updated lines back to the ADCIRC attribute file
    with open(outputAttrFile, 'w') as f:
        f.writelines(lines)

    print(
        f"Updated ADCIRC file '{outputAttrFile}' with new node data from mem results")

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Done interpolating tidal datums using IDW")
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")
    print("Job Finished ʕ •ᴥ•ʔ")
