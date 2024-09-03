#!/usr/bin/python3
# File: tidalDatums.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: Aug 22, 2024

# ----------------------------------------------------------
# F U N C T I O N    T I D A L D A T U M S
# ----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW),
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = function(domainIOFile,
#             inputHarmonicFreqFile,
#             outputHarmonicsFile,
#             tstep=1800.0)
# tstep (OPTIONAL 900-3600 [sec] would be good, DEFAULT = 1800 SECONDS)

##########################################################################
# --- Load internal modules ---
from .general_functions import *


def mean_high_water(data):
    """
    Calculates the mean high water (MHW) from a time series of water level data.
    :param data: 1D numpy array of water level data.
    :return: MHW, a scalar value representing the mean high water.
    """
    # Find the indices of local maxima in the data
    max_indices = (data[1:-1] > data[:-2]) & (data[1:-1] > data[2:])
    max_indices = np.where(max_indices)[0] + 1

    # Extract the values at the local maxima
    max_values = data[max_indices]

    # Calculate the mean of the maxima
    MHW = np.mean(max_values)

    return MHW


def mean_low_water(data):
    """
    Calculates the mean low water (MLW) from a time series of water level data.
    :param data: 1D numpy array of water level data.
    :return: MLW, a scalar value representing the mean low water.
    """
    # Find the indices of local minima in the data
    min_indices = (data[1:-1] < data[:-2]) & (data[1:-1] < data[2:])
    min_indices = np.where(min_indices)[0] + 1

    # Extract the values at the local minima
    min_values = data[min_indices]

    # Calculate the mean of the minima
    MLW = np.mean(min_values)

    return MLW


def tidaldatums(domainIOFile, inputHarmonicsFile,
                outputHarmonicsFile, tstep=1800.0):

    start_time = time.time()

    '''
    print ("")
    print ("")
    print ("----------------------------------------------------")
    print ("--- LAUNCHING HYDRO-MEM IN PYTHON - TIDAL DATUMS ---")
    print ("----------------------------------------------------")
    print ("")
    '''

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0  # No data value (ndv) using ADCIRC convention
    T = 30.0  # Tidal period [days]
    dt = tstep  # Time step [seconds]
    N = int((86400 * T / dt) + 1)  # Number of time steps
    t = np.cumsum(dt * np.ones((N, 1), dtype=float), axis=0)  # Time vector

    ##########################################################################
    # --- Read harmonics (inputHarmonicsFile) ---
    ##########################################################################
    print("   Processing scatter points...\n")

    # Read the CSV file
    df = pd.read_csv(domainIOFile)

    # Print the DataFrame
    print(df.columns)
    print(df.shape)
    numPoints = df.shape[0]  # Number of points within the domain

    # --- Read harmomics correspondent with scatter ---
    print("\n")
    print("   Processing harmonics...\n")
    lines = read_text_file(inputHarmonicsFile)

    numHarm = len(lines)
    freq = np.zeros((numHarm, 1), dtype=float)
    for i in range(numHarm):
        freq[i] = float(lines[i].split()[0])
    print(freq)

    # --- Calculate tidal datums for wetted zones ---
    print("\n")
    print("   Calculating tidal datums for wetted zones...\n")

    # Initialize mlw, msl, and mhw arrays with NoData values
    mlw = np.full((numPoints, 1), ndv, dtype=float)
    msl = np.full((numPoints, 1), ndv, dtype=float)
    mhw = np.full((numPoints, 1), ndv, dtype=float)

    # Get the column index for 'inundationtime', 'STEADY_amp', and 'STEADY_phase'
    inundationtime_index = df.columns.get_loc('inundationtime')
    amp_index = df.columns.get_loc('STEADY_amp')
    phase_index = df.columns.get_loc('STEADY_phase')
    print('amp_index:', amp_index, 'phase_index:', phase_index)

    for j in tqdm(range(numPoints)):
        if (df.iloc[j, inundationtime_index] >= 0.9999) & (df.iloc[j, inundationtime_index]
                                                           <= 1.0001):  # If inundationTime is 1.0, then the point is fully wetted
            # Tidal resynthesis
            wl = np.zeros((N, 1), dtype=float)
            for count, jj in enumerate(range(amp_index, amp_index + numHarm)):
                # Compute the water level a*Cos(wt - phi)
                wl += df.iloc[j,
                              jj] * np.cos(freq[count] * t - df.iloc[j,
                                           jj + numHarm] * np.pi / 180.0)

            # Calculate tidal datums
            mlw[j] = mean_low_water(wl)
            msl[j] = np.average(wl)
            mhw[j] = mean_high_water(wl)

            # Check for bogus values and replace with NoData
            # This can happen on the edges based on the raster resolution
            mhw[j][(mhw[j] < -99) | (mhw[j] > 99)] = ndv
            mhw[j][(mhw[j] < -99) | (mhw[j] > 99)] = ndv
            mhw[j][(mhw[j] < -99) | (mhw[j] > 99)] = ndv

        else:
            msl[j] = ndv
            mlw[j] = ndv
            mhw[j] = ndv

    # Add 'msl', 'mlw', and 'mhw' to the DataFrame (domain_inputs.csv)
    df['msl'] = msl.flatten()
    df['mlw'] = mlw.flatten()
    df['mhw'] = mhw.flatten()

    # Save the DataFrame to a CSV file
    df.to_csv(outputHarmonicsFile, index=False)

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # --- PRINT TO SCREEN ---
    print("\n----------------------------------------------------")
    print("------------- COMPLETED - TIDAL DATUMS -------------")
    print("----------------------------------------------------\n")
    print(f"Time to Compute: {elapsed_time} seconds")
    print("Job Finished ʕ •ᴥ•ʔ\n")
