#!/usr/bin/python3
# File: hydrodynamics.py
# Name: Jin Ikeda and Peter Bacopoulos
# Date: July 3, 2024

#--- Load modules ---
import numpy as np
import time
import pandas as pd

# def astronomic_tide_resynthesis(time, omega, phase, amplitude):
#     nc = len(amplitude)
#     signal = 0
#     for j in range(0, nc):
#         signal += amplitude[j] * np.cos(omega[j] * time - phase[j])
#     return signal
#
#
# # Thanks for the help, ChatGPT
# def mean_high_water(data):
#     """
#     Calculates the mean high water (MHW) from a time series of water level data.
#     :param data: 1D numpy array of water level data.
#     :return: MHW, a scalar value representing the mean high water.
#     """
#     # Find the indices of local maxima in the data
#     max_indices = (data[1:-1] > data[:-2]) & (data[1:-1] > data[2:])
#     max_indices = np.where(max_indices)[0] + 1
#
#     # Extract the values at the local maxima
#     max_values = data[max_indices]
#
#     # Calculate the mean of the maxima
#     MHW = np.mean(max_values)
#
#     return MHW
#
#
# # Thanks for the help, ChatGPT
# def mean_low_water(data):
#     """
#     Calculates the mean low water (MLW) from a time series of water level data.
#     :param data: 1D numpy array of water level data.
#     :return: MLW, a scalar value representing the mean low water.
#     """
#     # Find the indices of local minima in the data
#     min_indices = (data[1:-1] < data[:-2]) & (data[1:-1] < data[2:])
#     min_indices = np.where(min_indices)[0] + 1
#
#     # Extract the values at the local minima
#     min_values = data[min_indices]
#
#     # Calculate the mean of the minima
#     MLW = np.mean(min_values)
#
#     return MLW

def read_text_file(fileName):
    with open(fileName, "r") as f:
        lines = f.readlines()
    return lines

# ----------------------------------------------------------
# F U N C T I O N    T I D A L D A T U M S
# ----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW),
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = function(inputHarmonicsReadMe,inputRasterHarmonics,inputRasterHyControl,outputRasterControl,tstep=900.0)

#### need to modify the function Jin July 3, 2024

# tstep (OPTIONAL, DEFAULT = 900 SECONDS)
# ----------------------------------------------------------

#--- Initialize code ---
start=time.time()
print("\n"); print("LAUNCH: Launching script!\n");

########################################################################################################################
#--- Read harmonics (inputHarmonicsFile) ---
########################################################################################################################
print("   Processing scatter points...\n")
# Read the CSV file
df = pd.read_csv("domain_inputs.csv")

# Print the DataFrame
print(df.columns)
print(df.shape)

#--- Read harWmomics correspondent with scatter ---
print("\n"); print("   Processing harmonics...\n");
inputHarmonicsFile = "harmonics_Freq.txt" # Read the harmonics file
lines = read_text_file (inputHarmonicsFile)

numHarm = len(lines)
freq=np.zeros((numHarm,1),dtype=float)
for i in range(numHarm):
    freq[i]=float(lines[i].split()[0])
print (freq)

#--- Calculate tidal datums for wetted zones ---
print("\n"); print("   Calculating tidal datums for wetted zones...\n")

numPoints = df.shape[0] # Number of points within the domain

T=30.0;
dt=3600.0;
N=int((86400*T/dt)+1);
t=dt*np.ones((N,1),dtype=float);
t=np.cumsum(t,axis=0);
T2=2.0*T;
D=24.0*3600.0/dt;
D2=D/2.0;
wl2 = np.ones((int(D2), 1), dtype=float)
wl2L = np.ones((int(T2), 1), dtype=float)
wl2H = np.ones((int(T2), 1), dtype=float)
mlw = np.ones((numPoints, 1), dtype=float)
mhw = np.ones((numPoints, 1), dtype=float)
msl = np.ones((numPoints, 1), dtype=float)

amp_index = df.columns.get_loc('STEADY_amp') # Get the column index for 'STEADY_amp'
phase_index = df.columns.get_loc('STEADY_phase') # Get the column index for 'STEADY_phase'
print ('amp_index:', amp_index, 'phase_index:', phase_index)

for j in range(numPoints):
    if df.iloc[j,amp_index-1] > 0.0:
        # Tidal resynthesis
        wl = df.iloc[j, amp_index] * np.ones((N, 1), dtype=float)
        for count, jj in enumerate(range(amp_index, amp_index + numHarm)):
            wl += df.iloc[j, jj] * np.cos(freq[count] * t - df.iloc[j,jj + numHarm])

### Jin, July 3, 2024 Need to modify the code below
########################################################################################################################
        # Tidal datums
        wl2 = wl[int(D2) * np.arange(int(T2)) + np.arange(int(D2))[:, None]]
        wl2L = np.min(wl2, axis=0)
        wl2H = np.max(wl2, axis=0)
########################################################################################################################

        msl[j, 0] = np.mean(wl)
        mlw[j, 0] = np.mean(wl2L)
        mhw[j, 0] = np.mean(wl2H)
    else:
        msl[j, 0] = -99999.0
        mlw[j, 0] = -99999.0
        mhw[j, 0] = -99999.0

# Add 'msl', 'mlw', and 'mhw' to the DataFrame
df['msl'] = msl.flatten()
df['mlw'] = mlw.flatten()
df['mhw'] = mhw.flatten()

# Save the DataFrame to a CSV file
df.to_csv("domain_tide.csv", index=True)

#--- Exit script ---
print("\n"); print("EXIT: Existing script!\n");
end=time.time(); print ("Time elapsed (seconds):",end-start);


