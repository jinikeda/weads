#!/usr/bin/python3
# File: hydrodynamics.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: Jul 9, 2024

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

########################################################################################################################
#--- Load internal modules ---
from general_functions import *

#--- Initialize code ---
start_time = time.time()
print("\n"); print("LAUNCH: Launching script!\n")

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
lines = read_text_file(inputHarmonicsFile)

numHarm = len(lines)
freq = np.zeros((numHarm,1),dtype=float)
for i in range(numHarm):
    freq[i] = float(lines[i].split()[0])
print (freq)

#--- Calculate tidal datums for wetted zones ---
print("\n"); print("   Calculating tidal datums for wetted zones...\n")

numPoints = df.shape[0] # Number of points within the domain

T = 30.0 # Tidal period [days]
dt = 3600.0 # Time step [seconds]
N = int((86400*T/dt)+1) # Number of time steps
t = np.cumsum(dt * np.ones((N,1), dtype=float), axis=0) # Time vector
T2 = 2.0*T
D = 24.0*3600.0/dt
D2 = D/2.0
wl2 = np.ones((int(D2), 1), dtype=float)
wl2L = np.ones((int(T2), 1), dtype=float)
wl2H = np.ones((int(T2), 1), dtype=float)
mlw = np.ones((numPoints, 1), dtype=float)
mhw = np.ones((numPoints, 1), dtype=float)
msl = np.ones((numPoints, 1), dtype=float)

inundationtime_index = df.columns.get_loc('inundationtime') # Get the column index for 'inundationtime'
amp_index = df.columns.get_loc('STEADY_amp') # Get the column index for 'STEADY_amp'
phase_index = df.columns.get_loc('STEADY_phase') # Get the column index for 'STEADY_phase'
print ('amp_index:', amp_index, 'phase_index:', phase_index)

for j in range(numPoints):
    if (df.iloc[j,inundationtime_index] >= 0.9999) & (df.iloc[j,inundationtime_index] <= 1.0001): # If inundationTime is 1.0, then the point is fully wetted
        # Tidal resynthesis
        wl = df.iloc[j, amp_index] * np.ones((N, 1), dtype=float)
        for count, jj in enumerate(range(amp_index, amp_index + numHarm)):
            wl += df.iloc[j, jj] * np.cos(freq[count] * t - df.iloc[j,jj + numHarm])

### Jin, July 3, 2024 Need to modify the code below (cannot remember the reasons)
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
df.to_csv("domain_tide.csv", index=False)

########################################################################################################################
# Calculate the elapsed time
end_time = time.time()
elapsed_time = end_time - start_time

# Print the elapsed time
print("Done water level calculation")
print("Time to Compute: \t\t\t", elapsed_time, " seconds")
print("Job Finished ʕ •ᴥ•ʔ")


