#!/usr/bin/python3
# File: tidalDatums.py
# Name: Peter Bacopoulos, Jin Ikeda
# Last modified: March 8, 2024
# Command line: python tidalDatums.py -i harmonics.README -s harmonics.img
# -d HyControl.img -o tidalDatums.img
from osgeo import gdal
from osgeo import osr
import numpy as np
import src.basics

# import pybind11_functions as pybind # This is the C++ code that is
# called by the Python code
import pybind11_functions as pyb

# ----------------------------------------------------------

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


# ----------------------------------------------------------
# F U N C T I O N    T I D A L D A T U M S
# ----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW),
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = function(inputHarmonicsReadMe,inputRasterHarmonics,inputRasterHyControl,outputRasterControl,tstep=900.0)
# tstep (OPTIONAL, DEFAULT = 900 SECONDS)
# ----------------------------------------------------------
def tidaldatums(inputHarmonicsReadMe, inputRasterHarmonics,
                inputRasterHyControl, outputRasterControl, tstep=900.0):
    # ----------------------------------------------------------
    # Check that input files exist
    # ----------------------------------------------------------
    src.basics.fileexists(inputHarmonicsReadMe)
    src.basics.fileexists(inputRasterHarmonics)
    src.basics.fileexists(inputRasterHyControl)

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
    tdfin = 15.0   # Total length of resynthesis in units of days

    # --- READ INPUTS ---
    #print ("")
    #print ("Reading rasters")

    # htdfin=int(tdfin*2)
    nstep = int((86400 * tdfin / tstep) + 1)
    # dww=int(24*3600/tstep)
    # hdww=int(dww/2)
    t = tstep * np.ones((nstep, 1), dtype=float)
    t = np.cumsum(t, axis=0)

    # --- INITIALIZE ARRAYS ---
    # hdhlwt=np.zeros((hdww,1),dtype=float)
    # hdhlwu=np.zeros((int(htdfin),1),dtype=float)
    # hdhlwd=np.zeros((int(htdfin),1),dtype=float)

    # Open the file and read the lines
    with open(inputHarmonicsReadMe, 'r') as f:
        lines = f.readlines()

    # Initialize an empty dictionary
    harmonics = {}

    # Split each line into a key-value pair and add it to the dictionary
    for line in lines:
        key, value = line.split(':', 1)
        harmonics[key.strip()] = value.strip()

    # Get the number of constituents
    nc = int(harmonics['number.constituents'])
    print("Number of constituents:", nc)

    # Initialize a 1D numpy array with zeros
    omega = np.zeros(nc)

    # Fill the omega array with values from the harmonics dictionary
    for j, (k, v) in enumerate(harmonics.items(), start=1):
        if j > 1:
            omega[j - 3] = v

    # harmonics=dict()
    # f=open(inputHarmonicsReadMe,'r')
    # for line in f:
    #     fields=line.split(':',1)
    #     try:
    #         harmonics[fields[0].strip()]=fields[1].strip()
    #     except IndexError:
    #         continue
    # f.close()
    # #print ("  Dictionary of harmonics (README) read successfully")
    # nc=int(harmonics['number.constituents'])
    # print ("    Number of constituents:",nc);
    #
    # #print ("  Writing frequency information to internal memory...")
    # omega=tstep*np.zeros((nc,1),dtype=float)
    # j=0
    # for k, v in harmonics.items():
    #     j=j+1
    #     if j>1:
    #         omega[j-3][0]=v
    #print ("    DONE!")
    print(f"Size of omega: {omega.shape}, the values\t ", omega)

    rasterHARM = gdal.Open(inputRasterHarmonics)
    rasterHC = gdal.Open(inputRasterHyControl)

    band = rasterHC.GetRasterBand(1)
    hc = band.ReadAsArray()
    amp = np.zeros(
        (rasterHC.RasterYSize,
         rasterHC.RasterXSize,
         nc),
        dtype=float)
    pha = np.zeros(
        (rasterHC.RasterYSize,
         rasterHC.RasterXSize,
         nc),
        dtype=float)
    for j in range(0, nc):
        band = rasterHARM.GetRasterBand(2 * j + 1)
        amp[:, :, j] = band.ReadAsArray()

        band = rasterHARM.GetRasterBand(2 * j + 1 + 1)
        pha[:, :, j] = band.ReadAsArray()

    pha = pha * np.pi / 180.0

    # Stored in valid_indices of water regions (valid_rows, valid_cols) #
    # water regions are 2.0
    print(f"Size of t: {t.shape}, dtype: {t.dtype}")
    print(f"Size of omega: {omega.shape}, dtype: {omega.dtype}")
    print(f"Size of pha: {pha.shape}, dtype: {pha.dtype}")
    print(f"Size of amp: {amp.shape}, dtype: {amp.dtype}")

    # # Find valid indices of water regions
    # valid_rows, valid_cols = np.where((hc >= 1.5) & (hc <= 2.5))
    # print(f"Size of wl: {valid_rows.shape}, dtype: {valid_rows.dtype}")
    #
    # # Initialize mlw, msl, and mhw arrays with NoData values
    # mlw = np.full_like(hc, ndv, dtype=float)
    # msl = np.full_like(hc, ndv, dtype=float)
    # mhw = np.full_like(hc, ndv, dtype=float)
    #
    # # Calculate mlw, msl, and mhw only for valid hc values
    # for kk, k in zip(valid_rows, valid_cols):
    #     wl = astronomic_tide_resynthesis(t, omega, pha[kk, k, :], amp[kk, k, :]) # each row and column needs to calculate time series of astnomic tides (this part is time consuming)
    #     #wl = pybind.astronomic_tide_resynthesis(t, omega, pha[kk, k, :], amp[kk, k, :]) # This is the C++ code that is called by the Python code
    #     mlw[kk, k] = mean_low_water(wl)
    #     msl[kk, k] = np.average(wl)
    #     mhw[kk, k] = mean_high_water(wl)
    # print ("wl.shape", wl.shape)

    # --- CALCULATE TIDAL DATUMS ---
    # Call the C++ function
    ans = pyb.add(1, 2)
    print(f"ans: {ans}")
    varid_index = pyb.find_valid_indices(hc, hc.shape[0], hc.shape[1])
    # Convert varid_index to a numpy array
    varid_index_np = np.array(varid_index)
    valid_rows, valid_cols = np.where((hc >= 1.5) & (hc <= 2.5))
    print(f"Size of wl: {valid_rows.shape}, dtype: {valid_rows.dtype}")

    print(hc.shape[0], hc.shape[1], varid_index_np.shape[0])
    # Call the function
    mlw, msl, mhw = pyb.compute_mlmsmhw(hc, t, omega, pha, amp, ndv)

    # Print the results (replace this with your desired further processing)
    print("mlw:")
    print(mlw.min(), mlw.max())
    print("mlw size:", mlw.shape)
    print("msl:")
    print(msl.min(), msl.max())
    print("msl size:", msl.shape)
    print("mhw:")
    print(mhw.min(), mhw.max())
    print("mhw size:", mhw.shape)

    # Check for bogus values and replace with NoData
    # This can happen on the edges based on the raster resolution
    mlw[(mlw < -99) | (mlw > 99)] = ndv
    msl[(msl < -99) | (msl > 99)] = ndv
    mhw[(mhw < -99) | (mhw > 99)] = ndv

    # --- WRITE OUTPUTS ---
    #print ("")
    #print ("Writing output raster")
    driver = gdal.GetDriverByName('HFA')
    dst_datatype = gdal.GDT_Float32
    dst_geot = rasterHC.GetGeoTransform()
    dst_proj = osr.SpatialReference()
    dst_proj.ImportFromWkt(rasterHC.GetProjectionRef())
    dst_ds = driver.Create(
        outputRasterControl,
        rasterHC.RasterXSize,
        rasterHC.RasterYSize,
        3,
        dst_datatype)
    dst_ds.SetGeoTransform(dst_geot)
    dst_ds.SetProjection(dst_proj.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(1).WriteArray(mlw)
    dst_ds.GetRasterBand(2).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(2).WriteArray(msl)
    dst_ds.GetRasterBand(3).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(3).WriteArray(mhw)
    #print ("  Output raster written successfully")
    #print ("")

    # --- PRINT TO SCREEN ---
    '''
    print ("")
    print ("----------------------------------------------------")
    print ("--- COMPLETED HYDRO-MEM IN PYTHON - TIDAL DATUMS ---")
    print ("----------------------------------------------------")
    print ("")
    print ("")

    print ("")
    print ("")
    '''
