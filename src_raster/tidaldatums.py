#!/usr/bin/python3
# File: tidalDatums.py
# Developer:  Peter Bacopoulos & Jin Ikeda
# Last modified: Aug 22, 2024

# ----------------------------------------------------------
# F U N C T I O N    T I D A L D A T U M S
# ----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW),
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = function(inputHarmonicsReadMe,inputRasterHarmonics,inputRasterHyControl,outputRasterControl,tstep=1800.0)
# tstep (OPTIONAL 900-3600 [sec] would be good, DEFAULT = 1800 SECONDS)
# ----------------------------------------------------------

# Command line: python tidalDatums.py -i harmonics.README -s harmonics.img
# -d HyControl.img -o tidalDatums.img

##########################################################################
# --- Load internal modules ---
from osgeo import gdal
from osgeo import osr
import numpy as np
import src_raster.basics
import time

def astronomic_tide_resynthesis(time, omega, phase, amplitude):
    nc = len(amplitude)
    signal = 0
    for j in range(0, nc):
        signal += amplitude[j] * np.cos(omega[j] * time - phase[j])
    return signal


# Thanks for the help, ChatGPT
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


# Thanks for the help, ChatGPT
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


def tidaldatums(inputHarmonicsReadMe, inputRasterHarmonics,
                inputRasterHyControl, outputRasterControl, tstep=1800.0):

    start_time = time.time()

    # ----------------------------------------------------------
    # Check that input files exist
    # ----------------------------------------------------------
    src_raster.basics.fileexists(inputHarmonicsReadMe)
    src_raster.basics.fileexists(inputRasterHarmonics)
    src_raster.basics.fileexists(inputRasterHyControl)

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

    # Read Harmonic data
    harmonics = dict()
    f = open(inputHarmonicsReadMe, 'r')
    for line in f:
        fields = line.split(':', 1)
        try:
            harmonics[fields[0].strip()] = fields[1].strip()
        except IndexError:
            continue
    f.close()
    #print ("  Dictionary of harmonics (README) read successfully")
    nc = int(harmonics['number.constituents'])
    #print ("    Number of constituents:",nc);

    #print ("  Writing frequency information to internal memory...")
    omega = tstep * np.zeros((nc, 1), dtype=float)
    j = 0
    for k, v in harmonics.items():
        j = j + 1
        if j > 1:
            omega[j - 3][0] = v
    #print ("    DONE!")

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
    valid_rows, valid_cols = np.where((hc >= 1.5) & (hc <= 2.5))

    print("Size of t:", t.shape)
    print("Size of omega:", omega.shape)
    print("Size of pha:", pha.shape)
    print("Size of amp:", amp.shape)
    print("Size of wl:", valid_rows.shape)

    # Initialize mlw, msl, and mhw arrays with NoData values
    mlw = np.full_like(hc, ndv, dtype=float)
    msl = np.full_like(hc, ndv, dtype=float)
    mhw = np.full_like(hc, ndv, dtype=float)

    # Calculate mlw, msl, and mhw only for valid hc values
    for kk, k in zip(valid_rows, valid_cols):
        # each row and column needs to calculate time series of astnomic tides
        # (this part is time consuming)
        wl = astronomic_tide_resynthesis(
            t, omega, pha[kk, k, :], amp[kk, k, :])

        mlw[kk, k] = mean_low_water(wl)
        msl[kk, k] = np.average(wl)
        mhw[kk, k] = mean_high_water(wl)

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