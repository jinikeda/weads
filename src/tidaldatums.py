#!/usr/bin/python3
# File: tidalDatums.py
# Name: Peter Bacopoulos
# Date: August 16, 2021
# Command line: python tidalDatums.py -i harmonics.README -s harmonics.img -d HyControl.img -o tidalDatums.img
from osgeo import gdal
from osgeo import osr
import numpy as np
import src.basics

#----------------------------------------------------------

def astronomic_tide_resynthesis(time, omega, phase, amplitude):
    nc = len(amplitude)
    signal = 0
    for j in range(0,nc):
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


#----------------------------------------------------------
# F U N C T I O N    T I D A L D A T U M S           
#----------------------------------------------------------
#
# Computes tidal datums of mean low water (MLW), 
# mean sea level (MSL), and mean high water (MHW)
# for all nodes hydraulically connected to the ocean.
# result = function(inputHarmonicsReadMe,inputRasterHarmonics,inputRasterHyControl,outputRasterControl,tstep=900.0)
# tstep (OPTIONAL, DEFAULT = 900 SECONDS)
#----------------------------------------------------------
def tidaldatums(inputHarmonicsReadMe,inputRasterHarmonics,\
    inputRasterHyControl,outputRasterControl,tstep=900.0):

    #----------------------------------------------------------
    # Check that input files exist
    #----------------------------------------------------------
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
    ndv=-99999.0 # No data value (ndv) using ADCIRC convention
    tdfin=15.0   # Total length of resynthesis in units of days

    # --- READ INPUTS ---
    #print ("")
    #print ("Reading rasters")

    htdfin=int(tdfin*2)
    nstep=int((86400*tdfin/tstep)+1)
    dww=int(24*3600/tstep)
    hdww=int(dww/2)
    t=tstep*np.ones((nstep,1),dtype=float); t=np.cumsum(t,axis=0);

    # --- INITIALIZE ARRAYS ---
    hdhlwt=np.zeros((hdww,1),dtype=float)
    hdhlwu=np.zeros((int(htdfin),1),dtype=float)
    hdhlwd=np.zeros((int(htdfin),1),dtype=float)

    harmonics=dict()
    f=open(inputHarmonicsReadMe,'r')
    for line in f:
        fields=line.split(':',1)
        try:
            harmonics[fields[0].strip()]=fields[1].strip()
        except IndexError:
            continue
    f.close()
    #print ("  Dictionary of harmonics (README) read successfully")
    nc=int(harmonics['number.constituents'])
    #print ("    Number of constituents:",nc);

    #print ("  Writing frequency information to internal memory...")
    omega=tstep*np.zeros((nc,1),dtype=float); j=0;
    for k, v in harmonics.items():
        j=j+1
        if j>1:
            omega[j-3][0]=v
    #print ("    DONE!")

    rasterHARM=gdal.Open(inputRasterHarmonics)
    gt=rasterHARM.GetGeoTransform(); pixelSizeX=gt[1]; pixelSizeY=-gt[5];

    rasterHC=gdal.Open(inputRasterHyControl)
    gt=rasterHC.GetGeoTransform()
    pixelSizeX=gt[1]
    pixelSizeY=-gt[5];
         
    band=rasterHC.GetRasterBand(1); hc=band.ReadAsArray();
    amp=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize,nc),dtype=float)
    pha=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize,nc),dtype=float)
    for j in range(0,nc):
        band=rasterHARM.GetRasterBand(2*j+1)
        amp[:,:,j]=band.ReadAsArray();
        
        band=rasterHARM.GetRasterBand(2*j+1+1)
        pha[:,:,j]=band.ReadAsArray();
    
    pha=pha*np.pi/180.0

    mlw = np.ones((rasterHC.RasterYSize, rasterHC.RasterXSize), dtype=float) * ndv
    msl = np.ones((rasterHC.RasterYSize, rasterHC.RasterXSize), dtype=float) * ndv
    mhw = np.ones((rasterHC.RasterYSize, rasterHC.RasterXSize), dtype=float) * ndv

###########################################################################################################
# This part has to speed up (Jin -> Pete)
    for k in range(0, rasterHC.RasterXSize - 1):
        for kk in range(0, rasterHC.RasterYSize - 1):
            # Check if hc value is within a certain range
            if hc[kk][k] < 0.5 and hc[kk][k] > -0.5:
                wl = astronomic_tide_resynthesis(t, omega, pha[kk][k][:], amp[kk][k][:])
    
                # Calculate the mean low water, mean sea level, and mean high water
                mlw[kk][k] = mean_low_water(wl)
                msl[kk][k] = np.average(wl)
                mhw[kk][k] = mean_high_water(wl)

#######################################################################################
# ######here is the chatGPT suggestion (Please confirm it)

# # Create a mask based on the condition on hc values
# mask = (hc < 0.5) & (hc > -0.5)

# # Calculate wl for the entire array
# ## wl need background? If so, make an array list with 0 or np.nan (Jin's comment)
# wl = astronomic_tide_resynthesis(t, omega, pha[mask], amp[mask])

# ##############################################################
# # Especially please check this part (Jin's comments) #

# # Assign calculated values to mlw and mhw arrays using the mask
# mlw[mask] = mean_low_water(wl)
# mhw[mask]  = mean_high_water(wl)
        
###########################################################################################################                
                
        # Check for bogus values and replace with NoData
        # This can happen on the edges based on the raster resolution
        mlw[(mlw < -99) | (mlw > 99)] = ndv
        msl[(msl < -99) | (msl > 99)] = ndv
        mhw[(mhw < -99) | (mhw > 99)] = ndv
    

    # --- WRITE OUTPUTS ---
    #print ("")
    #print ("Writing output raster")
    driver=gdal.GetDriverByName('HFA'); dst_datatype=gdal.GDT_Float32;
    dst_geot=rasterHC.GetGeoTransform(); dst_proj=osr.SpatialReference(); dst_proj.ImportFromWkt(rasterHC.GetProjectionRef());
    dst_ds=driver.Create(outputRasterControl,rasterHC.RasterXSize,rasterHC.RasterYSize,3,dst_datatype)
    dst_ds.SetGeoTransform(dst_geot); dst_ds.SetProjection(dst_proj.ExportToWkt());
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv); dst_ds.GetRasterBand(1).WriteArray(mlw);
    dst_ds.GetRasterBand(2).SetNoDataValue(ndv); dst_ds.GetRasterBand(2).WriteArray(msl);
    dst_ds.GetRasterBand(3).SetNoDataValue(ndv); dst_ds.GetRasterBand(3).WriteArray(mhw);
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

