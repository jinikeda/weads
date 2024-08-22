#!/usr/bin/python3
# File: tidalDatumsIDW.py
# Name: Peter Bacopoulos, Matthew V. Bilskie
# Date: August 20, 2021
# Command line: python tidalDatumsIDW.py -i HyControl.img -s
# tidalDatums.img -o tidalDatumsIDW.img

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------
# ----------------------------------------------------------
from osgeo import gdal
from osgeo import osr
import numpy as np
import src.raster
import lidario as lio
from scipy.spatial import KDTree


# ----------------------------------------------------------
# F U N C T I O N	T I D A L D A T U M S I D W
# ----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums outside of the hydraulically connect wet areas.
# result = function(inputRasterHyControl,inputRasterTidalDatums,\
#                   outputRaster,numNeighbors = 12)
# ----------------------------------------------------------
def tidaldatumsidw(inputRasterHyControl, inputRasterTidalDatums,
                   outputRaster, numNeighbors=12):
    '''
    print ("")
    print ("")
    print ("--------------------------------------------------------")
    print ("--- LAUNCHING HYDRO-MEM IN PYTHON - TIDAL DATUMS IDW ---")
    print ("--------------------------------------------------------")
    print ("")
    '''

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0  # No data value (ndv) using ADCIRC convention

    # --- READ INPUTS ---
    rasterHC = gdal.Open(inputRasterHyControl)
    gt = rasterHC.GetGeoTransform()
    rasterTD = gdal.Open(inputRasterTidalDatums)
    gt = rasterTD.GetGeoTransform()

    '''
    print ("")
    print ("Reading rasters")

    print ("  Raster of HyControl (HC) read successfully")
    print ("    Size (x):",rasterHC.RasterXSize)
    print ("    Size (y):",rasterHC.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterHC.RasterCount)


    print ("  Raster of tidal datums (TD) read successfully")
    print ("    Size (x):",rasterTD.RasterXSize)
    print ("    Size (y):",rasterTD.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterTD.RasterCount)
    '''

    # --- INITIALIZE ARRAYS ---
    band = rasterHC.GetRasterBand(1)
    hc = band.ReadAsArray()
    mlwIDW = np.ones((rasterTD.RasterYSize, rasterTD.RasterXSize), dtype=float)
    mlwIDW = mlwIDW * ndv
    mslIDW = np.ones((rasterTD.RasterYSize, rasterTD.RasterXSize), dtype=float)
    mslIDW = mslIDW * ndv
    mhwIDW = np.ones((rasterTD.RasterYSize, rasterTD.RasterXSize), dtype=float)
    mhwIDW = mhwIDW * ndv
    X = np.zeros((rasterHC.RasterYSize, rasterHC.RasterXSize), dtype=float)
    x = np.ones((1, rasterHC.RasterXSize), dtype=float)
    x = np.cumsum(x, axis=1)
    X = X + x
    Y = np.zeros((rasterHC.RasterYSize, rasterHC.RasterXSize), dtype=float)
    y = np.ones((rasterHC.RasterYSize, 1), dtype=float)
    y = np.cumsum(y, axis=0)
    Y = Y + y

    # --- PERFORM IDW OF TIDAL DATUMS ---
    #print ("")
    #print ("Performing IDW of tidal datums")

    # supress/hide: invalid value encountered in true_divide
    np.seterr("ignore")

    # Gathering x,y,z point cloud from the raster
    translator = lio.Translator("geotiff", "np")
    point_cloud = translator.translate(
        inputRasterTidalDatums, no_data=ndv, band=1)

    pts = np.delete(point_cloud, np.s_[2], axis=1)  # Remove Z
    mlw = np.delete(point_cloud, np.s_[0, 1], axis=1)
    msl = translator.translate(inputRasterTidalDatums, no_data=ndv, band=2)
    msl = np.delete(msl, np.s_[0, 1], axis=1)
    mhw = translator.translate(inputRasterTidalDatums, no_data=ndv, band=3)
    mhw = np.delete(mhw, np.s_[0, 1], axis=1)

    tree = KDTree(pts)

    for k in range(0, rasterHC.RasterXSize):
        for kk in range(0, rasterHC.RasterYSize):
            # Only interpolate across the LAND (0) and intertidal (1) regions
            if (-0.5 <= hc[kk][k]) & (hc[kk][k] <= 1.5):
                px, py = src.raster.pixel2coord(
                    k, kk, rasterHC)  # x,y we want to interpolate to
                # Find the nearest neighbors
                dd, ii = tree.query([px, py], numNeighbors)

                w = np.zeros(numNeighbors)
                tmlw = np.zeros(numNeighbors)
                tmsl = np.zeros(numNeighbors)
                tmhw = np.zeros(numNeighbors)
                mlw_IDW = 0.0
                msl_IDW = 0.0
                mhw_IDW = 0.0

                for i in range(0, numNeighbors):
                    # Find the col,row for the ith nearest neighbor
                    c, r = src.raster.coord2pixel(
                        pts[ii[i]][0], pts[ii[i]][1], rasterTD)

                    # Get the value
                    tmlw[i] = mlw[ii[i]]
                    tmsl[i] = msl[ii[i]]
                    tmhw[i] = mhw[ii[i]]

                    w[i] = 1.0 / (dd[i]**2)  # Inverse distance of ith neighbor
                    # w[i] = 1.0/dd[i] # Inverse distance of ith neighbor

                    mlw_IDW = mlw_IDW + w[i] * tmlw[i]
                    msl_IDW = msl_IDW + w[i] * tmsl[i]
                    mhw_IDW = mhw_IDW + w[i] * tmhw[i]

                mlwIDW[kk][k] = mlw_IDW / np.nansum(w)
                mslIDW[kk][k] = msl_IDW / np.nansum(w)
                mhwIDW[kk][k] = mhw_IDW / np.nansum(w)

##########################################################################
    # Future consideration: do we need to reduce the IDW values for land
    # regions?
    land_mask = (-0.5 <= hc) & (hc <= 0.5)
    reduction_factor = 1.0
    mlwIDW[land_mask] = mlwIDW[land_mask] * reduction_factor
    mslIDW[land_mask] = mslIDW[land_mask] * reduction_factor
    mhwIDW[land_mask] = mhwIDW[land_mask] * reduction_factor
##########################################################################

    # --- WRITE OUTPUTS ---
    #print ("")
    #print ("Writing output raster")
    driver = gdal.GetDriverByName('HFA')
    dst_datatype = gdal.GDT_Float32
    dst_geot = rasterHC.GetGeoTransform()
    dst_proj = osr.SpatialReference()
    dst_proj.ImportFromWkt(rasterHC.GetProjectionRef())
    dst_ds = driver.Create(
        outputRaster,
        rasterHC.RasterXSize,
        rasterHC.RasterYSize,
        3,
        dst_datatype)
    dst_ds.SetGeoTransform(dst_geot)
    dst_ds.SetProjection(dst_proj.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(1).WriteArray(mlwIDW)
    dst_ds.GetRasterBand(2).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(2).WriteArray(mslIDW)
    dst_ds.GetRasterBand(3).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(3).WriteArray(mhwIDW)
    #print ("  Output raster written successfully")
    #print ("")

    # --- PRINT TO SCREEN ---
    '''
    print ("")
    print ("--------------------------------------------------------")
    print ("--- COMPLETED HYDRO-MEM IN PYTHON - TIDAL DATUMS IDW ---")
    print ("--------------------------------------------------------")
    print ("")
    print ("")
    '''
