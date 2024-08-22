#!/usr/bin/python3
# File: originally ponds.py by Karim Alizad (Aug 17, 2023)
# Modified: March 8, 2024 by Jin Ikeda

##########################################################################
# Development Notes:
""" This script is used to classify the land, ocean, and pond/lake areas based on everdried.tif. However, Jin Ikeda changes
the script to use the input inundationtime raster file (0,1 and -99999.0) to classify the subtidal,intertidal, land, and pond/lake areas.
The input raster file is the elevation file (e.g. tbathy.img). The script is used in the WEAD project.
"""
""""
Example of how to run script from command line: python hyConn.py -i tbathy.img -o HyControl.img -r everdried.img
-i = input elevation file e.g. tbathy.img
-o = output raster file e.g. HyControl.img
-r = input raster file e.g. everdried.img
"""

# ----------------------------------------------------------
# M O D U L E S
# ----------------------------------------------------------
# ----------------------------------------------------------

# ----------------------------------------------------------
# F U N C T I O N    M A I N
# ----------------------------------------------------------
#
# ----------------------------------------------------------




import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from osgeo import osr
from scipy import ndimage
def hydro_classify(inputElevation, inputRaster, outputRaster, DEBUG=False):
    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0
    # DEBUG = True

    # --- READ INPUTS ---
    # print ("")
    # print ("Reading raster")

    rasterinunT = gdal.Open(inputRaster)  # inundationtime.tif
    rasterTB = gdal.Open(inputElevation)

    print(" Raster of normalized inundation time (inpuT) read successfully")
    print("    Size (x):", rasterinunT.RasterXSize)
    print("    Size (y):", rasterinunT.RasterYSize)
    print("    Number of bands:", rasterinunT.RasterCount)

    print(" Raster of tbathy (TB) read successfully")
    print("    Size (x):", rasterTB.RasterXSize)
    print("    Size (y):", rasterTB.RasterYSize)
    print("    Number of bands:", rasterTB.RasterCount)

    print("")
    print("Reading Binary Rasters\n\n")
    inunTband = rasterinunT.GetRasterBand(1)
    inunTBN = inunTband.ReadAsArray()
    TBband = rasterTB.GetRasterBand(1)
    TBBN = TBband.ReadAsArray()
    # original_everdried = inunTBN

    if DEBUG:
        print('----- Original raster -----')
        plt.imshow(inunTBN)
        plt.show()
        plt.savefig('terminal_everdried.png')

        plt.imshow(TBBN)
        plt.show()
        plt.savefig('terminal_tbathy.png')

    # compare everdried and tbathy and Clean up values to only have 1 and
    # -99999
    accuracy = 1.0 * 10**-6

    """ inunTBN = -99999.0: nodata,
                0: land, 1: intertidal, 2: subtidal, 3: pond/lake"""

    mask_land = (0 - accuracy < inunTBN) & (inunTBN < 0 + accuracy)
    mask_intertidal = (accuracy < inunTBN) & (inunTBN < 1 - accuracy)
    mask_water = (1 - accuracy < inunTBN) & (inunTBN < 1 + accuracy)
    mask_outdomain = (TBBN == ndv)  # out of domain

    inunTBN[mask_land] = 0.0  # fully dried (land) region
    inunTBN[mask_intertidal] = 1.0  # intertidal region
    # temporary set to water region and will separate to ocean (subtidal zone)
    # and pond/lake
    inunTBN[mask_water] = 2.0
    # set nodata value to -99999.0 in the domain of TBBN
    inunTBN[inunTBN <= 0 - accuracy] = ndv

    # labeled connected region and provide the number of regions
    labeled, num_objects = ndimage.label(mask_water)

    if DEBUG:
        print('----- Raster from inundationtime and tbathy -----')
        plt.imshow(inunTBN)
        plt.axis('off')
        plt.show()

        print('----- Labeled raster -----')
        plt.imshow(labeled)
        plt.axis('off')
        plt.show()

    # print("\nCreating 3 raster bands\n")

    # sizes is array containing sizes of all the labeled objects
    # calculate the size of each regions
    sizes = ndimage.sum(mask_water, labeled, range(num_objects + 1))
    print('size max and min', sizes.max(), sizes.min(), ',number of areas',
          num_objects)  # max size should be ocean and min size should be land

    # maximum is tuple where first element is index of largest object in the
    # labeled image array # This approach may not work where ocean is not a
    # largest size
    maximum = np.where(sizes == sizes.max())[0]
    # create 1D mask arrays where each index that has max in corresponding labeled image is set to 1
    # max feature is 2D image array where  only the largest labels are
    # highlighted
    max_index = np.zeros(num_objects + 1, np.float32)
    max_index[maximum] = 1
    max_feature = max_index[labeled]

    ##########################################################################
    # This command eliminates small regions and replace as land

    mask_size = sizes < 10  # user defined
    remove_pixel = mask_size[labeled]
    remove_pixel.shape
    labeled[remove_pixel] = np.argmin(sizes)

    labels = np.unique(labeled)
    labeled1 = np.searchsorted(labels, labeled)

    ##########################################################################

    minimum = np.where(sizes == sizes.min())[0]
    min_index = np.zeros(num_objects + 1, np.float32)
    min_index[minimum] = 1
    min_feature = min_index[labeled]
    # print('maxf: ', max_feature,'minf: ', min_feature)
    OCEAN = max_feature
    LAND = min_feature - (TBBN == -99999.0)  # binary 1: True, 0: False
    # Use this classification to separate from water region
    POND = 1 - np.add(OCEAN, LAND) - (TBBN == -99999.0)
    if DEBUG:
        print('----- Ocean raster -----')
        plt.imshow(OCEAN)
        plt.axis('off')
        plt.savefig('OCEAN.png')
        plt.show()

        print('----- Land raster -----')
        plt.imshow(LAND)
        plt.axis('off')
        plt.savefig('LAND.png')
        plt.show()

        print('----- Pond/lake raster -----')
        plt.imshow(POND)
        plt.axis('off')
        plt.savefig('POND.png')
        plt.show()

    # create 2D array with 0 some locations are not completely filled in the
    # raster classification
    out_array = np.zeros((inunTBN.shape[0], inunTBN.shape[1]), np.float32)
    #out_array[mask_land] = 0
    out_array[mask_intertidal] = 1
    out_array[OCEAN == 1] = 2  # subtidal zone
    out_array[POND == 1] = 3  # pond/lake
    # set nodata value to -99999.0 in the domain of TBBN
    out_array[mask_outdomain] = ndv

    if DEBUG:

        print('----- 4 band output raster -----')
        plt.imshow(out_array)
        plt.axis('off')
        plt.savefig('Hydro_class.png')
        plt.show()

    # --- WRITE OUTPUTS ---
    # print ("")
    #print ("Writing output rasters")
    x_pixels = inunTBN.shape[1]
    y_pixels = inunTBN.shape[0]
    driver = gdal.GetDriverByName('HFA')
    dst_datatype = gdal.GDT_Float32
    dst_geo = rasterinunT.GetGeoTransform()
    dst_proj = osr.SpatialReference()
    dst_proj.ImportFromWkt(rasterinunT.GetProjectionRef())
    # dst_ds=driver.Create('HyCon_Binary.img',rasterinunT.RasterXSize, rasterinunT.RasterYSize, 3, dst_datatype)
    dst_ds = driver.Create(outputRaster, x_pixels, y_pixels, 1, dst_datatype)
    dst_ds.SetGeoTransform(dst_geo)
    dst_ds.SetProjection(dst_proj.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(1).WriteArray(out_array)

    # close file, deallocate memory
    dst_ds = None

   # print ("Output rasters written successfully")
   # print ("")


### RUN MAIN ###
# if __name__ == "__main__":
 #   main(sys.argv[1:])
