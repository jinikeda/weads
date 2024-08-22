#!/usr/bin/python3

# Crate a classification of bydro_type using tbathy and inundationtime. Since we didn't use raster (structure grid) data. So, probably we cannot distinguish water vs. pond/lake regions
# Modified: June 27, 2024 by Jin Ikeda

##########################################################################
# Development Notes:
""" This script is used to classify the land, ocean, (and pond/lake areas) based on inundationtime file (the value between 0,1 and -99999.0) to classify the subtidal,intertidal, land, and pond/lake areas.
The input file is the elevation file (e.g. ADCIRC_elevation.txt). The script is used in the WEAD project.
"""

""""
Example of how to run script from command line: python hyConn.py -i tbathy.img -o HyControl.img -r everdried.img
-i = input elevation file e.g. tbathy.txt
-o = output file e.g. HyControl.txt
-r = input file e.g. inundationtime.txt
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
def hydro_classify(inputElevation, inputInundationTFile, outputClass, DEBUG=False):
    # We may use csv file. But the code development will be conducted using
    # txt file for each step.

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0
    # DEBUG = True

    # --- READ INPUTS ---
    # print ("")
    # print ("Reading ADCIRC node")

    print(" Read an normalized inundation time (inpuT)")
    inunT = np.loadtxt(inputInundationTFile)  # domain_inundationtime.txt
    # inunTband = rasterinunT.GetRasterBand(1)
    inunTBN = xxxx  # need to check this part late (Jin June 27th, 2024)

    print(" Read an tbathy (TB)")
    TB = np.loadtxt(inputElevation)  # domain_elevation?.txt
    # TBband = rasterTB.GetRasterBand(1)
    # Here need to read column 4: Jin June 27th, 2024
    TBBN = xxxx

    # if DEBUG:
    #     print('----- Original raster -----')
    #     plt.imshow(inunTBN)
    #     plt.show()
    #     plt.savefig('terminal_everdried.png')
    #
    #     plt.imshow(TBBN)
    #     plt.show()
    #     plt.savefig('terminal_tbathy.png')

    # compare inundationtime and tbathy and Clean up values to only have 1 and
    # -99999
    accuracy = 1.0 * 10**-6

    """ inunTBN = -99999.0: nodata,
                0: land, 1: intertidal, 2: subtidal """  # , 3: pond/lake"""

    mask_land = (0 - accuracy < inunTBN) & (inunTBN < 0 + accuracy)
    mask_intertidal = (accuracy < inunTBN) & (inunTBN < 1 - accuracy)
    mask_water = (1 - accuracy < inunTBN) & (inunTBN < 1 + accuracy)
    # out of domain
    # mask_outdomain = (TBBN == ndv) This part will not be needed. However, we will keep at this stage

    inunTBN[mask_land] = int(0)  # fully dried (land) region
    inunTBN[mask_intertidal] = int(1)  # intertidal region
    # temporary set to water region and will separate to ocean (subtidal zone)
    # and pond/lake
    inunTBN[mask_water] = int(2)
    # set nodata value to -99999.0 in the domain of TBBN
    inunTBN[inunTBN <= 0 - accuracy] = int(ndv)

##########################################################################
# This block can be removed on point based WEADS
##########################################################################

    # labeled, num_objects = ndimage.label(mask_water) # labeled connected region and provide the number of regions
    #
    # if DEBUG:
    #     print('----- Raster from inundationtime and tbathy -----')
    #     plt.imshow(inunTBN)
    #     plt.axis('off')
    #     plt.show()
    #
    #     print('----- Labeled raster -----')
    #     plt.imshow(labeled)
    #     plt.axis('off')
    #     plt.show()
    #
    # print("\nCreating 3 raster bands\n")
    #
    # # sizes is array containing sizes of all the labeled objects
    # sizes = ndimage.sum(mask_water, labeled, range(num_objects+1)) # calculate the size of each regions
    # print('size max and min',sizes.max(),sizes.min(), ',number of areas', num_objects) #max size should be ocean and min size should be land
    #
    # # maximum is tuple where first element is index of largest object in the labeled image array # This approach may not work where ocean is not a largest size
    # maximum = np.where(sizes == sizes.max())[0]
    # # create 1D mask arrays where each index that has max in corresponding labeled image is set to 1
    # # max feature is 2D image array where  only the largest labels are highlighted
    # max_index = np.zeros(num_objects + 1, np.float32)
    # max_index[maximum] = 1
    # max_feature = max_index[labeled]
    #
    # ############################################################################################################
    # # This command eliminates small regions and replace as land
    #
    # mask_size = sizes < 10 # user defined
    # remove_pixel = mask_size[labeled]
    # remove_pixel.shape
    # labeled[remove_pixel]=np.argmin(sizes)
    #
    # labels = np.unique(labeled)
    # labeled1 = np.searchsorted(labels, labeled)
    #
    # ############################################################################################################
    #
    # minimum = np.where(sizes==sizes.min())[0]
    # min_index = np.zeros(num_objects + 1, np.float32)
    # min_index[minimum] = 1
    # min_feature = min_index[labeled]
    # # print('maxf: ', max_feature,'minf: ', min_feature)
    # OCEAN = max_feature
    # LAND = min_feature - (TBBN == -99999.0) # binary 1: True, 0: False
    # POND = 1 - np.add(OCEAN, LAND) - (TBBN == -99999.0) # Use this classification to separate from water region
    # if DEBUG:
    #     print('----- Ocean raster -----')
    #     plt.imshow(OCEAN)
    #     plt.axis('off')
    #     plt.savefig('OCEAN.png')
    #     plt.show()
    #
    #     print('----- Land raster -----')
    #     plt.imshow(LAND)
    #     plt.axis('off')
    #     plt.savefig('LAND.png')
    #     plt.show()
    #
    #     print('----- Pond/lake raster -----')
    #     plt.imshow(POND)
    #     plt.axis('off')
    #     plt.savefig('POND.png')
    #     plt.show()

##########################################################################

    # #Keep 1D array rather than 2D array.
    # out_array = np.zeros((inunTBN.shape[0], inunTBN.shape[1]), np.float32) # create 2D array with 0 some locations are not completely filled in the raster classification
    # out_array[mask_intertidal] = 1
    # out_array[OCEAN == 1] = 2 # subtidal zone
    # out_array[POND == 1] = 3 # pond/lake
    # out_array[mask_outdomain] = ndv # set nodata value to -99999.0 in the
    # domain of TBBN

    # if DEBUG:
    #
    #     print('----- 4 band output raster -----')
    #     plt.imshow(out_array)
    #     plt.axis('off')
    #     plt.savefig('Hydro_class.png')
    #     plt.show()
    #
    # # --- WRITE OUTPUTS ---
    # # print ("")
    # #print ("Writing output rasters")
    # x_pixels = inunTBN.shape[1]
    # y_pixels = inunTBN.shape[0]
    # driver= gdal.GetDriverByName('HFA')
    # dst_datatype=gdal.GDT_Float32
    # dst_geo=rasterinunT.GetGeoTransform()
    # dst_proj=osr.SpatialReference()
    # dst_proj.ImportFromWkt(rasterinunT.GetProjectionRef())
    # # dst_ds=driver.Create('HyCon_Binary.img',rasterinunT.RasterXSize, rasterinunT.RasterYSize, 3, dst_datatype)
    # dst_ds=driver.Create(outputRaster,x_pixels, y_pixels, 1, dst_datatype)
    # dst_ds.SetGeoTransform(dst_geo)
    # dst_ds.SetProjection(dst_proj.ExportToWkt())
    # dst_ds.GetRasterBand(1).SetNoDataValue(ndv)
    # dst_ds.GetRasterBand(1).WriteArray(out_array)
    #
    # #close file, deallocate memory
    # dst_ds = None

   # print ("Output rasters written successfully")
   # print ("")


### RUN MAIN ###
# if __name__ == "__main__":
 #   main(sys.argv[1:])
