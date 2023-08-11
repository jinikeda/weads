# File: ponds.py
# Name: Karim Alizad
# Date: Nov 21, 2021
# Edits: Eric Swanson

"""
Example of how to run script from command line: python hyConn.py -i tbathy.img -o HyControl.img -r everdried.img
-i = input elevation file e.g. tbathy.img
-o = output raster file e.g. HyControl.img
-r = input raster file e.g. everdried.img
"""

#----------------------------------------------------------
# M O D U L E S                                   
#----------------------------------------------------------
#----------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from osgeo import osr
from scipy import ndimage

#----------------------------------------------------------
# F U N C T I O N    M A I N
#----------------------------------------------------------
#
#----------------------------------------------------------
def hyconn(inputElevation,inputRaster,outputRaster,DEBUG=False):    

    # --- GLOBAL PARAMETERS ---
    ndv=-99999.0
    #DEBUG = True

    # --- READ INPUTS ---
   # print ("")
   # print ("Reading raster")
    
    rasterED = gdal.Open(inputRaster) #everdried
    rasterTB = gdal.Open(inputElevation)

   # print ("  Raster of everdried (ED) read successfully")
   # print ("    Size (x):",rasterED.RasterXSize)
   # print ("    Size (y):",rasterED.RasterYSize)
   # print ("    Number of bands:",rasterED.RasterCount)

   # print ("  Raster of tbathy (TB) read successfully")
   # print ("    Size (x):",rasterTB.RasterXSize)
   # print ("    Size (y):",rasterTB.RasterYSize)
   # print ("    Number of bands:",rasterTB.RasterCount)


   # print ("")
   # print ("Reading Binary Rasters\n\n")
    EDband=rasterED.GetRasterBand(1)
    EDBN=EDband.ReadAsArray()
    TBband=rasterTB.GetRasterBand(1)
    TBBN=TBband.ReadAsArray()
   # original_everdried = EDBN

    if DEBUG:

        print('----- Original raster -----')
        plt.imshow(EDBN)
        plt.show()
        plt.savefig('terminal_everdried.png')

        plt.imshow(TBBN)
        plt.show()
        plt.savefig('terminal_tbathy.png')

    #compare everdried and tbathy
    # for i in range(0, EDBN.shape[0]):
    #     for j in range(0, EDBN.shape[1]):
    #         #clean up values to only have 1 and -99999
    #         if (EDBN[i][j] > 0.8) and (EDBN[i][j] < 1.2):
    #             EDBN[i][j] == 1.0   
                
    #         if (EDBN[i][j] < 0):
    #             EDBN[i][j] = -99999.0

    #         #all nodes considered wet by elevation set to 1
    #         #if (TBBN[i][j] < 0) and (EDBN[i][j] == -99999.0):
    #             #EDBN[i][j] = 1.0

    #         if (TBBN[i][j] == -99999.0):
    #             EDBN[i][j] = -99999.0

    ##########################################################
    
    Here, we can significantly speed up (Jin 2023 Aug 11)
    # Clean up values to only have 1 and -99999
    mask = (EDBN > 0.8) & (EDBN < 1.2)
    EDBN[mask] = 1.0

    EDBN[EDBN < 0] = -99999.0

    # Set all nodes considered wet by elevation to 1
    # EDBN[(TBBN < 0) & (EDBN == -99999.0)] = 1.0

    EDBN[TBBN == -99999.0] = -99999.0

    ##########################################################
    
    labeled, num_objects = ndimage.label(EDBN > 0.5)
   
    if DEBUG:

        print('----- Raster from everdried and tbathy -----')
        plt.imshow(EDBN)
        plt.axis('off')
        plt.show()

        print('----- Labeled raster -----')
        plt.imshow(labeled)
        plt.axis('off')
        plt.show()

    #print("\nCreating 3 raster bands\n")

    #sizes is array containing sizes of all the labeled objects
    sizes = ndimage.sum(EDBN, labeled, range(0, num_objects))
    #maximum is tuple where first element is index of largest object in the labeled image array
    maximum = np.where(sizes==sizes.max())[0]
    #create 1D mask arrays where each index that has max in corresponding labeled image is set to 1
    #max feature is 2D image array where  only the largest labels are highlighted
    max_index = np.zeros(num_objects + 1, np.float32)
    max_index[maximum] = 1
    max_feature = max_index[labeled]

    minimum = np.where(sizes==sizes.min())[0]
    min_index = np.zeros(num_objects + 1, np.float32)
    min_index[minimum] = 1
    min_feature = min_index[labeled]

    OCEAN = max_feature
    LAND = min_feature - (TBBN == -99999.0)
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

    out_array = np.zeros((EDBN.shape[0], EDBN.shape[1]), np.float32)
    #write value to 0 for OCEAN, 1 for LAND, or 2 for POND to one array
    for x in range(0, out_array.shape[0]):
        for y in range(0, out_array.shape[1]):
            if (OCEAN[x][y] == 1):
                out_array[x][y] = 0
            elif (LAND[x][y] == 1):
                out_array[x][y] = 1
            elif (POND[x][y] == 1):
                out_array[x][y] = 2
            else:
                out_array[x][y] = -99999.0

    if DEBUG:

        print('----- 3 band output raster -----')            
        plt.imshow(out_array)
        plt.axis('off')
        plt.savefig('HyControl.png')
        plt.show()

    # --- WRITE OUTPUTS ---
   # print ("")
    #print ("Writing output rasters")
    x_pixels = EDBN.shape[1]
    y_pixels = EDBN.shape[0]
    driver= gdal.GetDriverByName('HFA')
    dst_datatype=gdal.GDT_Float32
    dst_geo=rasterED.GetGeoTransform()
    dst_proj=osr.SpatialReference()
    dst_proj.ImportFromWkt(rasterED.GetProjectionRef())
    # dst_ds=driver.Create('HyCon_Binary.img',rasterED.RasterXSize, rasterED.RasterYSize, 3, dst_datatype)
    dst_ds=driver.Create(outputRaster,x_pixels, y_pixels, 1, dst_datatype)
    dst_ds.SetGeoTransform(dst_geo)
    dst_ds.SetProjection(dst_proj.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv)
    dst_ds.GetRasterBand(1).WriteArray(out_array)

    #close file, deallocate memory
    dst_ds = None 


  #  print ("Output rasters written successfully")
   # print ("")


### RUN MAIN ###
#if __name__ == "__main__":
 #   main(sys.argv[1:])
