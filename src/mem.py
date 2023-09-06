#!/usr/bin/python3
# File: hydroMEM.py
# Name: Peter Bacopoulos, Jin Ikeda
# Date: August 30, 2023
from osgeo import gdal
from osgeo import osr
import numpy as np
import math

# --- GLOBAL PARAMETERS ---
ndv = -99999.0  # No data value (ndv) using ADCIRC convention

#----------------------------------------------------------
# MEM: Marsh Equilibrium Model and Mangrove Development
#----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote 
# tidal datums outside of the hydraulically connect wet areas.
# result = function(inputRasterHyControl,inputRasterTopoBathy,\
#                   inputRasterTidalDatumsIDW,outputRaster)
#----------------------------------------------------------

#----------------------------------------------------------
# GLOSSARY
#----------------------------------------------------------
#
# A Accretion [cm]?
# B Biomass density [g m-2 yr -1]
# Bl Biomass density left side
# Br Biomass density right side
# tbA modified topo-bathy
# qstar q might vary with each species and location, moved to if conditions
# qstar2  q might vary with each species and location, moved to if conditions
# P Productivity
# marsh marsh classification

#----------------------------------------------------------
# DATASETS
#----------------------------------------------------------
#
# al=1000; bl=-3718; cl=1021;     # Biomass curve coefficients LHS-NorthInlet
# ar=1000; br=-3718; cr=1021;     # Biomass curve coefficients RHS-NorthInlet
#
# al=1.975; bl=-0.987; cl=1999;   # Biomass curve coefficients LHS-Apalachicola
# ar=3.265; br=-1.633; cr=1998;   # Biomass curve coefficients RHS-Apalachicola
#
# al=73.8; bl=-1.14; cl=1587.1;   # Biomass curve coefficients LHS-Weeks Bay
# ar=73.8; br=-1.14; cr=1587.1;   # Biomass curve coefficients RHS-Weeks Bay
#
# al=32; bl=-3.2; cl=1920;        # Biomass curve coefficients LHS-Grand Bay  Alizad et al. (2018)
# ar=6.61; br=-0.661; cr=1983;    # Biomass curve coefficients RHS-Grand Bay
#
# al=24.96; bl=-0.193; cl=592.7;  # Biomass curve coefficients LHS-Plum Island
# ar=24.96; br=-0.193; cr=592.7;  # Biomass curve coefficients RHS-Plum Island
#----------------------------------------------------------
# F U N C T I O N
#----------------------------------------------------------
def calculate_biomass_parabola(D, DNonNeg, al, bl, cl, ar, br, cr):

    # --- BIOMASS CALCULATIONS ---
    Bl = al * DNonNeg + bl * DNonNeg * DNonNeg + cl
    Bl[Bl < 0.0] = 0.0
    Bl[D >= 0.0] = 0.0

    Br = ar * DNonNeg + br * DNonNeg * DNonNeg + cr
    Br[Br < 0.0] = 0.0
    Br[D < 0.0] = 0.0

    B = Bl + Br
    B[D == 0] = ndv
    B[B == 0] = ndv

    Bmax = np.amax(B)
    PL = Bmax / 3
    PH = Bmax * 2 / 3

    return B, PL, PH

def create_raster(file,rasterHC, zarray):
    # Create the output raster dataset
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(file, rasterHC.RasterXSize, rasterHC.RasterYSize, rasterHC.RasterCount, gdal.GDT_Int32)
    out_ds.SetProjection(rasterHC.GetProjection())
    out_ds.SetGeoTransform(rasterHC.GetGeoTransform())
    dst_band = out_ds.GetRasterBand(1)
    dst_band.WriteArray(zarray)
    dst_band.SetNoDataValue(ndv)  # Exclude nodata value
    dst_band.ComputeStatistics(0)
    out_ds = None
    print('Made a raster file')
    return

def mem(inputRasterHyControl, inputRasterTopoBathy, \
        inputRasterTidalDatumsIDW,vegetationFile, outputRaster):

    '''
    print ("")
    print ("")
    print ("-------------------------------------------------")
    print ("--- LAUNCHING HYDRO-MEM IN PYTHON - HYDRO MEM ---")
    print ("-------------------------------------------------")
    print ("")
    '''

    # ----------------------------------------------------------
    # Read biomass calculation parameters
    # ----------------------------------------------------------

    if vegetationFile==None:
        print('\nA monotypic species with no vegetation mapping\n')
        #subscript:sub-optimal(left) and super-optimal(right) branches that met at the parabola apex

        al=32; bl=-3.2; cl=1920;        # Biomass curve coefficients LHS-Grand Bay  Alizad et al. (2018)
        ar=6.61; br=-0.661; cr=1983;    # Biomass curve coefficients RHS-Grand Bay

        # Accretion coefficients
        # Morris et al. (2016) Contributions of organic and inorganic matter to sedimentvolume and accretion in tidal wetlands at steady state
        BDi=1.99; BDo=0.085; Kr=0.2; m_const=0.0001;
        dt=5.0

    else:
        print('\nMulti species vegetation mapping\n')

        print ("Read a vegetation raster")
        rasterVG = gdal.Open(vegetationFile)
        transformVG = rasterVG.GetGeoTransform()
        prjVG = rasterVG.GetProjection()  # Read projection
        band = rasterVG.GetRasterBand(1);
        RV_VG = band.ReadAsArray();

        print(np.max(RV_VG), np.min(RV_VG),RV_VG.shape)

        unique_values, counts = np.unique(RV_VG, return_counts=True)
        print('Original unique_value is : ', unique_values)
        print('Their count is : ', counts)

        #### here is the example ####
        # 8 = salt marsh(regularly flooded) follow with tidal cycle
        # 9 = mangrove
        # 20 = irregularly flooded marsh
        # 40 = water_mask
        # 55 = land_mask
        # 128 = NAN for byte data
        #############################

        values_to_remove = [40, 55, 128]  # Values to remove from the array
        updated_unique_values = np.delete(unique_values, np.where(np.isin(unique_values, values_to_remove)))

        print('Interest unique_value is :', updated_unique_values)

        ################# Jin -> Pete box #######################################################################
        # Please work on this box

        # Create a mask for different vegetation types
        mask_salt_marsh = (RV_VG == 8)
        mask_mangrove = (RV_VG == 9)
        mask_irregular_marsh = (RV_VG == 20)
        dt = 5.0

        print ('The biomass parameters for salt marsh(regularly flooded)')
        ### later make function ####
        al=32; bl=-3.2; cl=1920;        # Biomass curve coefficients probably not parabola equation
        ar=6.61; br=-0.661; cr=1983;    # Biomass curve coefficients probably not parabola equation

        # Accretion coefficients
        q=2.8; # inorganic sediment load
        # organic and inorganic accumulation generated by decomposing vegetation
        BDi=1.99; BDo=0.085; Kr=0.2; m_const=0.0001;

        print('The biomass parameters for mangrove')
        print('Mangrove parameters goes here')

        print('The biomass parameters for salt marsh(regularly flooded)')
        print('Irregularly flooded marsh parameter  goes here')
        ### Jin to Pete ###: here, we need max inundation depth info as well

        ################# Jin -> Pete #######################################################################
    
    # ----------------------------------------------------------
    # Read raster files
    # ----------------------------------------------------------
    #print ("")
    #print ("Reading rasters")
      
    rasterHC=gdal.Open(inputRasterHyControl)
    gt=rasterHC.GetGeoTransform()
    X = rasterHC.RasterXSize
    Y = rasterHC.RasterYSize
    band = rasterHC.GetRasterBand(1);
    hc = band.ReadAsArray();
    prj = rasterHC.GetProjection()  # Read projection

    rasterTB=gdal.Open(inputRasterTopoBathy)
    rasterTDIDW=gdal.Open(inputRasterTidalDatumsIDW)

    '''
    print ("  Raster of HyControl (HC) read successfully")
    print ("    Size (x):",rasterHC.RasterXSize)
    print ("    Size (y):",rasterHC.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterHC.RasterCount)
    
    '''

    # ----------------------------------------------------------
    # Check projection and raster shape between provided vegetation and calculated raster files
    # ----------------------------------------------------------
    if vegetationFile==None:
        print("Projection:", prj)

    else:
        # Normalize both projection strings
        normalized_prj = " ".join(prj.strip().split()).split(',')[0] ##### need to modify this part later Aug 29
        normalized_prjVG = " ".join(prjVG.strip().split()).split(',')[0]

        # Compare the normalized projection strings
        if normalized_prj != normalized_prjVG or transformVG != gt:
            print("Vegetation raster does not match with current raster files! Need to modify the VG raster")
            print("Projection:", prj,'\n', prjVG)
            print("Projection:", normalized_prj,'\n', normalized_prjVG)
        else:
            print("Jin will modify this part later")
            print("Projections and Raster shape match!")
            print("Projection:", prj)

    # --- Make masks for fully wet and dried regions ---
    water_mask = ((-0.5 < hc) & (hc < 0.5)) | (1.5 < hc) # water regions (water and pond in hyconn.tiff)
    land_mask = (hc >= 0.5) & (hc <= 1.5) # land regions in hy conn.tiff

    # --- Read topo-bathymetry and water surface elevation ---
    band = rasterTB.GetRasterBand(1); tb=band.ReadAsArray(); #tb=-1.0*tb;
    band = rasterTDIDW.GetRasterBand(1); mlwIDW=band.ReadAsArray();
    # band=rasterTDIDW.GetRasterBand(2); mslIDW=band.ReadAsArray();
    band = rasterTDIDW.GetRasterBand(3); mhwIDW=band.ReadAsArray();

    # --- DEPTH CALCULATIONS ---")
    D = 100.0 * (mhwIDW - tb); D[tb < 0] = ndv; # Relative depth [cm]
    DNonNeg = D.copy()
    DNonNeg[D < 0.0] = 0.0;
    Dt = 100.0 * (mhwIDW - mlwIDW); # Tidal range [cm]
    Dt[Dt == 0] = 0.0001;

    # --- PERFORM HYDRO-MEM CALCULATIONS ---
    #print ("")
    #print ("Starting Hydro-MEM Calculations")

    if vegetationFile==None:

        P = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)  # Create an array of default values (ndv)
        marsh = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        # print ("Biomass Calculations")
        B, PL, PH = calculate_biomass_parabola(D, DNonNeg, al, bl, cl, ar, br, cr)

        # --- ACCRETION CALCULATIONS ---
        #print ("Accretion calculations")
        qstar = q*(DNonNeg/Dt); qstar2= qstar; qstar2[qstar>=1.0] = 1.0;
        w = hc.copy();
        w[w < 0.5] = 0.0;
        A= m_const*qstar2*DNonNeg/(BDi*2)+ Kr*B/(BDo*10000); A[D<=0.0]=0.0; A=A*w; A=A/100; A[A<=0]=0; tbA=tb+A*dt; tbA=-1.0*tbA;
        D[hc<0.5]=ndv; D[tb==np.nan]=ndv;
        B[hc<0.5]=ndv; B[tb==np.nan]=ndv;
        marsh[hc<0.5]=ndv; marsh[tb==np.nan]=ndv;
        A[hc<0.5]=ndv; A[tb==np.nan]=ndv;
        tbA[hc<0.5]=ndv; tbA[tb==np.nan]=ndv;
        A[D==0]=ndv; tbA[D==0]=ndv;
        #print ("  DONE!")

        # --- PRODUCTIVITY CALCULATIONS ---
        # Create masks for different conditions
        mask1 = (B > 1) & (B <= 1000)
        mask2 = (B > PL) & (B < 1800)
        mask3 = (B >= 1800)

        # Assign values based on conditions
        # background raster
        P[land_mask] = 55
        P[water_mask] = 40

        # Marsh productivity
        P[mask1] = 16 # low productivity
        P[mask2] = 23 # medium productivity
        P[mask3] = 32 # high productivity

        # --- MARSH TYPE CALCULATIONS ---
        #print ("High-Low Marsh Calculations")
        Dmax = -(al/(2*bl));
        Dzero1 = (-al + math.sqrt(((al*al)-(4*bl*cl))))/(2*bl);
        Dzero2 = (-ar - math.sqrt(((ar*ar)-(4*br*cr))))/(2*br);
        DRange = abs(Dzero2-Dzero1);
        DHigh = Dmax + DRange*0.1;
        DLow = Dmax - DRange*0.1;

        condition_1 = (Dzero1 < D) & (D <= DLow)
        condition_2 = (DLow < D) & (D < DHigh)
        condition_3 = (DHigh <= D) & (D < Dzero2)

        marsh[condition_1] = 30
        marsh[condition_2] = 20
        marsh[condition_3] = 10

        #print ("  DONE!")

    else:

        P = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)  # Create an array of default values (ndv)
        marsh = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        # print ("Biomass Calculations")
        B, PL, PH = calculate_biomass_parabola(D, DNonNeg, al, bl, cl, ar, br, cr)

        # --- ACCRETION CALCULATIONS ---
        #print ("Accretion calculations")
        qstar = q*(DNonNeg/Dt); qstar2= qstar; qstar2[qstar>=1.0] = 1.0;
        w = hc.copy();
        w[w < 0.5] = 0.0;
        A= m_const*qstar2*DNonNeg/(BDi*2)+ Kr*B/(BDo*10000); A[D<=0.0]=0.0; A=A*w; A=A/100; A[A<=0]=0; tbA=tb+A*dt; tbA=-1.0*tbA;
        D[hc<0.5]=ndv; D[tb==np.nan]=ndv;
        B[hc<0.5]=ndv; B[tb==np.nan]=ndv;
        marsh[hc<0.5]=ndv; marsh[tb==np.nan]=ndv;
        A[hc<0.5]=ndv; A[tb==np.nan]=ndv;
        tbA[hc<0.5]=ndv; tbA[tb==np.nan]=ndv;
        A[D==0]=ndv; tbA[D==0]=ndv;
        #print ("  DONE!")

        # --- PRODUCTIVITY CALCULATIONS ---
        # Create a mask for different conditions
        mask_regular_1 = mask_salt_marsh & (B > 1) & (B <= 1000)
        mask_regular_2 = mask_salt_marsh & (B > PL) & (B < 1800)
        mask_regular_3 = mask_salt_marsh & (B >= 1800)
        mask_irregular = mask_irregular_marsh # 'TBD' &
        mask_mangrove = mask_mangrove # 'TBD' &

        # Assign values based on conditions
        # background raster
        P[land_mask] = 55
        P[water_mask] = 40

        # Marsh productivity
        P[mask_regular_1] = 16 # low productivity
        P[mask_regular_2] = 23 # medium productivity
        P[mask_regular_3] = 32 # high productivity
        P[mask_irregular] = 100 # irregular_marsh think about later
        P[mask_mangrove] = 110 # mask_mangrove think about later

        # --- MARSH TYPE CALCULATIONS ---
        #print ("High-Low Marsh Calculations")
        Dmax = -(al/(2*bl));
        Dzero1 = (-al + math.sqrt(((al*al)-(4*bl*cl))))/(2*bl);
        Dzero2 = (-ar - math.sqrt(((ar*ar)-(4*br*cr))))/(2*br);
        DRange = abs(Dzero2-Dzero1);
        DHigh = Dmax + DRange*0.1;
        DLow = Dmax - DRange*0.1;

        condition_1 = (Dzero1 < D) & (D <= DLow)
        condition_2 = (DLow < D) & (D < DHigh)
        condition_3 = (DHigh <= D) & (D < Dzero2)
        condition_4 = Dzero2 <= D # irregular flood marsh TBD

        marsh[condition_1] = 30
        marsh[condition_2] = 20
        marsh[condition_3] = 10
        marsh[condition_4] = 100 # TBD

        # how about mangrove?
        #marsh[B < 1] = ndv

        #print (" DONE!")

    # --- WRITE OUTPUTS ---
    #print ("")
    #print ("Writing output raster")
    driver=gdal.GetDriverByName('HFA'); dst_datatype=gdal.GDT_Float32;
    dst_geot=rasterHC.GetGeoTransform(); dst_proj=osr.SpatialReference(); dst_proj.ImportFromWkt(rasterHC.GetProjectionRef());
    dst_ds=driver.Create(outputRaster,rasterHC.RasterXSize,rasterHC.RasterYSize,6,dst_datatype)
    dst_ds.SetGeoTransform(dst_geot); dst_ds.SetProjection(dst_proj.ExportToWkt());
    dst_ds.GetRasterBand(1).SetNoDataValue(ndv); dst_ds.GetRasterBand(1).WriteArray(D);
    dst_ds.GetRasterBand(2).SetNoDataValue(ndv); dst_ds.GetRasterBand(2).WriteArray(B);
    dst_ds.GetRasterBand(3).SetNoDataValue(ndv); dst_ds.GetRasterBand(3).WriteArray(A);
    dst_ds.GetRasterBand(4).SetNoDataValue(ndv); dst_ds.GetRasterBand(4).WriteArray(tbA);
    dst_ds.GetRasterBand(5).SetNoDataValue(ndv); dst_ds.GetRasterBand(5).WriteArray(marsh);
    dst_ds.GetRasterBand(6).SetNoDataValue(ndv); dst_ds.GetRasterBand(6).WriteArray(P);
    dst_ds = None

    ###### Renew vegetation map ####################
    VM = P.copy()
    VM[(VM == 16) | (VM == 23) | (VM == 32)] = 8 # 8 = salt marsh(regularly flooded) follow with tidal cycle
    VM[(VM == 110)] = 9 # 9 = mangrove
    VM[(VM == 100)] = 20. # 20 = irregularly flooded marsh

    create_raster('new_NWI.tif', rasterHC, VM)

    ###### Here Jin will create a raster file for WATTE

    create_raster('Productivity.tif', rasterHC, P)
    ############################################################################################################
    '''
    print ("Output raster written successfully")
    print ("")
    
    # --- PRINT TO SCREEN ---
    print ("")
    print ("-------------------------------------------------")
    print ("--- COMPLETED HYDRO-MEM IN PYTHON - HYDRO MEM ---")
    print ("-------------------------------------------------")
    print ("")
    print ("")
    '''

