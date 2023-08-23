#!/usr/bin/python3
# File: hydroMEM.py
# Name: Peter Bacopoulos, Jin Ikeda
# Date: August 22, 2023
from osgeo import gdal
from osgeo import osr
import numpy as np
import math

#----------------------------------------------------------
# F U N C T I O N	M E M	   
#----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote 
# tidal datums outside of the hydraulically connect wet areas.
# result = function(inputRasterHyControl,inputRasterTopoBathy,\
#                   inputRasterTidalDatumsIDW,outputRaster)
#----------------------------------------------------------
def mem(inputRasterHyControl, inputRasterTopoBathy,\
        inputRasterTidalDatumsIDW,outputRaster):

    '''
    print ("")
    print ("")
    print ("-------------------------------------------------")
    print ("--- LAUNCHING HYDRO-MEM IN PYTHON - HYDRO MEM ---")
    print ("-------------------------------------------------")
    print ("")
    '''
    
    # --- GLOBAL PARAMETERS ---
    ndv=-99999.0                      # No data value (ndv) using ADCIRC convention
    #al=1000; bl=-3718; cl=1021;     # Biomass curve coefficients LHS-NorthInlet
    #ar=1000; br=-3718; cr=1021;     # Biomass curve coefficients RHS-NorthInlet
    
    #al=1.975; bl=-0.987; cl=1999;   # Biomass curve coefficients LHS-Apalachicola
    #ar=3.265; br=-1.633; cr=1998;   # Biomass curve coefficients RHS-Apalachicola
    
    #al=73.8; bl=-1.14; cl=1587.1;   # Biomass curve coefficients LHS-Weeks Bay
    #ar=73.8; br=-1.14; cr=1587.1;   # Biomass curve coefficients RHS-Weeks Bay
    
    al=32; bl=-3.2; cl=1920;        # Biomass curve coefficients LHS-Grand Bay
    ar=6.61; br=-0.661; cr=1983;    # Biomass curve coefficients RHS-Grand Bay
    
    #al=24.96; bl=-0.193; cl=592.7;  # Biomass curve coefficients LHS-Plum Island
    #ar=24.96; br=-0.193; cr=592.7;  # Biomass curve coefficients RHS-Plum Island
        
    q=2.8; BDi=1.99; BDo=0.085; Kr=0.2; m_const=0.0001; dt=5.0     # Accretion coefficients
    
    # --- READ INPUTS ---
    #print ("")
    #print ("Reading rasters")
      
    rasterHC=gdal.Open(inputRasterHyControl)
    gt=rasterHC.GetGeoTransform()
    rasterTB=gdal.Open(inputRasterTopoBathy)
    gt=rasterTB.GetGeoTransform()
    rasterTDIDW=gdal.Open(inputRasterTidalDatumsIDW)
    gt=rasterTDIDW.GetGeoTransform()
    
    '''
    print ("  Raster of HyControl (HC) read successfully")
    print ("    Size (x):",rasterHC.RasterXSize)
    print ("    Size (y):",rasterHC.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterHC.RasterCount)
    
    print ("  Raster of topo-bathy (TB) read successfully")
    print ("    Size (x):",rasterTB.RasterXSize)
    print ("    Size (y):",rasterTB.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterTB.RasterCount)
        
    print ("  Raster of tidal datums (TD) IDW read successfully")
    print ("    Size (x):",rasterTDIDW.RasterXSize)
    print ("    Size (y):",rasterTDIDW.RasterYSize)
    pixelSizeX=gt[1]; pixelSizeY=-gt[5];
    print ("    Pixel size (x):",pixelSizeX)
    print ("    Pixel size (y):",pixelSizeY)
    print ("    Number of bands:",rasterTDIDW.RasterCount)
    '''
    
    # --- INITIALIZE ARRAYS ---
    band = rasterHC.GetRasterBand(1); hc=band.ReadAsArray();
    prj =   rasterHC.GetProjection()  # Read projection
    print("Projection:", prj)

    #print (np.max(hc), np.min(hc))
    # For Wave Attenuation Tool, we need water and land rasters
    water_mask = ((-0.5 < hc) & (hc < 0.5)) | (1.5 < hc) # water regions (water and pond in hyconn.tiff)
    land_mask = (hc >= 0.5) & (hc <= 1.5) # land regions in hy conn.tiff
    
    band=rasterTB.GetRasterBand(1); tb=band.ReadAsArray(); #tb=-1.0*tb;
    band=rasterTDIDW.GetRasterBand(1); mlwIDW=band.ReadAsArray();
    band=rasterTDIDW.GetRasterBand(2); mslIDW=band.ReadAsArray();
    band=rasterTDIDW.GetRasterBand(3); mhwIDW=band.ReadAsArray();
    B=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    marsh=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    A=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    tbA=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    D=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    Dt=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    DNonNeg=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    qstar=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    qstar2=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    Bl=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    Br=np.zeros((rasterHC.RasterYSize,rasterHC.RasterXSize),dtype=float)
    P=np.full((rasterHC.RasterYSize,rasterHC.RasterXSize), ndv) # Create an array of default values (ndv)
    
    # --- PERFORM HYDRO-MEM CALCULATIONS ---
    #print ("")
    #print ("Starting Hydro-MEM Calculations")
    
    # --- BIOMASS CALCULATIONS ---
    #print ("Biomass Calculations")
    D = 100.0*(mhwIDW-tb); D[tb<0]=ndv;
    Bl=al*D+bl*D*D+cl; Bl[Bl<0.0]=0.0; Bl[D>=0.0]=0.0;
    Br=ar*D+br*D*D+cr; Br[Br<0.0]=0.0; Br[D<0.0]=0.0;
    B=Bl+Br; B[D==0]=ndv; B[B==0]=ndv;
    Bmax=np.amax(B)
    PL=Bmax/3; PH=Bmax*2/3;
    
    #print ("  DONE!")
    
    # --- ACCRETION CALCULATIONS ---
    #print ("Accretion calculations")
    DNonNeg=D; DNonNeg[D<0.0]=0.0; Dt=100.0*(mhwIDW-mlwIDW);Dt[Dt==0]=0.0001;
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
    X=rasterHC.RasterXSize
    Y=rasterHC.RasterYSize

    # Create a mask for different conditions
    mask1 = (B > 1) & (B <= 1000)
    mask2 = (B > PL) & (B < 1800)
    mask3 = (B >= 1800)

    # Assign values based on conditions
    # background raster
    P[land_mask] = 55
    P[water_mask] = 40

    # Marsh productivity
    P[mask1] = 16 # low productivity
    P[mask2] = 23 # midum productivity
    P[mask3] = 32 # high productivity
          
    # --- MARSH TYPE CALCULATIONS ---
    #print ("High-Low Marsh Calculations")
    Dmax = -(al/(2*bl));
    Dzero1 = (-al + math.sqrt(((al*al)-(4*bl*cl))))/(2*bl);
    Dzero2 = (-ar - math.sqrt(((ar*ar)-(4*br*cr))))/(2*br);
    DRange = abs(Dzero2-Dzero1);
    DHigh = Dmax + DRange*0.1;
    DLow = Dmax - DRange*0.1;
    for i in range(0, Y):
        for j in range(0, X):
            if (Dzero1<D[i][j] and D[i][j]<=DLow):
                marsh[i][j] = 30
            elif (DLow<D[i][j] and D[i][j]<DHigh):
                marsh[i][j] = 20
            elif (DHigh<=D[i][j] and D[i][j]<Dzero2):
                marsh[i][j] = 10
            else:
                marsh[i][j] = ndv
    marsh[B<1]=ndv
    #print ("  DONE!")
    
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

    ###### Here Jin will create raster file for WATTE

    gtiff_driver = gdal.GetDriverByName('GTiff')  # Use GeoTIFF driver
    out_ds = gtiff_driver.Create('Productivity.tif',  # Create a output file
                                 rasterHC.RasterXSize, rasterHC.RasterYSize, rasterHC.RasterCount, gdal.GDT_Int32)
    out_ds.SetProjection(prj)
    out_ds.SetGeoTransform(dst_geot)

    dst_band = out_ds.GetRasterBand(1)
    dst_band.WriteArray(P)
    dst_band = out_ds.GetRasterBand(1).SetNoDataValue(int(ndv))  # Exclude nodata value
    dst_band = out_ds.GetRasterBand(1).ComputeStatistics(0)
    out_ds = None

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

