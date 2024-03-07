#!/usr/bin/python3
# File: hydroMEM.py
# Name: Peter Bacopoulos, Jin Ikeda
# Date: March 7, 2024
from osgeo import gdal
from osgeo import osr
import numpy as np
import math

#----------------------------------------------------------
# MEM: Marsh Equilibrium Model and Mangrove Development
#----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums outside the hydraulically connect wet areas.
# result = function(inputRasterHyControl,inputRasterTopoBathy,\
#                   inputRasterTidalDatumsIDW,outputRaster)
#----------------------------------------------------------

#----------------------------------------------------------
# GLOSSARY OF VARIABLES # Jin -> Pete please update this list
#----------------------------------------------------------
#
# A: Accretion [cm]?
# B: Biomass density [g m-2 yr -1]
# Bl: Biomass density left side
# Br: Biomass density right side
# tbA: modified topo-bathy
# qstar: q might vary with each species and location, moved to if conditions
# qstar2:  q might vary with each species and location, moved to if conditions
# P: Productivity
# marsh: marsh classification
# Dmin, Dmax, Dopt
# Bmax

#----------------------------------------------------------
# DATASETS Jin -> Pete add States
#----------------------------------------------------------
#
# al=1000; bl=-3718; cl=1021;     # Biomass curve coefficients LHS-NorthInlet
# ar=1000; br=-3718; cr=1021;     # Biomass curve coefficients RHS-NorthInlet
#
# al=1.975; bl=-0.987; cl=1999;   # Biomass curve coefficients LHS-Apalachicola, FL
# ar=3.265; br=-1.633; cr=1998;   # Biomass curve coefficients RHS-Apalachicola, FL
#
# al=73.8; bl=-1.14; cl=1587.1;   # Biomass curve coefficients LHS-Weeks Bay, LA?
# ar=73.8; br=-1.14; cr=1587.1;   # Biomass curve coefficients RHS-Weeks Bay, LA
#
# al=32; bl=-3.2; cl=1920;        # Biomass curve coefficients LHS-Grand Bay  Alizad et al. (2018)
# ar=6.61; br=-0.661; cr=1983;    # Biomass curve coefficients RHS-Grand Bay
#
# al=24.96; bl=-0.193; cl=592.7;  # Biomass curve coefficients LHS-Plum Island
# ar=24.96; br=-0.193; cr=592.7;  # Biomass curve coefficients RHS-Plum Island

biomass_coefficients = {
    "NorthInlet": {"al": 1000, "bl": -3718, "cl": 1021, "ar": 1000, "br": -3718, "cr": 1021},
    "Apalachicola": {"al": 1.975, "bl": -0.987, "cl": 1999, "ar": 3.265, "br": -1.633, "cr": 1998},
    "WeeksBay": {"al": 73.8, "bl": -1.14, "cl": 1587.1, "ar": 73.8, "br": -1.14, "cr": 1587.1},
    "GrandBay": {"al": 32, "bl": -3.2, "cl": 1920, "ar": 6.61, "br": -0.661, "cr": 1983},
    "PlumIsland": {"al": 24.96, "bl": -0.193, "cl": 592.7, "ar": 24.96, "br": -0.193, "cr": 592.7}}

# Need to add Optimum Elevation

########################################################################################################################
print ("Input parameters for MEM")
########################################################################################################################
# --- GLOBAL PARAMETERS ---
ndv = -99999.0  # No data value (ndv) using ADCIRC convention
qmax = 2.8 # Maximum capture coefficient (-)
qmin = 1.0 # Minimum capture coefficient (-)
SSC = 25.0 # Suspended sediment concentration
FF = 353.0 # Flooding frequency (1/year)
k1 = 0.085  # Organic inputs # # Bulk density of organic matter (g/cm3)
k2 = 1.99  # Mineral inputs # Bulk density of inorganic matter (g/cm3)
dt = 5.0  # Time step (yr)

# --- LOCAL ACCRETION PARAMETERS ---

print('The accretion parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')
# Salt marsh (8)
Dmin1 = 2.0;
Dmax1 = 46.0;
Bmax1 = 2400.0;
Dopt1 = 22.0;

# Mangrove (9)
Dmin2 = 24.0
Dmax2 = 66.0
Bmax2 = 7800.0
Dopt2 = 45.0

# Irregularly flooded marsh (20)
Dmin3 = -21.0;
Dmax3 = 35.0;
Bmax3 = 1200.0;
Dopt3 = 7.0;  # 3: Irregularly flooded marsh (NWI = 20)


# ORGANIC: Morris MEM 2016
kr1 = 0.1 # Refractory fraction (-)
kr2 = 0.1 # Refractory fraction (-)
kr3 = 0.1 # Refractory fraction (-)
RRS1 = 2.0 # BG-to-shoot bio ratio (-)
RRS2 = 1.8 # BG-to-shoot bio ratio (-)
RRS3 = 1.5 # BG-to-shoot bio ratio (-)
BTR1 = 0.5 # BG turnover rate (1/year)
BTR2 = 0.25 # BG turnover rate (1/year) for mature mangrove
########################################################################################################################
### Jin to Pete ###: Check the BTR here and what is BG? is not used in the code ########################################
BTRmat = 0.22  # BG turnover rate (1/year) for mature mangroves
BTRjuv = 0.67  # BG turnover rate (1/year) for juvenile mangroves
########################################################################################################################
BTR3 = 0.5
Tmat = 30  # Time for pioneer mangroves to fully mature (yr)

#----------------------------------------------------------
# F U N C T I O N
#----------------------------------------------------------
def calculate_coefficients(Dmin, Dmax, Dopt, Bmax):
    a = -((-Dmin*Bmax-Dmax*Bmax)/((Dmin-Dopt)*(Dopt-Dmax))) # coefficient of DNonNeg
    b = -(Bmax/((Dmin-Dopt)*(Dopt-Dmax))) # coefficient of DNonNeg^2
    c = -Dmin*Bmax*Dmax/((Dmin-Dopt)*(Dopt-Dmax)) # constant term
    return a, b, c


# this is a problem
def calculate_biomass_parabola(D, DNonNeg,ndv, al, bl, cl, ar, br, cr):

    # --- BIOMASS CALCULATIONS ---
    Bl = al * DNonNeg + bl * DNonNeg * DNonNeg + cl
    Bl[Bl < 0.0] = ndv
    #Bl[Bl ==np.nan] = ndv
    #Bl[D >= 0.0] = 0.0
    print(Bl.min(), Bl.max())

    Br = ar * DNonNeg + br * DNonNeg * DNonNeg + cr
    Br[Br < 0.0] = ndv
    #Br[Br ==np.nan] = ndv

    B = Bl + Br
    B[D == 0] = 0.0
    B[B == 0] = 0.0
    B[B < 0] = ndv

    Bmax = B.max()
    PL = Bmax / 3
    PH = Bmax * 2 / 3

    return B, PL, PH

def calculate_biomass_parabola_beta(DNonNeg,Dopt,ndv, al, bl, cl, ar, br, cr):

    # --- BIOMASS CALCULATIONS ---
    Bl = al * DNonNeg + bl * DNonNeg * DNonNeg + cl
    Bl[Bl < 0.0] = ndv
    Bl[DNonNeg >= Dopt] = 0.0
    print(Bl.min(), Bl.max())

    Br = ar * DNonNeg + br * DNonNeg * DNonNeg + cr
    Br[Br < 0.0] = ndv
    Br[DNonNeg < Dopt] = 0.0
    #Br[Br ==np.nan] = ndv
    print(Br.min(), Br.max())

    B = Bl + Br
    # B[DNonNeg == 0] = 0.0
    B[B == 0] = 0.0
    B[B < 0] = ndv

    Bmax = B.max()
    print('check Bmax', Bmax)
    PL = Bmax / 3
    PH = Bmax * 2 / 3

    return B, PL, PH



def calculate_vertical_accretion(qmin, qmax,dt, B, Bmax, kr, RRS, BTR, SSC, FF, D, k1, k2):

    # --- ACCRETION CALCULATIONS ---
    q=qmin+(qmax-qmin)*B/Bmax
    q[B<0.0]=qmin
    Vorg=kr*RRS*BTR*B/(100.0**2) # organic matter
    Vmin=0.5*q*SSC*FF*D/(1000.0**2) # inorganic matter where is Jin -> Pete FIT?
    Vorg[B<0.0]=0.0
    Vmin[B<0.0]=0.0
    A=(Vorg/k1)+(Vmin/k2) # accretion rate
    tbA=A*dt/100.0 # accretion thickness

    return tbA, A

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


def set_ndv(array, conditions):
    for condition in conditions:
        array[condition] = ndv
    return array

def mem(inputRasterHyControl, inputRasterTopoBathy, inputRasterTidalDatumsIDW,vegetationFile, outputRaster):

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
    global dt,ndv, q, B, PL, PH, A, tbA, marsh, qstar, qstar2, w, DNonNeg, D, Dmin, Dmax, Dopt, Bmax, k1, k2, kr, RRS, BTR, SSC, FF

    if vegetationFile==None:
        print('\nA monotypic species with no vegetation mapping\n')
        #subscript:sub-optimal(left) and super-optimal(right) branches that met at the parabola apex

        # Select reference sites to get Biomass coefficients
        coefficients = biomass_coefficients["NorthInlet"]
        al = coefficients["al"]
        bl = coefficients["bl"]
        cl = coefficients["cl"]
        ar = coefficients["ar"]
        br = coefficients["br"]
        cr = coefficients["cr"]

        # Accretion coefficients
        # Morris et al. (2016) Contributions of organic and inorganic matter to sedimentvolume and accretion in tidal wetlands at steady state
        q = 2.8;  # inorganic sediment load
        BDi=1.99; BDo=0.085; Kr=0.2; m_const=0.0001;

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
        mask_irregular = (RV_VG == 20)
        dt = 5.0

        print ('The biomass parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')

        a1, b1, c1 = calculate_coefficients(Dmin1, Dmax1, Dopt1, Bmax1) # 1: Salt marsh (NWI = 8)
        a2, b2, c2 = calculate_coefficients(Dmin2, Dmax2, Dopt2, Bmax2)  # 2: Mangrove (NWI = 9) (mangrove (mature))
        a3, b3, c3 = calculate_coefficients(Dmin3, Dmax3, Dopt3, Bmax3) # 3: Irregularly flooded marsh (NWI = 20)
        print('check: abc', a1, b1, c1)

########################################################################################################################
        # Jin to Pete: how to incorporate juvenile mangrove here

        # Dmin2 = 36
        # Dmax2 = 55
        # Bmax2 = 500
        # Dopt2 = 45
        # Dmin2 = DminL + dD + dDmin
        # Dmax2 = DmaxR + dD + dDmax
        # Dopt2 = ((DoptL + DoptR) / 2) + dD
########################################################################################################################

        ### Jin to Pete ###: here, we need max inundation depth info as well
        ### Pete to Jin (Why, for WATTE?) Jin comments to max inundation on Feb 15 2024
        ### The code doesn't include max inundation depth info, but the file may specify the range and inundation depth of irregularly flooded marshes—also for WATTE.
        # The rasterize file was created through hydromem.py. But Jin will determine how to incorporate it into mem.py later (no action is needed at this moment).
    
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
    print(DNonNeg.min(),DNonNeg.max())
    Dt = 100.0 * (mhwIDW - mlwIDW); # Tidal range [cm]
    Dt[Dt == 0] = 0.0001;

    # --- PERFORM HYDRO-MEM CALCULATIONS ---
    #print ("")
    #print ("Starting Hydro-MEM Calculations")

    if vegetationFile==None:

        P = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)  # Create an array of default values (ndv)
        marsh = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        print ("Biomass Calculations")
        B, PL, PH = calculate_biomass_parabola(D, DNonNeg,ndv, al, bl, cl, ar, br, cr)

        # --- ACCRETION CALCULATIONS ---
        print ("Accretion calculations")
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
        mangrove = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)
        irregular = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        # Salt marsh (8)
        B1, PL1, PH1 = calculate_biomass_parabola_beta(DNonNeg, Dopt1, ndv, a1, b1, c1, a1, b1, c1)
        print('Salt marsh ','PH = ', PH1, ', PL =' ,PL1, 'Bmin',B1.min())

        # Mangrove (9)
        B2, PL2, PH2 = calculate_biomass_parabola_beta(DNonNeg, Dopt2, ndv, a2, b2, c2, a2, b2, c2)
        print('Mangrove ','PH = ', PH2, ', PL =',PL2,'Bmin',B2.min())

        # Irregularly marsh (20)
        B3, PL3, PH3 = calculate_biomass_parabola_beta(DNonNeg, Dopt3, ndv, a3, b3, c3, a3, b3, c3)
        print('Irregularly marsh ','PH = ', PH3, ', PL =',PL3,'Bmin',B3.min())

        # --- ACCRETION CALCULATIONS ---
        w=hc.copy(); w[w<0.5]=0.0;


        tbA, A =  calculate_vertical_accretion(qmin, qmax, dt, B1, Bmax1, kr1, RRS1, BTR1, SSC, FF, D, k1, k2)
        tbA1=tbA*w; A1=A*w;

        tbA, A = calculate_vertical_accretion(qmin, qmax, dt, B2, Bmax2, kr2, RRS2, BTR2, SSC, FF, D, k1, k2)
        tbA2=tbA*w; A2=A*w;

        tbA, A = calculate_vertical_accretion(qmin, qmax, dt, B3, Bmax3, kr3, RRS3, BTR3, SSC, FF, D, k1, k2)
        tbA3=tbA*w; A3=A*w;

        # # set_non values
        # ndv_conditions = [hc < 0.5,tb == np.nan] #, tb == np.nan # hc < 0.5 is ocean region
        # #ndv_conditions2 = [hc < 0.5, D ==  0,tb==np.nan]  # ,
        # D = set_ndv(D, ndv_conditions)
        # B1 = set_ndv(B1, ndv_conditions)
        # B2 = set_ndv(B2, ndv_conditions)
        # B3 = set_ndv(B3, ndv_conditions)
        # marsh = set_ndv(marsh, ndv_conditions)
        # mangrove = set_ndv(mangrove, ndv_conditions)
        # irregular = set_ndv(irregular, ndv_conditions)
        # A1 = set_ndv(A1, ndv_conditions)
        # A2 = set_ndv(A2, ndv_conditions)
        # A3 = set_ndv(A3, ndv_conditions)
        # tbA1 = set_ndv(tbA1, ndv_conditions)
        # tbA2 = set_ndv(tbA2, ndv_conditions)
        # tbA3 = set_ndv(tbA3, ndv_conditions)

        # --- PRODUCTIVITY CALCULATIONS ---
        # Create a mask for different conditions
        mask_regular_1 = mask_salt_marsh & (B1 <= PL1) & (B1 >= 0) #& (B1 <= PL1)
        mask_regular_2 = mask_salt_marsh & (B1 > PL1) & (B1 < PH1)
        mask_regular_3 = mask_salt_marsh & (B1 >= PH1)
        mask_mangrove = mask_mangrove #& (B2 > 0) & (B2 <= PH2)
        mask_irregular = mask_irregular  #& (B3 > 0) & (B3 <= PH3)
        # mask_mangrove = mask_mangrove #& (B2 > 0) & (B2 <= PL2)
        # mask_mangrove = mask_mangrove #& (B2 > PL2) & (B2 < PH2)
        # mask_mangrove = mask_mangrove #& (B2 >= PH2)
        # mask_irregular = mask_irregular #& (B3 > 0) & (B3 <= PL3)
        # mask_irregular = mask_irregular #& (B3 > PL3) & (B3 < PH3)
        # mask_irregular = mask_irregular #& (B3 >= PH3)

        # Assign values based on conditions
        # background raster
        P[land_mask] = 55
        P[water_mask] = 40

        # Marsh productivity
        P[mask_regular_1] = 16 # low productivity
        P[mask_regular_2] = 23 # medium productivity
        P[mask_regular_3] = 32 # high productivity
        P[mask_mangrove] = 109 # mangrove
        P[mask_irregular] = 120 # irregularly flooded marsh

        # --- MARSH TYPE CALCULATIONS ---
########################################################################################################################
        # Jin to Pete: We need to modify this part later
        print ("High-Low Marsh Calculations")
        ## here how to determine the range of the Dmax? mannuarryly? or from the raster file?
        ### Need to discuss with Pete!

        al = a1
        ar = a1
        bl = b1
        br = b1
        cl = c1
        cr = c1

        Dmax = -(al/(2*bl));
        Dzero1 = (-al + math.sqrt(((al*al)-(4*bl*cl))))/(2*bl);
        Dzero2 = (-ar - math.sqrt(((ar*ar)-(4*br*cr))))/(2*br);
        DRange = abs(Dzero2-Dzero1);
        DHigh = Dmax + DRange*0.1;
        DLow = Dmax - DRange*0.1;

        condition_1 = (Dzero1 < D) & (D <= DLow)
        condition_2 = (DLow < D) & (D < DHigh)
        condition_3 = (DHigh <= D) & (D < Dzero2)
        condition_4 = Dzero2 <= D

        marsh[condition_1] = 30
        marsh[condition_2] = 20
        marsh[condition_3] = 10

        irregular[condition_4] = 120

        # mangrove juvenile and mature
        # mangrove [xxx]

########################################################################################################################
        ### need to discuss with Pete about this part
        # Mask ndv values in the arrays
        B1_masked = np.ma.masked_equal(B1, ndv)
        B2_masked = np.ma.masked_equal(B2, ndv)
        B3_masked = np.ma.masked_equal(B3, ndv)

        A1_masked = np.ma.masked_equal(A1, ndv)
        A2_masked = np.ma.masked_equal(A2, ndv)
        A3_masked = np.ma.masked_equal(A3, ndv)

        tbA1_masked = np.ma.masked_equal(tbA1, ndv)
        tbA2_masked = np.ma.masked_equal(tbA2, ndv)
        tbA3_masked = np.ma.masked_equal(tbA3, ndv)

        # Perform the addition
        B = B1_masked + B2_masked + B3_masked # Biomass density [g m-2 yr -1]
        A = A1_masked + A2_masked + A3_masked # accretion rate
        tbA = tbA1_masked + tbA2_masked + tbA3_masked # accretion thickness
    ####################################################################################################################
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
    VM[(VM == 109)] = 9 # 9 = mangrove
    VM[(VM == 120)] = 20 # 20 = irregularly flooded marsh

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

