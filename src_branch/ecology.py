#!/usr/bin/python3
# File: ecology.py
# Name: Peter Bacopoulos, Jin Ikeda
# last modify: July 4, 2024

#--- Load modules ---
import numpy as np
import math
import time
import pandas as pd

#--- Initialize code ---
start=time.time()
print("\n")
print("LAUNCH: Launching script!\n")

#----------------------------------------------------------
# MEM: Marsh Equilibrium Model and Mangrove Development
#----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums outside fully wet areas.
### Need to modify the following function

# result = function(inputRasterHyControl,inputRasterTopoBathy,\
#                   inputRasterTidalDatumsIDW,outputRaster)
#----------------------------------------------------------

# Probably make a dictionary python file for the following
########################################################################################################################
#----------------------------------------------------------
# GLOSSARY OF VARIABLES #
#----------------------------------------------------------
#
# A: Accretion [cm]?
# B: Biomass density [g m-2 yr -1]
# Bl: Biomass density left side
# Br: Biomass density right side
# tb_update: modified topo-bathy
# qstar: q might vary with each species and location, moved to if conditions
# qstar2:  q might vary with each species and location, moved to if conditions
# P: Productivity
# marsh: marsh classification
# Dmin, Dmax, Dopt
# Bmax

#----------------------------------------------------------
# DATASETS Jin -> Pete Please add States and Optimum March 7, 2024
#----------------------------------------------------------
#
# al=1000; bl=-3718; cl=1021;     # Biomass curve coefficients LHS-NorthInlet, where?
# ar=1000; br=-3718; cr=1021;     # Biomass curve coefficients RHS-NorthInlet, where?
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

interest_reference = "Texas_Coastal_Bend"

# Use a dictionary to store the coefficients
biomass_coefficients = {
    "NorthInlet": {"al": 1000, "bl": -3718, "cl": 1021, "ar": 1000, "br": -3718, "cr": 1021, "Dopt": 22}, # can be include, "Bmax": 2400},
    "Apalachicola": {"al": 1.975, "bl": -0.987, "cl": 1999, "ar": 3.265, "br": -1.633, "cr": 1998},
    "WeeksBay": {"al": 73.8, "bl": -1.14, "cl": 1587.1, "ar": 73.8, "br": -1.14, "cr": 1587.1},
    "GrandBay": {"al": 32, "bl": -3.2, "cl": 1920, "ar": 6.61, "br": -0.661, "cr": 1983, "Dopt": 5.0},
    "PlumIsland": {"al": 24.96, "bl": -0.193, "cl": 592.7, "ar": 24.96, "br": -0.193, "cr": 592.7},
    "Texas_Coastal_Bend": {"al": 240.0, "bl": -5.0, "cl": -460.0, "ar": 240.0, "br": -5.0, "cr": -460.0, "Dopt": 22.0},
    "Texas_Coastal_Bend_mangrove": {"al": 1600, "bl": -17.7, "cl": -28016.0, "ar": 1600, "br": -17.7, "cr": -28016.0, "Dopt": 45.0}}

# Need to add Optimum Elevation at least to run calculate_biomass_parabola
# until here
########################################################################################################################


########################################################################################################################
print ("Input parameters for MEM")
########################################################################################################################
# --- GLOBAL PARAMETERS ---
ndv = -99999.0  # No data value (ndv) using ADCIRC conversion
ndv_byte = 128
qmax = 2.8    # Maximum capture coefficient (-)
qmin = 1.0    # Minimum capture coefficient (-)
SSC = 25.0    # Suspended sediment concentration
FF = 353.0    # Flooding frequency (1/year)
BDo = 0.085   # Organic inputs # # Bulk density of organic matter (g/cm3)
BDi = 1.99      # Mineral inputs # Bulk density of inorganic matter (g/cm3)
Kr = 0.2      # Refractory fraction (-)

# --- LOCAL ACCRETION PARAMETERS ---
# when vegetationFile == True: modify the following part
print('The accretion parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')

########################################################################################################################
vegetation_parameters = {
    "SaltMarsh": {
        "Dmin": 2.0,
        "Dmax": 46.0,
        "Bmax": 2400.0,
        "Dopt": 22.0,
        "Kr": 0.1,
        "RRS": 2.0,
        "BTR": 0.5
    },
    "Mangrove": {
        "Dmin": 24.0,
        "Dmax": 66.0,
        "Bmax": 7800.0,
        "Dopt": 45.0,
        "Kr": 0.1,
        "RRS": 1.8,
        "BTR": 0.25
    },
    "IrregularMarsh": {
        "Dmin": -70.0,  # -21.0 -> 0.5*std below the mean of the marsh
        "Dmax": -10.0,  # 35.0 -> 0.5*std above the mean of the marsh
        "Bmax": 1200.0,
        "Dopt": -40.0,  # 7.0 # modified based on tb and vegetation mapping
        "Kr": 0.1,
        "RRS": 1.5,
        "BTR": 0.5
    }
}

########################################################################################################################
### Jin to Pete ###: Check the BTR here and what is BG? is not used in the code March 7, 24 ############################
BTRmat = 0.22  # BG turnover rate (1/year) for mature mangroves
BTRjuv = 0.67  # BG turnover rate (1/year) for juvenile mangroves
########################################################################################################################
Tmat = 30  # Time for pioneer mangroves to fully mature (yr)


# ----------------------------------------------------------
# F U N C T I O N
# ----------------------------------------------------------
# def read_band_value(raster_file,
#                     band_num):  # Future modify the function read transform and projection information March 11, 2024 (Jin)
#     raster = gdal.Open(raster_file)
#     band = raster.GetRasterBand(band_num).ReadAsArray()
#     return band


def calculate_coefficients(Dmin, Dmax, Dopt, Bmax):
    a = -((-Dmin * Bmax - Dmax * Bmax) / ((Dmin - Dopt) * (Dopt - Dmax)))  # coefficient of DNonNeg
    b = -(Bmax / ((Dmin - Dopt) * (Dopt - Dmax)))  # coefficient of DNonNeg^2
    c = -Dmin * Bmax * Dmax / ((Dmin - Dopt) * (Dopt - Dmax))  # constant term
    return a, b, c


def calculate_biomass_parabola(DNonNeg, Dopt, mask, al, bl, cl, ar, br, cr):
    global ndv
    # Jin -> Pete please think about using D or DNonNeg

    # --- BIOMASS CALCULATIONS ---
    # Create a mask for the condition
    mask_left_parabora = (DNonNeg <= Dopt) & mask
    mask_right_parabora = (Dopt < DNonNeg) & mask

    # Perform the calculations only where the condition is true
    Bl = np.where(mask_left_parabora, al * DNonNeg + bl * DNonNeg * DNonNeg + cl, ndv)
    Br = np.where(mask_right_parabora, ar * DNonNeg + br * DNonNeg * DNonNeg + cr, ndv)

    print('BL:', Bl.min(), Bl.max())
    print('BR:', Br.min(), Br.max())

    # Set negative values to negative small value
    Bl[Bl < 0.0] = -0.001
    Br[Br < 0.0] = -0.001

    B = Bl + Br
    B[B < 0] = ndv

    Bmax = B.max()
    print('B max and min', Bmax, B.min())
    PL = Bmax / 3
    PH = Bmax * 2 / 3

    return B, PL, PH


def calculate_vertical_accretion(qmin, qmax, dt, B, Bmax, kr, RRS, BTR, SSC, FF, D, BDo, BDi, tb):
    # --- ACCRETION CALCULATIONS ---
    q = qmin + (qmax - qmin) * B / Bmax
    q[B < 0.0] = qmin
    Vorg = kr * RRS * BTR * B / (100.0 ** 2)  # organic matter
    Vmin = 0.5 * q * SSC * FF * D / (1000.0 ** 2)  # inorganic matter where is Jin -> Pete FIT?
    Vorg[B <= 0.0] = 0.0
    Vmin[
        B <= 0.0] = 0.0  # If we also evaluate the accretion rate for mineral matter, we will change the code using mask such as (mhwIDW != ndv)
    A = np.where(tb != ndv, ((Vorg / BDo) + (Vmin / BDi)) / 100.0, ndv)  # accretion rate per year [m]
    tb_update = np.where(tb != ndv, tb + (A * dt), ndv)  # accretion thickness [m]

    return tb_update, A


# def create_raster(file, rasterHC, zarray, dtype, no_data_value, stats_flag=False):
#     # Create the output raster dataset
#     gtiff_driver = gdal.GetDriverByName('GTiff')
#     out_ds = gtiff_driver.Create(file, rasterHC.RasterXSize, rasterHC.RasterYSize, rasterHC.RasterCount,
#                                  dtype)  # dtype is e.g. gdal.GDT_Int32 and gdal.GDT_Float32
#     out_ds.SetProjection(rasterHC.GetProjection())
#     out_ds.SetGeoTransform(rasterHC.GetGeoTransform())
#     dst_band = out_ds.GetRasterBand(1)
#     dst_band.WriteArray(zarray)
#     dst_band.SetNoDataValue(no_data_value)  # Exclude nodata value
#     stats = dst_band.ComputeStatistics(0)
#     min_val, max_val, mean_val, std_dev_val = stats
#     if stats_flag:
#         print(
#             f'Made a raster file. Statistics:\n\tMinimum: {min_val}, Maximum: {max_val}, Mean: {mean_val}, Standard Deviation: {std_dev_val}')
#     else:
#         print('Made a raster file')
#     out_ds = None
#     return


#def mem(inputRasterHyControl, inputRasterTopoBathy, inputRasterTidalDatumsIDW, vegetationFile, outputRaster, deltaT=5):
# define from hydroMEM.py


### manual inputs
deltaT = 5  # Time step (yr)
vegetationFile = None  # Vegetation file (Early development)

dt = deltaT  # Time step (yr)

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
# Turn on the global when make a function
#global ndv, ndv_byte, qmax, qmin, SSC, FF, BDo, BDi, Kr, biomass_coefficients, vegetation_parameters, interest_reference

if vegetationFile == None:
    print('\nA monotypic species with no vegetation mapping\n')
    # subscript:sub-optimal(left) and super-optimal(right) branches that met at the parabola apex

    # Select reference sites to get Biomass coefficients
    coefficients = biomass_coefficients[interest_reference]
    al = coefficients["al"]
    bl = coefficients["bl"]
    cl = coefficients["cl"]
    ar = coefficients["ar"]
    br = coefficients["br"]
    cr = coefficients["cr"]
    Dopt = coefficients["Dopt"]

    # Accretion coefficients
    # Morris et al. (2016) Contributions of organic and inorganic matter to sedimentvolume and accretion in tidal wetlands at steady state

    m_const = 0.0001;

else:
    print('\nMulti species vegetation mapping\n')

    # print("Read a vegetation raster")
    # rasterVG = gdal.Open(vegetationFile)
    # transformVG = rasterVG.GetGeoTransform()
    # prjVG = rasterVG.GetProjection()  # Read projection
    # band = rasterVG.GetRasterBand(1);
    # RV_VG = band.ReadAsArray();
    #
    # print(np.max(RV_VG), np.min(RV_VG), RV_VG.shape)
    #
    # unique_values, counts = np.unique(RV_VG, return_counts=True)
    # print('Original unique_value is : ', unique_values)
    # print('Their count is : ', counts)
    #
    # #### here is the example ####
    # # 8 = salt marsh(regularly flooded) follow with tidal cycle
    # # 9 = mangrove
    # # 20 = irregularly flooded marsh
    # # 40 = water_mask
    # # 55 = land_mask
    # # 128 = NAN for byte data
    # #############################
    #
    # values_to_remove = [40, 55, 128]  # Values to remove from the array
    # updated_unique_values = np.delete(unique_values, np.where(np.isin(unique_values, values_to_remove)))
    #
    # print('Interest unique_value is :', updated_unique_values)
    #
    # ################# Jin -> Pete box #######################################################################
    # # Please work on this box
    #
    # # Create a mask for different vegetation types
    # mask_salt_marsh = (RV_VG == 8)
    # mask_mangrove = (RV_VG == 9)
    # mask_irregular = (RV_VG == 20) | (RV_VG == 8)  # try
    #
    # print('The biomass parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')
    #
    # # Salt marsh
    # vegetation_type = "SaltMarsh"
    # params = vegetation_parameters[vegetation_type]
    # Dmin_marsh = params["Dmin"]
    # Dmax_marsh = params["Dmax"]
    # Dopt_marsh = params["Dopt"]
    # Bmax_marsh = params["Bmax"]
    # kr_marsh = params["Kr"]
    # RRS_marsh = params["RRS"]
    # BTR_masrh = params["BTR"]
    #
    # a1, b1, c1 = calculate_coefficients(Dmin_marsh, Dmax_marsh, Dopt_marsh, Bmax_marsh)  # 1: Salt marsh (NWI = 8)
    # print('Check: SaltMarsh abc', a1, b1, c1)
    #
    # # Mangrove
    # vegetation_type = "Mangrove"
    # params = vegetation_parameters[vegetation_type]
    # Dmin_matmang = params["Dmin"]
    # Dmax_matmang = params["Dmax"]
    # Dopt_matmang = params["Dopt"]
    # Bmax_matmang = params["Bmax"]
    # kr_matmang = params["Kr"]
    # RRS_matmang = params["RRS"]
    # BTR_matmang = params["BTR"]
    #
    # a2, b2, c2 = calculate_coefficients(Dmin_matmang, Dmax_matmang, Dopt_matmang,
    #                                     Bmax_matmang)  # 2: Mangrove (NWI = 9) (mangrove (mature))
    # print('Check: Mangrove abc', a2, b2, c2)
    #
    # # Irregularly flooded marsh
    # vegetation_type = "IrregularMarsh"
    # params = vegetation_parameters[vegetation_type]
    # Dmin_fmarsh = params["Dmin"]
    # Dmax_fmarsh = params["Dmax"]
    # Dopt_fmarsh = params["Dopt"]
    # Bmax_fmarsh = params["Bmax"]
    # kr_fmarsh = params["Kr"]
    # RRS_fmarsh = params["RRS"]
    # BTR_fmarsh = params["BTR"]
    #
    # a3, b3, c3 = calculate_coefficients(Dmin_fmarsh, Dmax_fmarsh, Dopt_fmarsh,
    #                                     Bmax_fmarsh)  # 20: Irregularly flooded marsh (NWI = 20)

########################################################################################################################
# Jin to Pete: how to incorporate juvenile mangrove here

# Dmin_matmang = 36
# Dmax_matmang = 55
# Bmax_matmang = 500
# Dopt_matmang = 45
# Dmin_matmang = DminL + dD + dDmin
# Dmax_matmang = DmaxR + dD + dDmax
# Dopt_matmang = ((DoptL + DoptR) / 2) + dD
########################################################################################################################

# ----------------------------------------------------------
# Read point files
# ----------------------------------------------------------

# ndv = -99999.0  # No data value (ndv) using ADCIRC convention

inEPSG = 4269  # GCS_North_American_1983
outEPSG = 26914
#################################################################################
PRJ = "EPSG:" + str(outEPSG)

# --- READ INPUTS ---
# Read the CSV file
df = pd.read_csv("tidal_prj.csv")
print("  Read HyControl (HC) and tidal datums (TD) successfully")

print(df.shape, df.columns, df.dtypes)



    # print ("Reading rasters")

    # rasterHC = gdal.Open(inputRasterHyControl)
    # gt = rasterHC.GetGeoTransform()
    # X = rasterHC.RasterXSize
    # Y = rasterHC.RasterYSize
    # prj = rasterHC.GetProjection()  # Read projection
    # hc = read_band_value(inputRasterHyControl, 1)

hc = df['HydroClass'] # 'HydroClass' 0: land, 1: intertidal, 2: subtidal(water)

# ----------------------------------------------------------
# Check projection and raster shape between provided vegetation and calculated raster files
# ----------------------------------------------------------
########################################################################################################################
# gdal prj will be different from the one that we have. Need to modify this part July 4, 2024
########################################################################################################################
if vegetationFile == None:
    prj = PRJ
    print("Projection:", prj)

else:
    pass
#     # Normalize both projection strings
#     normalized_prj = " ".join(prj.strip().split()).split(',')[0]  ##### need to modify this part later Aug 29
#     normalized_prjVG = " ".join(prjVG.strip().split()).split(',')[0]
#
#     # Compare the normalized projection strings
#     if normalized_prj != normalized_prjVG or transformVG != gt:
#         print("Vegetation raster does not match with current raster files! Need to modify the VG raster")
#         print("Projection:", prj, '\n', prjVG)
#         print("Projection:", normalized_prj, '\n', normalized_prjVG)
#     else:
#         print("Jin will modify this part later")
#         print("Projections and Raster shape match!")
#         print("Projection:", prj)

# --- Read topo-bathymetry and water surface elevation ---
tb = df['z']  # tb=-1.0*tb;
print('tb: min&max', np.min(tb), np.max(tb))

# --- Read interpolated tidal datums ---

mlwIDW = df['MLW_IDW']
mslIDW = df['MSL_IDW']
mhwIDW = df['MHW_IDW']

# --- DEPTH CALCULATIONS ---")
# Make masks for fully wet and dried regions ---
land_mask = (-0.5 < hc) & (hc <= 0.5)  # land region (hc = 0.0)
water_mask = (1.5 < hc) & (hc < 2.5)  # fully wet region (hc = 2.0)
intertidal_mask = (0.5 < hc) & (hc < 1.5)  # intertidal region (hc = 1.0)
submergence_mask = water_mask | intertidal_mask | (
            2.5 < hc)  # submergence region (hc = 1.0, 2.0, 3.0 including lake and pond) # these values are not existed current point-based code
above_subtidal_zone = (mhwIDW != ndv)  # above subtidal zone ()
np.savetxt('above_subtidal_zone.csv', above_subtidal_zone, delimiter=',')

# Perform the calculation only where the condition is true
D = np.where(above_subtidal_zone, 100.0 * (mhwIDW - tb),
             -ndv)  # [cm] Relative depth at nan should be positive value in this case

D[D > -ndv / 2] = ndv;  # Relative depth [cm] re-modify to negative value for ndv

DNonNeg = D.copy()
print(DNonNeg.min(), DNonNeg.max())

Dt = 100.0 * (mhwIDW - mlwIDW);  # Tidal range [cm]
print('Tidal range: min and max [cm]\t',Dt.min(), Dt.max())
Dt[Dt == 0] = 0.0001;

# --- PERFORM HYDRO-MEM CALCULATIONS ---
print("Starting Hydro-MEM Calculations")

if vegetationFile == None:

    P = np.full((df.shape[0], 1), ndv,
                dtype=float)  # Create an array of default values (ndv)
    marsh = np.full((df.shape[0], 1), ndv, dtype=float)

    # --- BIOMASS CALCULATIONS ---
    print("Biomass Calculations")
    B, PL, PH = calculate_biomass_parabola(DNonNeg, Dopt, above_subtidal_zone, al, bl, cl, ar, br, cr)
    np.where((B <= 0.0), 0.0, B)
    print('Salt marsh ', 'PH = ', PH, ', PL = ', PL, 'B max = ', B.max())

    # --- ACCRETION CALCULATIONS ---
    print("Accretion calculations")

    qstar2 = np.where(qmax * (DNonNeg / Dt) >= 1.0, 1.0, qmax * (DNonNeg / Dt))

    A = np.where((tb != ndv) & (above_subtidal_zone) & (B != ndv),
                 (m_const * qstar2 * SSC * FF * DNonNeg / (BDi * 2) + Kr * B / (BDo * 10000))/100, 0)  # accretion rate per year [m]

    B[tb == ndv] = ndv
    A[tb == ndv] = ndv

    tb_update = np.where(above_subtidal_zone, tb + (A * dt), tb)  # accretion thickness [m]
    tb_update[tb == ndv] = ndv

    # --- PRODUCTIVITY CALCULATIONS ---
    # Create masks for different conditions
    mask1 = (B > 0) & (B <= 1000)
    mask2 = (B > PL) & (B < 1800)
    mask3 = (B >= 1800)

    # Assign values based on conditions
    # background raster
    P[land_mask] = 55
    P[submergence_mask] = 40

    # Marsh productivity
    P[mask1] = 16  # low productivity
    P[mask2] = 23  # medium productivity
    P[mask3] = 32  # high productivity

    # --- MARSH TYPE CALCULATIONS ---
    # print ("High-Low Marsh Calculations")
    Dmax = -(al / (2 * bl));
    Dzero1 = (-al + math.sqrt(((al * al) - (4 * bl * cl)))) / (2 * bl);
    Dzero2 = (-ar - math.sqrt(((ar * ar) - (4 * br * cr)))) / (2 * br);
    DRange = abs(Dzero2 - Dzero1);
    DHigh = Dmax + DRange * 0.1;
    DLow = Dmax - DRange * 0.1;

    condition_1 = (Dzero1 < D) & (D <= DLow)
    condition_2 = (DLow < D) & (D < DHigh)
    condition_3 = (DHigh <= D) & (D < Dzero2)

    marsh[condition_1] = 30
    marsh[condition_2] = 20
    marsh[condition_3] = 10

else:

    pass
#
#     ##### Here check the elevation of mask irregularly flooded marsh ###############################################
#     # This part is for tuning the parameters for irregularly flooded marsh (not mandatory
#     NWI_original = read_band_value('/work/jinikeda/ETC/TCB/pyHydroMEM_dev/Resampled_raster_domain.tif',
#                                    1)  # Before dilation
#     mask_irregular_NWI = (NWI_original == 20)
#     irregular_elevation = tb[mask_irregular_NWI & above_subtidal_zone]
#     irregular_NWI_depth = D[mask_irregular_NWI & above_subtidal_zone]
#     print('Irregular_depth [cm], min, max, mean, std', irregular_NWI_depth.min(), irregular_NWI_depth.max(),
#           irregular_NWI_depth.mean(), irregular_NWI_depth.std())
#     ################################################################################################################
#
#     P = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv,
#                 dtype=float)  # Create an array of default values (ndv)
#     marsh = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)
#     mangrove = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)
#     irregular = np.full((rasterHC.RasterYSize, rasterHC.RasterXSize), ndv, dtype=float)
#
#     # --- BIOMASS CALCULATIONS ---
#     # Salt marsh (8)
#     B1, PL1, PH1 = calculate_biomass_parabola(DNonNeg, Dopt_marsh, above_subtidal_zone, a1, b1, c1, a1, b1, c1)
#     print('Salt marsh ', 'PH = ', PH1, ', PL = ', PL1)
#
#     # Mangrove (9)
#     B2, PL2, PH2 = calculate_biomass_parabola(DNonNeg, Dopt_matmang, above_subtidal_zone, a2, b2, c2, a2, b2, c2)
#     # print('Mangrove ','PH = ', PH2, ', PL = ',PL2)
#
#     # Irregularly marsh (20)
#     B3, PL3, PH3 = calculate_biomass_parabola(DNonNeg, Dopt_fmarsh, above_subtidal_zone, a3, b3, c3, a3, b3, c3)
#     # print('Irregularly marsh ','PH = ', PH3, ', PL =',PL3)
#
#     # --- ACCRETION CALCULATIONS ---
#     _, A1 = calculate_vertical_accretion(qmin, qmax, dt, B1, Bmax_marsh, kr_marsh, RRS_marsh, BTR_masrh, SSC, FF, D,
#                                          BDo, BDi, tb)
#     _, A2 = calculate_vertical_accretion(qmin, qmax, dt, B2, Bmax_matmang, kr_matmang, RRS_matmang, BTR_matmang,
#                                          SSC, FF, D, BDo, BDi, tb)
#     _, A3 = calculate_vertical_accretion(qmin, qmax, dt, B3, Bmax_fmarsh, kr_fmarsh, RRS_fmarsh, BTR_fmarsh, SSC,
#                                          FF, D, BDo, BDi, tb)
#
#     # --- PRODUCTIVITY CALCULATIONS ---
#     # Create a mask for different conditions
#     mask_regular_1 = mask_salt_marsh & (0 < B1) & (B1 <= PL1)
#     mask_regular_2 = mask_salt_marsh & (PL1 < B1) & (B1 < PH1)
#     mask_regular_3 = mask_salt_marsh & (PH1 <= B1)
#     mask_mangrove = mask_mangrove & (
#                 0 < B2)  # & (B2 <= PH2) # This part keep for future revision (0 < B2) & (B2 <= PH2)
#     mask_irregular = mask_irregular & (0 < B3)  # & (B3 <= PH3)
#     # mask_mangrove_mat = mask_mangrove & future revision
#     # mask_mangrove_juv = mask_mangrove & future revision
#     # mask_mangrove = mask_mangrove #& (B2 >= PH2)
#
#     # Assign values based on conditions
#     # From background raster first. here ndv value will not overwrite the background
#     P[land_mask] = 55
#     P[submergence_mask] = 40
#
#     # Marsh productivity (Caution for the oder of the mask)
#     P[mask_irregular] = 120  # irregularly flooded marsh
#     P[mask_regular_1] = 16  # low productivity
#     P[mask_regular_2] = 23  # medium productivity
#     P[mask_regular_3] = 32  # high productivity
#     P[mask_mangrove] = 109  # mangrove
#
#     # --- MARSH TYPE CALCULATIONS ---
#     ########################################################################################################################
#     # Jin to Pete: We need to check the following part
#     print("High-Low Marsh Calculations")
#
#     al = a1
#     ar = a1
#     bl = b1
#     br = b1
#     cl = c1
#     cr = c1
#
#     Dmax = -(al / (2 * bl));
#     Dzero1 = (-al + math.sqrt(((al * al) - (4 * bl * cl)))) / (2 * bl);
#     Dzero2 = (-ar - math.sqrt(((ar * ar) - (4 * br * cr)))) / (2 * br);
#     DRange = abs(Dzero2 - Dzero1);
#     DHigh = Dmax + DRange * 0.1;
#     DLow = Dmax - DRange * 0.1;
#
#     condition_1 = (Dzero1 < D) & (D <= DLow)
#     condition_2 = (DLow < D) & (D < DHigh)
#     condition_3 = (DHigh <= D) & (D < Dzero2)
#     condition_4 = Dzero2 <= D
#
#     marsh[condition_1] = 30
#     marsh[condition_2] = 20
#     marsh[condition_3] = 10
#
#     irregular[condition_4] = 120
#
#     # mangrove juvenile and mature
#     # mangrove [xxx]
#
#     # Mask ndv values in the arrays
#     B1_masked = np.where((mask_salt_marsh | mask_mangrove), B1,
#                          ndv)  # mask_marsh may compete with mangrove due to the dilution order, so evaluate both masks here.
#     B2_masked = np.where(mask_mangrove, B2, ndv)
#     B3_masked = np.where((mask_salt_marsh | mask_mangrove | mask_irregular), B3, ndv)
#
#     A1_masked = np.where((mask_salt_marsh | mask_mangrove), A1, 0)  # background of entire grid is set at 0 [m]
#     A2_masked = np.where(mask_mangrove, A2, 0)
#     A3_masked = np.where((mask_salt_marsh | mask_mangrove | mask_irregular), A3, 0)
#
#     print('Amax', A1_masked.max(), A2_masked.max(), A3_masked.max())
#
#     # Merge the arrays (used maximum value for each cell)
#     B = np.maximum.reduce([B1_masked, B2_masked, B3_masked])  # Maximum biomass density on each grid [g m-2 yr -1]
#     A = np.maximum.reduce([A1_masked, A2_masked, A3_masked])  # Maximum accretion rate per year on each grid [m]
#     A[tb == ndv] = ndv  # Overwrite the outside domain as  the ndv values
#
#     # Update topo-bathy data
#     tb_update = np.where((mask_salt_marsh | mask_mangrove | mask_irregular), tb + (A * dt), tb)
#     tb_update[tb == ndv] = ndv  # update topo-bathy data due to some holes are existed inside the domain
#
####################################################################################################################
# --- WRITE OUTPUTS ---
# print ("")
# print ("Writing point base output raster")
df['D'] = D.flatten() # mhwIDW - tb [cm]
df['B'] = B.flatten()
df['A'] = A.flatten()
df['tb_update'] = tb_update.flatten()
df['marsh'] = marsh.flatten()
df['P'] = P.flatten()

df.to_csv('mem.csv', index=False)



# driver = gdal.GetDriverByName('HFA');
# dst_datatype = gdal.GDT_Float32;
# dst_geot = rasterHC.GetGeoTransform();
# dst_proj = osr.SpatialReference();
# dst_proj.ImportFromWkt(rasterHC.GetProjectionRef());
# dst_ds = driver.Create(outputRaster, rasterHC.RasterXSize, rasterHC.RasterYSize, 6, dst_datatype)
# dst_ds.SetGeoTransform(dst_geot);
# dst_ds.SetProjection(dst_proj.ExportToWkt());
# dst_ds.GetRasterBand(1).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(1).WriteArray(D);
# dst_ds.GetRasterBand(2).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(2).WriteArray(B);
# dst_ds.GetRasterBand(3).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(3).WriteArray(A);
# dst_ds.GetRasterBand(4).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(4).WriteArray(tb_update);
# dst_ds.GetRasterBand(5).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(5).WriteArray(marsh);
# dst_ds.GetRasterBand(6).SetNoDataValue(ndv);
# dst_ds.GetRasterBand(6).WriteArray(P);
# dst_ds = None
#
###### Renew vegetation map ####################
if vegetationFile is not None:
    VM = P.copy()
    VM[(VM == 16) | (VM == 23) | (VM == 32)] = 8  # 8 = salt marsh(regularly flooded) follow with tidal cycle
    VM[(VM == 109)] = 9  # 9 = mangrove
    VM[(VM == 120)] = 20  # 20 = irregularly flooded marsh
    # VM[(VM == 55)] = ndv_byte # 55 = land_mask
    # VM[(VM == 40)] = ndv_byte # 40 = water_mask

    mask = (VM != 8) & (VM != 9) & (VM != 20)
    VM[mask] = ndv_byte
#
# ####################################################################################################################
# print("\n----------------------------------------------------------------------")
# print("Create WATTE input files")
# print("----------------------------------------------------------------------\n")
# ####################################################################################################################
#
# # Output
# Watte_bio_level = P.copy()  # Create an array of default values (ndv) for WATTE
# Watte_ndv = 255  # No data value for WATTE
# Watte_bio_level[Watte_bio_level == ndv] = Watte_ndv  # 255 = No data value
# print('WATTE bio level min and max:', Watte_bio_level.min(), Watte_bio_level.max())
# create_raster('Productivity.tif', rasterHC, Watte_bio_level, gdal.GDT_Byte,
#               Watte_ndv)  # Productivity level gdal.GDT_Byte is unsigned 8 bit integer (0 to 255)
#
# # Inundation_depth = tb_update.copy()
# Inundation_depth = np.where((mhwIDW != ndv) & ((mhwIDW - tb_update) > 0) & (tb_update > -0.5), mhwIDW - tb_update,
#                             0)  # due to rasterization even bathymetry region yields inundation depth so add (tb_update > -0.5) to remove some bugs. However, this is an arbitary number! Need to consider further (Jin June 19, 2024)
# Inundation_depth[mhwIDW == ndv] = ndv
# # Inundation_depth = np.where((mhwIDW != ndv) & ((mhwIDW - tb_update) > 0), mhwIDW - tb_update, 0) # due to rasterization even bathymetry region yields inundation depth so add (tb_update > -0.5) to remove some bugs. # Grand Bay
#
# create_raster('Inundation_depth.tif', rasterHC, Inundation_depth, gdal.GDT_Float32, ndv,
#               stats_flag=True)  # Inundation depth gdal.GDT_Float32 is 32 bit floating point
# ############################################################################################################
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

########################################################################################################################
#--- Read data as scatter points ---
# print("   Processing scatter points...\n")
# myFile=open("filteredTidalDatumsIDW.pts","r")
# myLines=myFile.readlines();
# numPoints=len(myLines);
# pts=np.ones((numPoints,21),dtype=float);
# myFile.seek(0);
# for j in range(numPoints):
#     myLine=myFile.readline();
#     myRow=myLine.split();
#     pts[j][0]=int(myRow[0]);
#     pts[j][1]=int(myRow[1]);
#     pts[j][2]=float(myRow[2]);
#     pts[j][3]=float(myRow[3]);
#     pts[j][4]=float(myRow[4]);
#     pts[j][5]=float(myRow[5]);
#     pts[j][6]=float(myRow[6]);
#     pts[j][7]=float(myRow[7]);
#     pts[j][8]=float(myRow[8]);
#     pts[j][9]=float(myRow[9]);
#     pts[j][10]=float(myRow[10]);
#     pts[j][11]=float(myRow[11]);
#     pts[j][12]=float(myRow[12]);
#     pts[j][13]=float(myRow[13]);
#     pts[j][14]=float(myRow[14]);
#     pts[j][15]=float(myRow[15]);
#     pts[j][16]=float("NaN");
#     pts[j][17]=float("NaN");
#     pts[j][18]=float("NaN");
#     pts[j][19]=float("NaN");
#     pts[j][20]=float("NaN")
# myFile.close()

# #--- Read mesh with z-values as NWI indices ---
# print("\n"); print("   Processing mesh with z-values as NWI indices...\n");
# myFile=open("fort.NWI.14","r")
# myLine=myFile.readline()
# myLine=myFile.readline();
# myRow=myLine.split();
# numNodes=int(myRow[1]);
# nn=np.ones((numNodes,1),dtype=int);
# nx=np.ones((numNodes,1),dtype=float);
# ny=np.ones((numNodes,1),dtype=float);
# nw=np.ones((numNodes,1),dtype=float);
# for j in tqdm(range(numNodes)):
#     myLine=myFile.readline();
#     myRow=myLine.split();
#     N=int(myRow[0]);
#     X=float(myRow[1]);
#     Y=float(myRow[2]);
#     W=-float(myRow[3]);
#     nn[j][0]=N;
#     nx[j][0]=X;
#     ny[j][0]=Y;
#     nw[j][0]=W;
# myFile.close()

# #--- Biomass and accretion ---
# print("\n"); print("   Performing arithmetic calculations for biomass and accretion...\n");
#
# a1l=32.0; b1l=-3.2; c1l=1920.0; a1r=6.61; b1r=-0.661; c1r=1983.0; D01=5.0; # Example of Grand Bay (Alizad et al., 2018; PLOS ONE)
# a2l=73.8; b2l=-1.14; c2l=1587.1; a2r=a2l; b2r=b2l; c2r=c2l; D02=0.0; # Example of Weeks Bay (Alizad et al., 2018; PLOS ONE)
# a3l=197.5; b3l=-9870.0; c3l=1999.0; a3r=326.5; b3r=-1633.0; c3r=1998.0; D03=0.1; # Example of Apalachicola Bay (Alizad et al., 2016; Earth's Future)
# q=0.0018; k=0.000015; dt=10.0; # Coefficients taken from Morris et al. (2002; Ecology)
#
# D=np.ones((numPoints,1),dtype=float);
# Dcm=np.ones((numPoints,1),dtype=float);
# B1l=np.ones((numPoints,1),dtype=float);
# B1r=np.ones((numPoints,1),dtype=float);
# B1=np.ones((numPoints,1),dtype=float);
# B2l=np.ones((numPoints,1),dtype=float);
# B2r=np.ones((numPoints,1),dtype=float);
# B2=np.ones((numPoints,1),dtype=float);
# B3l=np.ones((numPoints,1),dtype=float);
# B3r=np.ones((numPoints,1),dtype=float);
# B3=np.ones((numPoints,1),dtype=float);
# A1=np.ones((numPoints,1),dtype=float);
# A2=np.ones((numPoints,1),dtype=float);
# A3=np.ones((numPoints,1),dtype=float);
#
# D[:,0]=pts[:,15]-pts[:,6];
# D[D<0.0]=0.0;
# Dcm=100.0*D;
# B1l=a1l*Dcm+b1l*Dcm**2+c1l;
# B1l[B1l<0.0]=0.0; B1l[Dcm>=D01]=0.0;
# B1r=a1r*Dcm+b1r*Dcm**2+c1r;
# B1r[B1r<0.0]=0.0;
# B1r[Dcm<D01]=0.0;
# B1=B1l+B1r;
#
# B2l=a2l*Dcm+b2l*Dcm**2+c2l;
# B2l[B2l<0.0]=0.0; B2l[Dcm>=D02]=0.0;
# B2r=a2r*Dcm+b2r*Dcm**2+c2r;
# B2r[B2r<0.0]=0.0; B2r[Dcm<D02]=0.0;
# B2=B2l+B2r;
# B3l=a3l*D+b3l*D**2+c3l;
# B3l[B3l<0.0]=0.0;
# B3l[D>=D03]=0.0;
# B3r=a3r*D+b3r*D**2+c3r;
# B3r[B3r<0.0]=0.0;
# B3r[D<D03]=0.0;
# B3=B3l+B3r;
# A1=q*D+k*B1*D;
# A1[D<0.0]=0.0;
# A1=A1*dt;
# A2=q*D+k*B2*D;
# A2=A2*dt;
# A2[D<0.0]=0.0;
# A3=q*D+k*B3*D;
# A3[D<0.0]=0.0;
# A3=A3*dt;
#
# #--- Assign values per wetland type ---
# print("   Assigning values per wetland type...\n")
# NWI1=8;
# NWI2=9;
# NWI3=20; # [8; 9; 20] = [salt marsh (regularly flooded); mangrove; irregularly flooded marsh]
#
# for j in range(numPoints):
#     idx=int(pts[j][1]);
#     pts[j][16]=nw[idx-1][0];
#     if pts[j][16]==NWI1:
#         pts[j][17]=B1[j][0];
#         pts[j][18]=A1[j][0];
#     elif pts[j][16]==NWI2:
#         pts[j][17]=B2[j][0];
#         pts[j][18]=A2[j][0];
#     elif pts[j][16]==NWI3:
#         pts[j][17]=B3[j][0];
#         pts[j][18]=A3[j][0];
#     else:
#         pts[j][17]=B1[j][0];
#         pts[j][18]=A1[j][0];
#
# #--- Update z- and n-values ---
# print("   Updating z- and n-values...\n")
# BL=1000.0;
# BU=1800.0; # [BL; BU] = [lower-mid bound of biomass; mid-upper bound of biomass]
# nL=0.035;
# nM=0.05;
# nU=0.07; # [nL; nM; nU] = [lower Manning; middle Manning; upper Manning]
#
# for j in range(numPoints):
#     pts[j][19]=pts[j][6]+pts[j][18];
#     pts[j][20]=pts[j][7];
#     if pts[j][17]>0.001:
#         pts[j][20]=nL
#     if pts[j][17]>BL:
#         pts[j][20]=nM
#     if pts[j][17]>BU:
#         pts[j][20]=nU
#
# #--- Write output ---
# print("   Writing output...\n")
# myOut=open("filteredBiomassAccretion.pts","w")
# for j in range(numPoints):
#     myOut.write(str(int(pts[j][0]))+" ");
#     myOut.write(str(int(pts[j][1]))+" ");
#     myOut.write(str(pts[j][2])+" ");
#     myOut.write(str(pts[j][3])+" ");
#     myOut.write(str(pts[j][4])+" ");
#     myOut.write(str(pts[j][5])+" ");
#     myOut.write(str(pts[j][6])+" ");
#     myOut.write(str(pts[j][7])+" ");
#     myOut.write(str(pts[j][8])+" ");
#     myOut.write(str(pts[j][9])+" ");
#     myOut.write(str(pts[j][10])+" ");
#     myOut.write(str(pts[j][11])+" ");
#     myOut.write(str(pts[j][12])+" ");
#     myOut.write(str(pts[j][13])+" ");
#     myOut.write(str(pts[j][14])+" ");
#     myOut.write(str(pts[j][15])+" ");
#     myOut.write(str(int(pts[j][16]))+" ");
#     myOut.write(str(pts[j][17])+" ");
#     myOut.write(str(pts[j][18])+" ");
#     myOut.write(str(pts[j][19])+" ");
#     myOut.write(str(pts[j][20])+"\n")
# myOut.close()

#--- Exit script ---
print("EXIT: Existing script!\n")
end=time.time(); print ("Time elapsed (seconds):",end-start);


