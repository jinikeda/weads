#!/usr/bin/env python
# File: Ecology.py
# Developer: Jin Ikeda & Peter Bacopoulos
# Last modified: April 24, 2025


# ----------------------------------------------------------
# MEM: Marsh Equilibrium Model and Mangrove Development
# ----------------------------------------------------------
#
# Uses Inverse Distance Weighting to interpolate/extrapolote
# tidal datums outside fully wet areas.

# result = function(interpolateHarmonicsFile, inputvegetationFile, outputMEMFile + '.csv',inEPSG,outEPSG, deltaT=deltaT)
# ----------------------------------------------------------

##########################################################################
# --- Load internal modules ---
from .general_functions import *
from dataclasses import dataclass, field

# --- Initialize code ---
start_time = time.time()
print("\n")
print("LAUNCH: Launching script!\n")

# ----------------------------------------------------------
# GLOSSARY OF VARIABLES #
# ----------------------------------------------------------
#
# A: Annual Accretion [m yr -1]
# B: Biomass density [g m-2 yr -1]
# Bl: Biomass density left side
# Br: Biomass density right side
# Bmax: Maximum biomass production [g m-2 yr -1]
# tb_update: modified topo-bathy [m]
# P: Productivity [-], 16: low, 23: medium, 32: high for tidal marsh
# Dmin, Dmax, Dopt: Minimum, maximum, and optimal depth for biomass production in parabola equation

# ----------------------------------------------------------
# DATASETS Jin -> Pete Add States and Optimum March 7, 2024
# ----------------------------------------------------------
# Biomass curve coefficients Grand Bay  Alizad et al. (2018)


interest_reference = "Texas_Coastal_Bend"

# Use a dictionary to store the coefficients
biomass_coefficients = {
    # can be include, "Bmax": 2400},
    "NorthInlet": {"al": 1000, "bl": -3718, "cl": 1021, "ar": 1000, "br": -3718, "cr": 1021, "Dopt": 22},
    "Apalachicola": {"al": 1.975, "bl": -0.987, "cl": 1999, "ar": 3.265, "br": -1.633, "cr": 1998},
    "WeeksBay": {"al": 73.8, "bl": -1.14, "cl": 1587.1, "ar": 73.8, "br": -1.14, "cr": 1587.1},
    "GrandBay": {"al": 32, "bl": -3.2, "cl": 1920, "ar": 6.61, "br": -0.661, "cr": 1983, "Dopt": 5.0},
    "PlumIsland": {"al": 24.96, "bl": -0.193, "cl": 592.7, "ar": 24.96, "br": -0.193, "cr": 592.7},
    "Texas_Coastal_Bend": {"al": 240.0, "bl": -5.0, "cl": -460.0, "ar": 240.0, "br": -5.0, "cr": -460.0, "Dopt": 22.0},
    "Texas_Coastal_Bend_mangrove": {"al": 1600, "bl": -17.7, "cl": -28016.0, "ar": 1600, "br": -17.7, "cr": -28016.0, "Dopt": 45.0}}

# Need to add Optimum Elevation at least to run calculate_biomass_parabola

##########################################################################
print("Input parameters for MEM")
##########################################################################

# --- LOCAL ACCRETION PARAMETERS ---
# when vegetationFile == True: modify the following part
# print('The accretion parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')

##########################################################################
# Under consideration part Jin Aug 22, 2024
# We may need to modify the following part using the class
##########################################################################

# class Vegetation:
#     def __init__(self, Dmin, Dmax, Dopt, Bmax, Kr, RRS, BTR):
#         self.Dmin = Dmin
#         self.Dopt = Dopt
#         self.Dmax = Dmax
#         self.Bmax = Bmax
#         self.Kr = Kr
#         self.RRS = RRS
#         self.BTR = BTR
#
# SaltMarsh = Vegetation(2.0, 46.0, 2400.0, 22.0, 0.1, 2.0, 0.5)
# Mangrove = Vegetation(24.0, 66.0, 7800.0, 45.0, 0.1, 1.8, 0.25)
# IrregularMarsh = Vegetation(-70.0, -10.0, 1200.0, -40.0, 0.1, 1.5, 0.5)

##########################################################################
# Locarized parameters
##########################################################################
vegetation_parameters = {
    "SaltMarsh": {
        "Dmin": 2.0,
        "Dopt": 22.0,
        "Dmax": 46.0,
        "Bmax": 2400.0,
        "Kr": 0.1,
        "RRS": 2.0,
        "BTR": 0.5
    },
    "Mangrove": {
        "Dmin": 24.0,
        "Dopt": 45.0,
        "Dmax": 66.0,
        "Bmax": 7800.0,
        "Kr": 0.1,
        "RRS": 1.8,
        "BTR": 0.25
    },
    "IrregularMarsh": {
        "Dmin": -70.0,  # -21.0 -> 0.5*std below the mean of the marsh
        "Dopt": -40.0,  # 7.0 # modified based on tb and vegetation mapping
        "Dmax": -10.0,  # 35.0 -> 0.5*std above the mean of the marsh
        "Bmax": 1200.0,
        "Kr": 0.1,
        "RRS": 1.5,
        "BTR": 0.5
    },
    ##########################################################################
    # Need future Revision Aug 21, 2024, Jin
    ##########################################################################
    # "juvenile_mangrove": {
    #     "Dmin": 26.0,
    #     "Dopt": 34.0,
    #     "Dmax": 47.0,
    #     "Bmax": 1200.0,
    #     "Kr": 0.1,
    #     "RRS": 1.8,
    #     "BTR": 0.25
    # },
    ##########################################################################
}

BTRmat = 0.22  # BG turnover rate (1/year) for mature mangroves
BTRjuv = 0.67  # BG turnover rate (1/year) for juvenile mangroves
Tmat = 30  # Time for pioneer mangroves to fully mature (yr)

##########################################################################

@dataclass
class MEMConfig:
    ndv: float =  -99999.0  # No data value (ndv) using ADCIRC conversion
    ndv_byte: int = 128
    qmax: float =2.8     # Maximum capture coefficient (-)
    qmin: float = 1.0    # Minimum capture coefficient (-)
    SSC: float = 25.     # Suspended sediment concentration (mg/L) default 25.0
    FF: float = 353.0    # Flooding frequency (1/year) default 353.0
    BDo: float = 0.085   # Organic inputs # # Bulk density of organic matter (g/cm3)
    BDi: float = 1.99    # Mineral inputs # Bulk density of inorganic matter (g/cm3)
    Kr: float = 0.1      # Refractory fraction (g/g) default 0.1:
    RRS: float = 2.0     # Below Ground Bio to Shoot Ratio (g/g) default 2.0
    BTR: float = 0.5     # Below Ground turnover rate (/yr) default 0.5
    biomass_coefficients: dict = field(default_factory=lambda: biomass_coefficients)
    vegetation_parameters: dict = field(default_factory=lambda: vegetation_parameters)
    interest_reference: str = interest_reference

# ----------------------------------------------------------
# F U N C T I O N
# ----------------------------------------------------------


def calculate_coefficients(Dmin, Dmax, Dopt, Bmax):
    a = -((-Dmin * Bmax - Dmax * Bmax) /
          ((Dmin - Dopt) * (Dopt - Dmax)))  # coefficient of D
    b = -(Bmax / ((Dmin - Dopt) * (Dopt - Dmax)))  # coefficient of D^2
    c = -Dmin * Bmax * Dmax / ((Dmin - Dopt) * (Dopt - Dmax))  # constant term
    return a, b, c


def calculate_biomass_parabola(D, Dopt, mask, al, bl, cl, ar, br, cr):

    config = MEMConfig()  # call the configuration class

    # --- BIOMASS CALCULATIONS ---
    # Create a mask for the condition
    mask_left_parabora = (D <= Dopt) & mask
    mask_right_parabora = (Dopt < D) & mask

    # Perform the calculations only where the condition is true
    Bl = np.where(mask_left_parabora, al * D + bl * D * D + cl, config.ndv)
    Br = np.where(mask_right_parabora, ar * D + br * D * D + cr, config.ndv)

    print('BL:', Bl.min(), Bl.max())
    print('BR:', Br.min(), Br.max())

    # Set negative values to negative small value
    Bl[Bl < 0.0] = -0.001
    Br[Br < 0.0] = -0.001

    B = Bl + Br
    B[B < 0] = config.ndv

    Bmax = B.max()
    print('B max and min', Bmax, B.min())
    PL = Bmax / 3
    PH = Bmax * 2 / 3

    return B, PL, PH


def calculate_vertical_accretion(
        qmin, qmax, dt, B, Bmax, Kr, RRS, BTR, SSC, FF, D, Dt, BDo, BDi, tb):

    config = MEMConfig()  # call the configuration class

    # --- ACCRETION CALCULATIONS ---
    q = qmin + (qmax - qmin) * B / Bmax
    q[B < 0.0] = qmin
    Vorg = Kr * RRS * BTR * B / (100.0**2)  # organic matter

    # Calculate tidal flood time (FIT)  (Vahsen et al., 2024).
    FIT = np.full(D.shape,1,dtype=float)
    FIT[(D > 0) & (Dt > 0.0)] = D[(D > 0) & (Dt > 0.0)] / Dt[(D > 0) & (Dt > 0.0)]
    FIT[(FIT>=1)] = 1.0  # Tidal range is less than few mm should not calculate the FIT
    FIT[D <= 0] = 0
    print(f" FIT min:\t{np.min (FIT)}, max\t{np.max(FIT)}")
    # print(f" Dt min:\t{np.min(Dt[(D > 0) & (Dt > 0.0)])}, max\t{np.max(Dt[(D > 0) & (Dt > 0.0)])}")
    # print(f" D min:\t{np.min(D[(D > 0) & (Dt > 0.0)])}, max\t{np.max(D[(D > 0) & (Dt > 0.0)])}")

    # Assert that all values in FIT should be range of 0 and 1
    assert np.all((FIT >= 0) & (FIT <= 1)), "Some values in FIT are not between 0 and 1"

    Vmin = 0.5 * FIT * q * SSC * FF * D / \
        (1000.0**2)  # inorganic matter q2 add FIT
    Vorg[B <= 0.0] = 0.0
    # If we also evaluate the accretion rate for mineral matter, we will
    # change the code using mask such as (mhwIDW != ndv)
    Vmin[np.logical_or(B <= 0.0, D <= 0.0)] = 0.0
    A = np.where(tb != config.ndv, ((Vorg / BDo) + (Vmin / BDi)) /
                 100.0, config.ndv)  # accretion rate per year [m]
    # accretion thickness [m]
    tb_update = np.where(tb != config.ndv, tb + (A * dt), config.ndv)

    return tb_update, A


def calculate_manning(P):

    config = MEMConfig()  # call the configuration class

    if 10 <= P < 20:  # 16: low productivity
        return 0.035
    elif 20 <= P < 28:  # 23 medium productivity
        return 0.05
    elif 28 <= P < 38:  # 32 high productivitiy
        return 0.07
    elif 100 <= P < 110:  # 109 mangrove (juvenile)
        return 0.07  # Estuarine Scrub/Shrub Wetland 0.07 (juvenile), Estuarine Forested Wetland n =0.15 (mature), Passeri et al 2018
    elif 120 <= P < 130:  # 120 irregularly flooded marsh
        return 0.07 # Estuarine Scrub/Shrub Wetland 0.07
    else:
        return config.ndv

##########################################################################


def mem(interpolateHarmonicsFile, vegetationFile,
        outputMEMFile, inEPSG, outEPSG, deltaT=5):

    config = MEMConfig()  # call the configuration class

    # define from WEADS_Point.py
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

    if vegetationFile is None:
        print('\nA monotypic species with no vegetation mapping\n')
        # subscript:sub-optimal(left) and super-optimal(right) branches that
        # met at the parabola apex

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
        # Morris et al. (2016) Contributions of organic and inorganic matter to
        # sedimentvolume and accretion in tidal wetlands at steady state

    else:
        print('\nMulti species vegetation mapping\n')

        print("Read a vegetation file")
        df_VG = pd.read_csv(vegetationFile)
        Point_VG = df_VG['NWI'].values

        print(np.max(Point_VG), np.min(Point_VG), df_VG.shape)

        unique_values, counts = np.unique(Point_VG, return_counts=True)
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
        updated_unique_values = np.delete(
            unique_values,
            np.where(
                np.isin(
                    unique_values,
                    values_to_remove)))

        print('Interest unique_value is :', updated_unique_values)

        # Create a mask for different vegetation types
        mask_salt_marsh = (Point_VG == 8)
        mask_mangrove = (Point_VG == 9)
        mask_irregular = (Point_VG == 20) | (Point_VG == 8)  # Since overwrap by salt marshes, include mask_salt_marsh into this mask

        print('The biomass parameters for salt marsh (8), mangrove (9) and irregularly flooded marsh (20)')

        # 1: Salt marsh (NWI = 8)
        vegetation_type = "SaltMarsh"
        params = vegetation_parameters[vegetation_type]
        Dmin_marsh = params["Dmin"]
        Dmax_marsh = params["Dmax"]
        Dopt_marsh = params["Dopt"]
        Bmax_marsh = params["Bmax"]
        kr_marsh = params["Kr"]
        RRS_marsh = params["RRS"]
        BTR_masrh = params["BTR"]

        a1, b1, c1 = calculate_coefficients(
            Dmin_marsh, Dmax_marsh, Dopt_marsh, Bmax_marsh)  # 1: Salt marsh (NWI = 8)
        print('Check: SaltMarsh abc', a1, b1, c1)

        # 2: Mangrove (NWI = 9) (mangrove (mature))
        vegetation_type = "Mangrove"
        params = vegetation_parameters[vegetation_type]
        Dmin_matmang = params["Dmin"]
        Dmax_matmang = params["Dmax"]
        Dopt_matmang = params["Dopt"]
        Bmax_matmang = params["Bmax"]
        kr_matmang = params["Kr"]
        RRS_matmang = params["RRS"]
        BTR_matmang = params["BTR"]

        a2, b2, c2 = calculate_coefficients(
            Dmin_matmang, Dmax_matmang, Dopt_matmang, Bmax_matmang)
        print('Check: Mangrove abc', a2, b2, c2)

        # 3: Irregularly flooded marsh (NWI = 20)
        vegetation_type = "IrregularMarsh"
        params = vegetation_parameters[vegetation_type]
        Dmin_fmarsh = params["Dmin"]
        Dmax_fmarsh = params["Dmax"]
        Dopt_fmarsh = params["Dopt"]
        Bmax_fmarsh = params["Bmax"]
        kr_fmarsh = params["Kr"]
        RRS_fmarsh = params["RRS"]
        BTR_fmarsh = params["BTR"]

        a3, b3, c3 = calculate_coefficients(
            Dmin_fmarsh, Dmax_fmarsh, Dopt_fmarsh, Bmax_fmarsh)

##########################################################################
        # Juvenile mangroves
### Jin is pending the following part to modify  #########################

        # vegetation_type = "Mangrove_Juvenile"
        # params = vegetation_parameters[vegetation_type]
        # Dmin_matmang = params["Dmin"]
        # Dmax_matmang = params["Dmax"]
        # Dopt_matmang = params["Dopt"]
        # Bmax_matmang = params["Bmax"]
        # kr_matmang = params["Kr"]
        # RRS_matmang = params["RRS"]
        # BTR_matmang = params["BTR"]
        #
        # a4, b4, c4 = calculate_coefficients(Dmin_juvmang, Dmax_juvmang, Dopt_juvmang, Bmax_juvmang)

##########################################################################

    # ----------------------------------------------------------
    # Read point files
    # ----------------------------------------------------------

    PRJ = "EPSG:" + str(outEPSG)

    # --- READ INPUTS ---
    # Read the CSV file
    df = pd.read_csv(interpolateHarmonicsFile)
    print("  Read HyControl (HC) and tidal datums (TD) successfully")
    print(df.shape, df.columns, df.dtypes)
    if vegetationFile:
        assert Point_VG.shape != df_VG.shape[0], "The shape of the vegetation file is not matched with the tidal file"

    # 'HydroClass' 0: land, 1: intertidal, 2: subtidal(water)
    hc = df['HydroClass']

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
    # submergence region (hc = 1.0, 2.0, 3.0 including lake and pond) # these
    # values are not existed current point-based code
    submergence_mask = water_mask | intertidal_mask | (2.5 < hc)
    above_subtidal_zone = (mhwIDW != config.ndv)  # above subtidal zone ()
    # np.savetxt('above_subtidal_zone.csv', above_subtidal_zone, delimiter=',')

    # Perform the calculation only where the condition is true
    # [cm] Relative depth at nan should be positive value in this case
    D = np.where(above_subtidal_zone, 100.0 * (mhwIDW - tb), -config.ndv)
    D[D > -config.ndv / 2] = config.ndv  # Relative depth [cm] re-modify to negative value for ndv
    # Create a mask for the elements of D that are not equal to ndv
    D_mask = (D != config.ndv)
    # Print the minimum and maximum of the elements of D that are not equal to
    # ndv
    print(f"Depth Min:\t{D[D_mask].min()}, Max:\t{D[D_mask].max()}")

    Dt = 100.0 * (mhwIDW - mlwIDW)  # Tidal range [cm]
    Dt[np.logical_or(mhwIDW == config.ndv, mlwIDW == config.ndv)] = config.ndv
    print('Tidal range: min and max [cm]\t', Dt.min(), Dt.max())

    # --- PERFORM HYDRO-MEM CALCULATIONS ---
    print("Starting Hydro-MEM Calculations")

    if vegetationFile is None:

        P = np.full((df.shape[0], 1), config.ndv_byte,
                    dtype=int)  # Create an array of default values (ndv)
        marsh = np.full((df.shape[0], 1), config.ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        print("Biomass Calculations")
        B, PL, PH = calculate_biomass_parabola(
            D, Dopt, above_subtidal_zone, al, bl, cl, ar, br, cr)
        print('Salt marsh ', 'PH = ', PH, ', PL = ', PL, 'B max = ', B.max())

        # --- ACCRETION CALCULATIONS ---
        print("Accretion calculations")
        tb_update, A = calculate_vertical_accretion(
            config.qmin, config.qmax, dt, B, B.max(), config.Kr, config.RRS, config.BTR, config.SSC, config.FF, D, Dt, config.BDo, config.BDi, tb)

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
        Dmax = -(al / (2 * bl))
        Dzero1 = (-al + np.sqrt(((al * al) - (4 * bl * cl)))) / (2 * bl)
        Dzero2 = (-ar - np.sqrt(((ar * ar) - (4 * br * cr)))) / (2 * br)
        DRange = abs(Dzero2 - Dzero1)
        DHigh = Dmax + DRange * 0.1
        DLow = Dmax - DRange * 0.1

        condition_1 = (Dzero1 < D) & (D <= DLow)
        condition_2 = (DLow < D) & (D < DHigh)
        condition_3 = (DHigh <= D) & (D < Dzero2)

        marsh[condition_1] = 30
        marsh[condition_2] = 20
        marsh[condition_3] = 10

    else:

        # ##### Here check the elevation of mask irregularly flooded marsh ######
        # # This part is for tuning the parameters for irregularly flooded marsh (not mandatory
        # # # Before dilation
        # df_VG_Org = pd.read_csv('domain_nwi_original.csv')  # Before dilation
        # NWI_original = df_VG_Org['NWI'].values
        # mask_irregular_NWI = (NWI_original == 20)
        # irregular_elevation = tb[mask_irregular_NWI & above_subtidal_zone]
        # irregular_NWI_depth = D[mask_irregular_NWI & above_subtidal_zone]
        # print(
        #     'Irregular_depth [cm], min, max, mean, std',
        #     irregular_NWI_depth.min(),
        #     irregular_NWI_depth.max(),
        #     irregular_NWI_depth.mean(),
        #     irregular_NWI_depth.std())
        # ######################################################################

        P = np.full((df.shape[0], 1), config.ndv,
                    dtype=float)  # Create an array of default values (ndv)
        marsh = np.full((df.shape[0], 1), config.ndv, dtype=float)
        mangrove = np.full((df.shape[0], 1), config.ndv, dtype=float)
        irregular = np.full((df.shape[0], 1), config.ndv, dtype=float)

        # --- BIOMASS CALCULATIONS ---
        # Salt marsh (8)
        B1, PL1, PH1 = calculate_biomass_parabola(
            D, Dopt_marsh, above_subtidal_zone, a1, b1, c1, a1, b1, c1)
        print('Salt marsh ', 'PH = ', PH1, ', PL = ', PL1)

        # Mangrove (9)
        B2, PL2, PH2 = calculate_biomass_parabola(
            D, Dopt_matmang, above_subtidal_zone, a2, b2, c2, a2, b2, c2)
        print('Mangrove ', 'PH = ', PH2, ', PL = ', PL2)

        # Irregularly marsh (20)
        B3, PL3, PH3 = calculate_biomass_parabola(
            D, Dopt_fmarsh, above_subtidal_zone, a3, b3, c3, a3, b3, c3)
        print('Irregularly marsh ', 'PH = ', PH3, ', PL =', PL3)

        # --- ACCRETION CALCULATIONS ---
        _, A1 = calculate_vertical_accretion(
            config.qmin, config.qmax, dt, B1, Bmax_marsh, kr_marsh, RRS_marsh, BTR_masrh, config.SSC, config.FF, D, Dt, config.BDo, config.BDi, tb)
        _, A2 = calculate_vertical_accretion(
            config.qmin, config.qmax, dt, B2, Bmax_matmang, kr_matmang, RRS_matmang, BTR_matmang, config.SSC, config.FF, D, Dt, config.BDo, config.BDi, tb)
        _, A3 = calculate_vertical_accretion(
            config.qmin, config.qmax, dt, B3, Bmax_fmarsh, kr_fmarsh, RRS_fmarsh, BTR_fmarsh, config.SSC, config.FF, D, Dt, config.BDo, config.BDi, tb)

        # --- PRODUCTIVITY CALCULATIONS ---
        # Create a mask for different conditions
        mask_regular_1 = mask_salt_marsh & (0 < B1) & (B1 <= PL1)
        mask_regular_2 = mask_salt_marsh & (PL1 < B1) & (B1 < PH1)
        mask_regular_3 = mask_salt_marsh & (PH1 <= B1)
        # & (B2 <= PH2) # This part keep for future revision (0 < B2) & (B2 <= PH2)
        mask_mangrove = mask_mangrove & (0 < B2)
        mask_irregular = mask_irregular & (0 < B3)  # & (B3 <= PH3)
        # mask_mangrove_mat = mask_mangrove & future revision
        # mask_mangrove_juv = mask_mangrove & future revision
        # mask_mangrove = mask_mangrove #& (B2 >= PH2)

        # Assign values based on conditions
        # From background raster first. here ndv value will not overwrite the
        # background
        P[land_mask] = 55
        P[submergence_mask] = 40

        # Marsh productivity (Caution for the oder of the mask)
        P[mask_irregular] = 120  # irregularly flooded marsh
        P[mask_regular_1] = 16  # low productivity
        P[mask_regular_2] = 23  # medium productivity
        P[mask_regular_3] = 32  # high productivity
        P[mask_mangrove] = 109  # mangrove

        # --- MARSH TYPE CALCULATIONS ---
        print("High-Low Marsh Calculations")

        al = a1
        ar = a1
        bl = b1
        br = b1
        cl = c1
        cr = c1

        Dmax = -(al / (2 * bl))
        Dzero1 = (-al + np.sqrt(((al * al) - (4 * bl * cl)))) / (2 * bl)
        Dzero2 = (-ar - np.sqrt(((ar * ar) - (4 * br * cr)))) / (2 * br)
        DRange = abs(Dzero2 - Dzero1)
        DHigh = Dmax + DRange * 0.1
        DLow = Dmax - DRange * 0.1

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

        # Mask ndv values in the arrays
        # mask_marsh may compete with mangrove due to the dilution order, so
        # evaluate both masks here.
        B1_masked = np.where((mask_salt_marsh | mask_mangrove), B1, config.ndv)
        B2_masked = np.where(mask_mangrove, B2, config.ndv)
        B3_masked = np.where(
            (mask_salt_marsh | mask_mangrove | mask_irregular), B3, config.ndv)

        # background of entire grid is set at 0 [m]
        A1_masked = np.where((mask_salt_marsh | mask_mangrove), A1, 0)
        A2_masked = np.where(mask_mangrove, A2, 0)
        A3_masked = np.where(
            (mask_salt_marsh | mask_mangrove | mask_irregular), A3, 0)

        print('Amax', A1_masked.max(), A2_masked.max(), A3_masked.max())

        # Merge the arrays (used maximum value for each cell)
        # Maximum biomass density on each grid [g m-2 yr -1]
        B = np.maximum.reduce([B1_masked, B2_masked, B3_masked])
        # Maximum accretion rate per year on each grid [m]
        A = np.maximum.reduce([A1_masked, A2_masked, A3_masked])
        A[tb == config.ndv] = config.ndv  # Overwrite the outside domain as  the ndv values

        # Update topo-bathy data
        tb_update = np.where(
            (mask_salt_marsh | mask_mangrove | mask_irregular), tb + (A * dt), tb)
        # update topo-bathy data due to some holes are existed inside the
        # domain
        tb_update[tb == config.ndv] = config.ndv

    ##########################################################################
    # --- WRITE OUTPUTS ---
    print("")
    print("Writing point base output")

    df['D'] = D.flatten()  # mhwIDW - tb [cm]
    df['B'] = B.flatten()
    df['A'] = A.flatten()
    df['tb_update'] = tb_update.flatten()
    df['marsh'] = marsh.flatten()
    df['P'] = P.flatten()
    df['manning'] = df['P'].apply(calculate_manning)

    ###### Renew vegetation map ####################
    print("Renew the vegetation map from the PRODUCTIVITY CALCULATIONS")
    VM = P.copy()
    if vegetationFile is not None:
        # 8 = salt marsh(regularly flooded) follow with tidal cycle
        VM[(VM == 16) | (VM == 23) | (VM == 32)] = 8
        VM[(VM == 109)] = 9  # 9 = mangrove
        VM[(VM == 120)] = 20  # 20 = irregularly flooded marsh
        # VM[(VM == 55)] = ndv_byte # 55 = land_mask
        # VM[(VM == 40)] = ndv_byte # 40 = water_mask

    df['new_NWI'] = VM.flatten()



    ####################################################################################################################
    # --- WRITE OUTPUTS (optional) ---
    print("\n----------------------------------------------------------------------")
    print("Create WATTE input files")
    print("----------------------------------------------------------------------\n")
    ####################################################################################################################

    # Output
    Watte_bio_level = VM.copy()  # Create an array of default values (ndv) for WATTE
    Watte_ndv = 255  # No data value for WATTE
    Watte_bio_level[Watte_bio_level == config.ndv] = Watte_ndv  # 255 = No data value
    print('WATTE bio level min and max:', Watte_bio_level.min(), Watte_bio_level.max())

    df['bio_level'] = Watte_bio_level.flatten()

    # Inundation_depth = tb_update.copy()
    Inundation_depth = np.where((mhwIDW != config.ndv) & ((mhwIDW - tb_update) > 0) & (tb_update > -0.5), mhwIDW - tb_update,
                                0)  # due to rasterization even bathymetry region yields inundation depth so add (tb_update > -0.5) to remove some bugs. However, this is an arbitary number! Need to consider further (Jin June 19, 2024)
    Inundation_depth[mhwIDW == config.ndv] = config.ndv
    # Inundation_depth = np.where((mhwIDW != ndv) & ((mhwIDW - tb_update) > 0), mhwIDW - tb_update, 0) # due to rasterization even bathymetry region yields inundation depth so add (tb_update > -0.5) to remove some bugs. # Grand Bay
    df['inun_depth'] = Inundation_depth.flatten()

    df.to_csv(outputMEMFile, index=False)
    ############################################################################################################
    '''
    print ("Output file written successfully")
    print ("")

    # --- PRINT TO SCREEN ---
    print ("")
    print ("-------------------------------------------------")
    print ("--- COMPLETED HYDRO-MEM IN PYTHON - HYDRO MEM ---")
    print ("-------------------------------------------------")
    print ("")
    print ("")
    '''

    ##########################################################################
    # Calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the elapsed time
    print("Done ecological response")
    print("Time to Compute: \t\t\t", elapsed_time, " seconds")
