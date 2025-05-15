# File: src_point/mem_delft3d4.py
# Author: Shabnam
# Modified: June 2025

import time
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from .basics import fileexists

# --- Configuration and Constants ---

interest_reference = "Texas_Coastal_Bend"

biomass_coefficients = {
    "Texas_Coastal_Bend": {
        "al": 240.0, "bl": -5.0, "cl": -460.0,
        "ar": 240.0, "br": -5.0, "cr": -460.0,
        "Dopt": 22.0
    },
    "Texas_Coastal_Bend_mangrove": {
        "al": 1600, "bl": -17.7, "cl": -28016.0,
        "ar": 1600, "br": -17.7, "cr": -28016.0,
        "Dopt": 45.0
    }
}

vegetation_parameters = {
    "SaltMarsh": {
        "Dmin": 2.0, "Dopt": 22.0, "Dmax": 46.0,
        "Bmax": 2400.0, "Kr": 0.1, "RRS": 2.0, "BTR": 0.5
    },
    "Mangrove": {
        "Dmin": 24.0, "Dopt": 45.0, "Dmax": 66.0,
        "Bmax": 7800.0, "Kr": 0.1, "RRS": 1.8, "BTR": 0.25
    },
    "IrregularMarsh": {
        "Dmin": -70.0, "Dopt": -40.0, "Dmax": -10.0,
        "Bmax": 1200.0, "Kr": 0.1, "RRS": 1.5, "BTR": 0.5
    }
}

@dataclass
class MEMConfig:
    ndv: float = -99999.0
    ndv_byte: int = 128
    qmax: float = 2.8
    qmin: float = 1.0
    SSC: float = 25.0
    FF: float = 353.0
    BDo: float = 0.085
    BDi: float = 1.99
    Kr: float = 0.1
    RRS: float = 2.0
    BTR: float = 0.5
    biomass_coefficients: dict = field(default_factory=lambda: biomass_coefficients)
    vegetation_parameters: dict = field(default_factory=lambda: vegetation_parameters)
    interest_reference: str = interest_reference

# --- Core Functions ---

def calculate_coefficients(Dmin, Dmax, Dopt, Bmax):
    a = -((-Dmin * Bmax - Dmax * Bmax) / ((Dmin - Dopt) * (Dopt - Dmax)))
    b = -(Bmax / ((Dmin - Dopt) * (Dopt - Dmax)))
    c = -Dmin * Bmax * Dmax / ((Dmin - Dopt) * (Dopt - Dmax))
    return a, b, c

def calculate_biomass_parabola(D, Dopt, mask, al, bl, cl, ar, br, cr):
    config = MEMConfig()
    Bl = np.where((D <= Dopt) & mask, al * D + bl * D**2 + cl, config.ndv)
    Br = np.where((D > Dopt) & mask, ar * D + br * D**2 + cr, config.ndv)
    Bl[Bl < 0.0] = -0.001
    Br[Br < 0.0] = -0.001
    B = Bl + Br
    B[B < 0] = config.ndv

    Dmax = -(al / (2 * bl))
    Dzero1 = (-al + np.sqrt(al**2 - 4 * bl * cl)) / (2 * bl)
    Dzero2 = (-ar - np.sqrt(ar**2 - 4 * br * cr)) / (2 * br)
    DRange = abs(Dzero2 - Dzero1)
    DHigh = Dmax + DRange * 0.1
    DLow = Dmax - DRange * 0.1

    B_valid = B[B > 0]
    PL_val = np.percentile(B_valid, 25) if B_valid.size > 0 else config.ndv
    PH_val = np.percentile(B_valid, 75) if B_valid.size > 0 else config.ndv

    PL = np.where((D > 0) & (D <= DLow), PL_val, config.ndv)
    PH = np.where((D > 0) & (D >= DHigh), PH_val, config.ndv)

    return B, PL, PH

def calculate_vertical_accretion(qmin, qmax, dt, B, Bmax, Kr, RRS, BTR, SSC, FF, D, Dt, BDo, BDi, tb):
    config = MEMConfig()
    q = qmin + (qmax - qmin) * B / Bmax
    q[B < 0.0] = qmin
    Vorg = Kr * RRS * BTR * B / 10000.0
    FIT = np.full(D.shape, 1.0, dtype=float)
    FIT[(D > 0) & (Dt > 0)] = D[(D > 0) & (Dt > 0)] / Dt[(D > 0) & (Dt > 0)]
    FIT[FIT >= 1] = 1.0
    FIT[D <= 0] = 0
    Vmin = 0.5 * FIT * q * SSC * FF * D / 1e6
    Vorg[B <= 0.0] = 0.0
    Vmin[np.logical_or(B <= 0.0, D <= 0.0)] = 0.0
    A = np.where(tb != config.ndv, (Vorg / BDo + Vmin / BDi) / 100.0, config.ndv)
    tb_update = np.where(tb != config.ndv, tb + A * dt, config.ndv)
    return tb_update, A

def calculate_manning(P):
    config = MEMConfig()
    if 10 <= P < 20:
        return 0.035
    elif 20 <= P < 28:
        return 0.05
    elif 28 <= P < 38:
        return 0.07
    elif 100 <= P < 110:
        return 0.07
    elif 120 <= P < 130:
        return 0.07
    else:
        return config.ndv

# --- Main MEM Delft3D Function ---

def mem_delft(domainIOFile, vegetationFile, outputMEMFile, inEPSG, outEPSG, deltaT=5):
    start_time = time.time()
    print("\nLAUNCH: Running MEM step for Delft3D coupling\n")

    config = MEMConfig()
    dt = deltaT

    fileexists(domainIOFile)
    df = pd.read_csv(domainIOFile)
    print("Read domain input file:", df.shape)

    hc = df['HydroClass']
    tb = df['z']
    mhw = df['MHW_IDW']
    mlw = df['MLW_IDW']

    land_mask = (hc == 0)
    intertidal_mask = (hc == 1)
    water_mask = (hc == 2)
    submergence_mask = water_mask | intertidal_mask
    above_subtidal_zone = (mhw != config.ndv)

    D = np.where(above_subtidal_zone, 100.0 * (mhw - tb), -config.ndv)
    D[D > -config.ndv / 2] = config.ndv
    D_mask = (D != config.ndv)

    Dt = 100.0 * (mhw - mlw)
    Dt[np.logical_or(mhw == config.ndv, mlw == config.ndv)] = config.ndv

    print(f"Depth range [cm]: min={D[D_mask].min()}, max={D[D_mask].max()}")

    if vegetationFile is None:
        print("\nMonotypic case: No vegetation file provided")
        coeffs = config.biomass_coefficients[config.interest_reference]
        al, bl, cl = coeffs['al'], coeffs['bl'], coeffs['cl']
        ar, br, cr = coeffs['ar'], coeffs['br'], coeffs['cr']
        Dopt = coeffs['Dopt']

        B, PL, PH = calculate_biomass_parabola(D, Dopt, above_subtidal_zone, al, bl, cl, ar, br, cr)
        tb_update, A = calculate_vertical_accretion(config.qmin, config.qmax, dt, B, B.max(), config.Kr, config.RRS, config.BTR, config.SSC, config.FF, D, Dt, config.BDo, config.BDi, tb)

        P = np.full(df.shape[0], config.ndv_byte, dtype=int)
        P[land_mask] = 55
        P[submergence_mask] = 40
        P[(B > 0) & (B <= PL)] = 16
        P[(B > PL) & (B < PH)] = 23
        P[B >= PH] = 32

        marsh = np.full(df.shape[0], config.ndv, dtype=float)
        Dmax = -(al / (2 * bl))
        Dzero1 = (-al + np.sqrt(al**2 - 4 * bl * cl)) / (2 * bl)
        Dzero2 = (-ar - np.sqrt(ar**2 - 4 * br * cr)) / (2 * br)
        DRange = abs(Dzero2 - Dzero1)
        DHigh = Dmax + DRange * 0.1
        DLow = Dmax - DRange * 0.1

        marsh[(Dzero1 < D) & (D <= DLow)] = 30
        marsh[(DLow < D) & (D < DHigh)] = 20
        marsh[(DHigh <= D) & (D < Dzero2)] = 10

    else:
        raise NotImplementedError("Vegetation raster input not yet implemented.")

    df['D'] = D
    df['B'] = B
    df['A'] = A
    df['tb_update'] = tb_update
    df['marsh'] = marsh
    df['P'] = P
    df['manning'] = df['P'].apply(calculate_manning)

    VM = P.copy()
    VM[(VM == 16) | (VM == 23) | (VM == 32)] = 8
    df['new_NWI'] = VM

    df.to_csv(outputMEMFile, index=False)
    print("\nMEM step completed.")
    print("Time elapsed: {:.2f} seconds".format(time.time() - start_time))
