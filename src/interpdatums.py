#!/usr/bin/python3
# File: interpdatums.py
# Name: Peter Bacopoulos
# Date: August 16, 2023

#--- Load modules ---
import numpy as np
import math
import sys
from tqdm import tqdm
import time

#--- Initialize code ---
start=time.time()
print("\n"); print("LAUNCH: Launching script!\n");

#--- Read data as scatter points ---
print("Processing scatter points...\n")
myFile=open("filteredTidalDatums.pts","r")
myLines=myFile.readlines(); numPoints=len(myLines);
pts=np.ones((numPoints,13),dtype=float); myFile.seek(0);
for j in tqdm(range(numPoints)):
    myLine=myFile.readline(); myRow=myLine.split();
    pts[j][0]=int(myRow[0]); pts[j][1]=int(myRow[1]);
    pts[j][2]=float(myRow[2]); pts[j][3]=float(myRow[3]);
    pts[j][4]=float(myRow[4]); pts[j][5]=float(myRow[5]);
    pts[j][6]=float(myRow[6]); pts[j][7]=float(myRow[7]);
    pts[j][8]=float(myRow[8]); pts[j][9]=float(myRow[9]);
    pts[j][10]=float(myRow[10]); pts[j][11]=float(myRow[11]);
    pts[j][12]=float(myRow[12])
myFile.close()

#--- Interpolate datums ---
print("\n"); print("Interpolating tidal datums from fully wetted zones to intertidal zones...\n");
itmax=np.max(pts,axis=0); itmax=itmax[9];
itmin=np.min(pts,axis=0); itmin=itmin[9];
Xint=[]; Yint=[]; Xwet=[]; Ywet=[];
for j in range(numPoints):
    if pts[j][9]>itmin and pts[j][9]<itmax:
        Xint.append(pts[j][4]); Yint.append(pts[j][5]);
    if pts[j][8]>0.0:
        Xwet.append(pts[j][4]); Ywet.append(pts[j][5]);

#--- Exit script ---
print("\n"); print("EXIT: Existing script!\n")
end=time.time(); print ("Time elapsed (seconds):",end-start);
print("\n"); sys.exit();


