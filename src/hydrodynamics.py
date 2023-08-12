#!/usr/bin/python3
# File: hydrodynamics.py
# Name: Peter Bacopoulos
# Date: August 11, 2023

#--- Load modules ---
import numpy as np
import math
import sys

#--- Read data as scatter points ---
print("\n"); print("Processing scatter points...\n");
myFile=open("filteredData.pts","r")
myLines=myFile.readlines(); numPoints=len(myLines);
pts=np.ones((numPoints,9),dtype=float); myFile.seek(0);
for j in range(numPoints):
    myLine=myFile.readline(); myRow=myLine.split();
    pts[j][0]=int(myRow[0]); pts[j][1]=int(myRow[1]);
    pts[j][2]=float(myRow[2]); pts[j][3]=float(myRow[3]);
    pts[j][4]=float(myRow[4]); pts[j][5]=float(myRow[5]);
    pts[j][6]=float(myRow[6]); pts[j][7]=float(myRow[7]);
    pts[j][8]=float(myRow[8])
myFile.close()

#--- Read harmomics correspondent with scatter ---
print("Processing harmonics...\n")
myFile1=open("filteredHarmAmp.pts","r")
myFile2=open("filteredHarmPha.pts","r")
myFile3=open("harmonics.freq","r")
myLine=myFile1.readline(); myRow=myLine.split(); numHarm=len(myRow); myFile1.seek(0);
harmAMP=np.ones((numPoints,numHarm),dtype=float); harmPHA=np.ones((numPoints,numHarm),dtype=float);
freq=np.ones((numHarm,1),dtype=float)
for j in range(numHarm):
    myLine3=myFile3.readline(); myRow3=myLine3.split(); freq[j][0]=float(myRow3[0]);
for j in range(numPoints):
    myLine1=myFile1.readline(); myLine2=myFile2.readline();
    myRow1=myLine1.split(); myRow2=myLine2.split();
    for jj in range(numHarm):
        harmAMP[j][jj]=float(myRow1[jj]); harmPHA[j][jj]=float(myRow2[jj]); harmPHA[j][jj]=math.radians(harmPHA[j][jj]);
myFile1.close(); myFile2.close(); myFile3.close();

#--- Calculate tidal datums for fully wetted zones ---
print("Calculating tidal datums for fully wetted zones...\n")
T=30.0; dt=900.0; N=int((86400*T/dt)+1); t=dt*np.ones((N,1),dtype=float); t=np.cumsum(t,axis=0);
T2=2.0*T; D=24.0*3600.0/dt; D2=D/2.0; wl2=np.ones((int(D2),1),dtype=float); wl2L=np.ones((int(T2),1),dtype=float); wl2H=np.ones((int(T2),1),dtype=float);
mlw=np.ones((numPoints,1),dtype=float); mhw=np.ones((numPoints,1),dtype=float); msl=np.ones((numPoints,1),dtype=float);
for j in range(numPoints):
    if pts[j][8]>0.0:
        #--- Tidal resynthesis ---
        wl=harmAMP[j][0]*np.ones((N,1),dtype=float);
        for jj in range(1,numHarm):
            wl=wl+harmAMP[j][jj]*np.cos(freq[jj][0]*t-harmPHA[j][jj])
        #--- Tidal datums ---
        for k in range(int(T2)):
            for kk in range(int(D2)):
                wl2[kk][0]=wl[int(D2)*k+kk+1][0]
            wl2L[k][0]=min(wl2); wl2H[k][0]=max(wl2);
        msl[j][0]=np.mean(wl); mlw[j][0]=np.mean(wl2L); mhw[j][0]=np.mean(wl2H);
    else:
        msl[j][0]=-99999.0; mlw[j][0]=-99999.0; mhw[j][0]=-99999.0;

myOut=open("filteredTidalDatums.pts","w")
for j in range(numPoints):
    myOut.write(str(pts[j][0])+" "); myOut.write(str(pts[j][1])+" ");
    myOut.write(str(pts[j][2])+" "); myOut.write(str(pts[j][3])+" ");
    myOut.write(str(pts[j][4])+" "); myOut.write(str(pts[j][5])+" ");
    myOut.write(str(pts[j][6])+" "); myOut.write(str(pts[j][7])+" ");
    myOut.write(str(pts[j][8])+" "); myOut.write(str(mlw[j][0])+" ");
    myOut.write(str(msl[j][0])+" "); myOut.write(str(mhw[j][0])+"\n");
myOut.close()

#--- Interpolate tidal datums to intertidal zones ---
print("Interpolating tidal datums to intertidal zones...\n")

#--- Exit script ---
print("Exiting script...\n")
sys.exit()


