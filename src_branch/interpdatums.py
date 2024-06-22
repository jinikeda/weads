#!/usr/bin/python3
# File: interpdatums.py
# Name: Peter Bacopoulos
# Date: August 16, 2023

#--- Load modules ---
import numpy as np
import math
from tqdm import tqdm
import time

#--- Initialize code ---
start=time.time()
print("\n"); print("LAUNCH: Launching script!\n");

#--- Read data as scatter points ---
print("   Processing scatter points...\n")
myFile=open("filteredTidalDatums.pts","r")
myLines=myFile.readlines();
numPoints=len(myLines);
pts=np.ones((numPoints,16),dtype=float);
myFile.seek(0);
for j in tqdm(range(numPoints)):
    myLine=myFile.readline();
    myRow=myLine.split();
    pts[j][0]=int(myRow[0]);
    pts[j][1]=int(myRow[1]);
    pts[j][2]=float(myRow[2]);
    pts[j][3]=float(myRow[3]);
    pts[j][4]=float(myRow[4]);
    pts[j][5]=float(myRow[5]);
    pts[j][6]=float(myRow[6]);
    pts[j][7]=float(myRow[7]);
    pts[j][8]=float(myRow[8]);
    pts[j][9]=float(myRow[9]);
    pts[j][10]=float(myRow[10]);
    pts[j][11]=float(myRow[11]);
    pts[j][12]=float(myRow[12]);
    pts[j][13]=float(myRow[13]);
    pts[j][14]=float("NaN");
    pts[j][15]=float("NaN");
myFile.close()

#--- Interpolate datums ---
print("\n"); print("   Interpolating tidal datums from fully wetted zones to intertidal zones...\n");
itmax=np.max(pts,axis=0);
itmax=itmax[9];
itmin=np.min(pts,axis=0);
itmin=itmin[9];
listI=[];
listW=[];
for j in range(numPoints):
    if pts[j][10]>0.0:
        listI.append(pts[j][:])
    if pts[j][8]>0.0:
        listW.append(pts[j][:])
arrayI=np.array(listI);
arrayW=np.array(listW);
distI=np.ones((len(listI),3),dtype=float);
distW=np.ones((len(listW),3),dtype=float);
mlwI=np.ones((len(listI),3),dtype=float);
mhwI=np.ones((len(listI),3),dtype=float);
for j in tqdm(range(len(listI))):
    distW[:,0]=(arrayI[j][4]-arrayW[:,4])**2.0;
    distW[:,1]=(arrayI[j][5]-arrayW[:,5])**2.0;
    distW[:,2]=(distW[:,0]+distW[:,1])**0.5;
    arrayW[:,14]=distW[:,2];
    arrayW=arrayW[np.argsort(arrayW[:,14])]
    distI[j][0]=arrayW[0][14];
    distI[j][1]=arrayW[1][14];
    distI[j][2]=arrayW[2][14];
    distI[j][0]=1.0/distI[j][0]**2.0;
    distI[j][1]=1.0/distI[j][1]**2.0;
    distI[j][2]=1.0/distI[j][2]**2.0;
    #--- Mean low water ---
    mlwI[j][0]=arrayW[0][11];
    mlwI[j][1]=arrayW[1][11];
    mlwI[j][2]=arrayW[2][11];
    numer=distI[j][0]*mlwI[j][0]+distI[j][1]*mlwI[j][1]+distI[j][2]*mlwI[j][2]
    denom=distI[j][0]+distI[j][1]+distI[j][2]
    arrayI[j][14]=numer/denom
    #--- Mean high water ---
    mhwI[j][0]=arrayW[0][13];
    mhwI[j][1]=arrayW[1][13];
    mhwI[j][2]=arrayW[2][13];
    numer=distI[j][0]*mhwI[j][0]+distI[j][1]*mhwI[j][1]+distI[j][2]*mhwI[j][2]
    denom=distI[j][0]+distI[j][1]+distI[j][2]
    arrayI[j][15]=numer/denom

#--- Write output ---
print("\n"); print("   Writing output...\n");
myOut=open("filteredTidalDatumsIDW.pts","w")

# Jin's note! The following for-loops may something wrong for indent
for j in range(len(listI)):

#    myOut.write(str(int(arrayI[j][0]))+" "); myOut.write(str(int(arrayI[j][1]))+" ");

    myOut.write(str(j+1)+" ");
    myOut.write(str(int(arrayI[j][1]))+" ");
    myOut.write(str(arrayI[j][2])+" ");
    myOut.write(str(arrayI[j][3])+" ");
    myOut.write(str(arrayI[j][4])+" ");
    myOut.write(str(arrayI[j][5])+" ");
    myOut.write(str(arrayI[j][6])+" ");
    myOut.write(str(arrayI[j][7])+" ");
    myOut.write(str(arrayI[j][8])+" ");
    myOut.write(str(arrayI[j][9])+" ");
    myOut.write(str(arrayI[j][10])+" ");
    myOut.write(str(arrayI[j][11])+" ");
    myOut.write(str(arrayI[j][12])+" ");
    myOut.write(str(arrayI[j][13])+" ");
    myOut.write(str(arrayI[j][14])+" ");
    myOut.write(str(arrayI[j][15])+"\n")
myOut.close()

#--- Exit script ---
print("EXIT: Existing script!\n")
end=time.time(); print ("Time elapsed (seconds):",end-start);


