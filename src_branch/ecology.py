#!/usr/bin/python3
# File: ecology.py
# Name: Peter Bacopoulos
# Date: August 21, 2023

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
myFile=open("filteredTidalDatumsIDW.pts","r")
myLines=myFile.readlines(); numPoints=len(myLines);
pts=np.ones((numPoints,21),dtype=float); myFile.seek(0);
for j in tqdm(range(numPoints)):
    myLine=myFile.readline(); myRow=myLine.split();
    pts[j][0]=int(myRow[0]); pts[j][1]=int(myRow[1]);
    pts[j][2]=float(myRow[2]); pts[j][3]=float(myRow[3]);
    pts[j][4]=float(myRow[4]); pts[j][5]=float(myRow[5]);
    pts[j][6]=float(myRow[6]); pts[j][7]=float(myRow[7]);
    pts[j][8]=float(myRow[8]); pts[j][9]=float(myRow[9]);
    pts[j][10]=float(myRow[10]); pts[j][11]=float(myRow[11]);
    pts[j][12]=float(myRow[12]); pts[j][13]=float(myRow[13]);
    pts[j][14]=float(myRow[14]); pts[j][15]=float(myRow[15]);
    pts[j][16]=float("NaN"); pts[j][17]=float("NaN");
    pts[j][18]=float("NaN"); pts[j][19]=float("NaN");
    pts[j][20]=float("NaN")
myFile.close()

#--- Read mesh with z-values as NWI indices ---
print("\n"); print("   Processing mesh with z-values as NWI indices...\n");
myFile=open("fort.NWI.14","r")
myLine=myFile.readline()
myLine=myFile.readline(); myRow=myLine.split(); numNodes=int(myRow[1]);
nn=np.ones((numNodes,1),dtype=int);
nx=np.ones((numNodes,1),dtype=float);
ny=np.ones((numNodes,1),dtype=float);
nw=np.ones((numNodes,1),dtype=float);
for j in tqdm(range(numNodes)):
    myLine=myFile.readline(); myRow=myLine.split();
    N=int(myRow[0]); X=float(myRow[1]); Y=float(myRow[2]); W=-float(myRow[3]);
    nn[j][0]=N; nx[j][0]=X; ny[j][0]=Y; nw[j][0]=W;
myFile.close()

#--- Biomass and accretion ---
print("\n"); print("   Performing arithmetic calculations for biomass and accretion...\n");

a1l=32.0; b1l=-3.2; c1l=1920.0; a1r=6.61; b1r=-0.661; c1r=1983.0; D01=5.0; # Example of Grand Bay (Alizad et al., 2018; PLOS ONE)
a2l=73.8; b2l=-1.14; c2l=1587.1; a2r=a2l; b2r=b2l; c2r=c2l; D02=0.0; # Example of Weeks Bay (Alizad et al., 2018; PLOS ONE)
a3l=197.5; b3l=-9870.0; c3l=1999.0; a3r=326.5; b3r=-1633.0; c3r=1998.0; D03=0.1; # Example of Apalachicola Bay (Alizad et al., 2016; Earth's Future)
q=0.0018; k=0.000015; dt=10.0; # Coefficients taken from Morris et al. (2002; Ecology)

D=np.ones((numPoints,1),dtype=float); Dcm=np.ones((numPoints,1),dtype=float);
B1l=np.ones((numPoints,1),dtype=float); B1r=np.ones((numPoints,1),dtype=float); B1=np.ones((numPoints,1),dtype=float);
B2l=np.ones((numPoints,1),dtype=float); B2r=np.ones((numPoints,1),dtype=float); B2=np.ones((numPoints,1),dtype=float);
B3l=np.ones((numPoints,1),dtype=float); B3r=np.ones((numPoints,1),dtype=float); B3=np.ones((numPoints,1),dtype=float);
A1=np.ones((numPoints,1),dtype=float); A2=np.ones((numPoints,1),dtype=float); A3=np.ones((numPoints,1),dtype=float);

D[:,0]=pts[:,15]-pts[:,6]; D[D<0.0]=0.0; Dcm=100.0*D;
B1l=a1l*Dcm+b1l*Dcm**2+c1l; B1l[B1l<0.0]=0.0; B1l[Dcm>=D01]=0.0; B1r=a1r*Dcm+b1r*Dcm**2+c1r; B1r[B1r<0.0]=0.0; B1r[Dcm<D01]=0.0; B1=B1l+B1r;
B2l=a2l*Dcm+b2l*Dcm**2+c2l; B2l[B2l<0.0]=0.0; B2l[Dcm>=D02]=0.0; B2r=a2r*Dcm+b2r*Dcm**2+c2r; B2r[B2r<0.0]=0.0; B2r[Dcm<D02]=0.0; B2=B2l+B2r;
B3l=a3l*D+b3l*D**2+c3l; B3l[B3l<0.0]=0.0; B3l[D>=D03]=0.0; B3r=a3r*D+b3r*D**2+c3r; B3r[B3r<0.0]=0.0; B3r[D<D03]=0.0; B3=B3l+B3r;
A1=q*D+k*B1*D; A1[D<0.0]=0.0; A1=A1*dt; A2=q*D+k*B2*D; A2=A2*dt; A2[D<0.0]=0.0; A3=q*D+k*B3*D; A3[D<0.0]=0.0; A3=A3*dt;

#--- Assign values per wetland type ---
print("   Assigning values per wetland type...\n")
NWI1=8; NWI2=9; NWI3=20; # [8; 9; 20] = [salt marsh (regularly flooded); mangrove; irregularly flooded marsh]
for j in range(numPoints):
    idx=int(pts[j][1]); pts[j][16]=nw[idx-1][0];
    if pts[j][16]==NWI1:
        pts[j][17]=B1[j][0]; pts[j][18]=A1[j][0];
    elif pts[j][16]==NWI2:
        pts[j][17]=B2[j][0]; pts[j][18]=A2[j][0];
    elif pts[j][16]==NWI3:
        pts[j][17]=B3[j][0]; pts[j][18]=A3[j][0];
    else:
        pts[j][17]=B1[j][0]; pts[j][18]=A1[j][0];

#--- Update z- and n-values ---
print("   Updating z- and n-values...\n")
BL=1000.0; BU=1800.0; # [BL; BU] = [lower-mid bound of biomass; mid-upper bound of biomass]
nL=0.035; nM=0.05; nU=0.07; # [nL; nM; nU] = [lower Manning; middle Manning; upper Manning]
for j in range(numPoints):
    pts[j][19]=pts[j][6]+pts[j][18]; pts[j][20]=pts[j][7];
    if pts[j][17]>0.001:
        pts[j][20]=nL
    if pts[j][17]>BL:
        pts[j][20]=nM
    if pts[j][17]>BU:
        pts[j][20]=nU

#--- Write output ---
print("   Writing output...\n")
myOut=open("filteredBiomassAccretion.pts","w")
for j in range(numPoints):
    myOut.write(str(int(pts[j][0]))+" "); myOut.write(str(int(pts[j][1]))+" ");
    myOut.write(str(pts[j][2])+" "); myOut.write(str(pts[j][3])+" ");
    myOut.write(str(pts[j][4])+" "); myOut.write(str(pts[j][5])+" ");
    myOut.write(str(pts[j][6])+" "); myOut.write(str(pts[j][7])+" ");
    myOut.write(str(pts[j][8])+" "); myOut.write(str(pts[j][9])+" ");
    myOut.write(str(pts[j][10])+" "); myOut.write(str(pts[j][11])+" ");
    myOut.write(str(pts[j][12])+" "); myOut.write(str(pts[j][13])+" ");
    myOut.write(str(pts[j][14])+" "); myOut.write(str(pts[j][15])+" ");
    myOut.write(str(int(pts[j][16]))+" "); myOut.write(str(pts[j][17])+" ");
    myOut.write(str(pts[j][18])+" "); myOut.write(str(pts[j][19])+" ");
    myOut.write(str(pts[j][20])+"\n")
myOut.close()

#--- Exit script ---
print("EXIT: Existing script!\n")
end=time.time(); print ("Time elapsed (seconds):",end-start);


