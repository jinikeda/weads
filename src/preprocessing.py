#!/usr/bin/python3
# File: preprocessing.py
# Name: Peter Bacopoulos
# Date: August 10, 2023

#--- Load modules ---
import numpy as np
import math
import sys
from tqdm import tqdm
import time

#--- Initialize code ---
start=time.time()
print("\n"); print("LAUNCH: Launching script!\n");

#--- Coordinate conversions (GEO-CPP) ---
print("   Processing coordinate conversion...\n")
re=6378206.4; lat=30.0;
kx=math.pi*re*math.cos(lat*math.pi/180.0)/180.0; kx=kx**1.0;
ky=math.pi*re/180.0; ky=ky**1.0;

#--- Read bounding box ---
print("   Processing bounding box...\n")
bbox=[None]*4; bbox[0]=-97.7; bbox[1]=27.3; bbox[2]=-96.8; bbox[3]=28.4;

#--- Read mesh ---
print("   Processing mesh...\n")
myFile=open("fort.14","r")
myLine=myFile.readline()
myLine=myFile.readline(); myRow=myLine.split(); numNodes=int(myRow[1]);
nn=np.ones((numNodes,1),dtype=int);
nx=np.ones((numNodes,1),dtype=float);
ny=np.ones((numNodes,1),dtype=float);
nz=np.ones((numNodes,1),dtype=float);
nxcpp=np.ones((numNodes,1),dtype=float);
nycpp=np.ones((numNodes,1),dtype=float);
for j in tqdm(range(numNodes)):
    myLine=myFile.readline(); myRow=myLine.split();
    N=int(myRow[0]); X=float(myRow[1]); Y=float(myRow[2]); Z=-float(myRow[3]);
    nn[j][0]=N; nx[j][0]=X; ny[j][0]=Y; nz[j][0]=Z;
    nxcpp[j][0]=X*kx; nycpp[j][0]=Y*ky;
myFile.close()

#--- Read attributes ---
print("\n"); print("   Processing attributes...\n");
myFile=open("fort.13","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline(); myRow=myLine.split(); numAttributes=int(myRow[0]);
for j in range(numAttributes):
    myLine=myFile.readline(); myRow=myLine.split();
    if myRow[0]=="mannings_n_at_sea_floor":
        myLine=myFile.readline()
        myLine=myFile.readline()
        myLine=myFile.readline(); myRow=myLine.split(); mann0=float(myRow[0]);
    else:
        myLine=myFile.readline()
        myLine=myFile.readline()
        myLine=myFile.readline()
mann=np.ones((numNodes,1),dtype=float); mann=mann*mann0;
for j in tqdm(range(numAttributes)):
    myLine=myFile.readline(); myRow=myLine.split();
    if myRow[0]=="mannings_n_at_sea_floor":
        myLine=myFile.readline(); myRow=myLine.split(); numMann=int(myRow[0]);
        for jj in range(numMann):
            myLine=myFile.readline(); myRow=myLine.split(); NN=int(myRow[0]); mann[NN-1][0]=float(myRow[1]);
    else:
        myLine=myFile.readline(); myRow=myLine.split(); numAttr=int(myRow[0]);
        for jj in range(numAttr):
            myLine=myFile.readline()
myFile.close()

#--- Read everdried ---
print("\n"); print("   Processing everdried...\n");
myFile=open("everdried.63","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline()
ed=np.ones((numNodes,1),dtype=float)
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split(); ed[j][0]=float(myRow[1]);
myFile.close()

#--- Read inundationtime ---
print("   Processing inundationtime...\n")
myFile=open("inundationtime.63","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline()
it=np.ones((numNodes,1),dtype=float)
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split(); it[j][0]=float(myRow[1]);
myFile.close()

#--- Read maxinundepth ---
print("   Processing maxinundepth...\n")
myFile=open("maxinundepth.63","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline()
mi=np.ones((numNodes,1),dtype=float)
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split(); mi[j][0]=float(myRow[1]);
myFile.close()

# --- Read harmonics ---
print("   Processing harmonics...\n")
myFile=open("fort.53","r"); myOut=open("harmonics.freq","w");
myLine=myFile.readline(); myRow=myLine.split(); numHarm=int(myRow[0]);
harmAMP=np.ones((numNodes,numHarm),dtype=float); harmPHA=np.ones((numNodes,numHarm),dtype=float);
for j in range(numHarm):
    myLine=myFile.readline(); myRow=myLine.split(); myOut.write(str(myRow[0])+"\n");
myLine=myFile.readline()
for j in tqdm(range(numNodes)):
    myLine=myFile.readline()
    for jj in range(numHarm):
        myLine=myFile.readline(); myRow=myLine.split();
        harmAMP[j][jj]=float(myRow[0]); harmPHA[j][jj]=float(myRow[1]);
myFile.close()

#--- Filter for bounding box ---
print("\n"); print("   Filtering data for bounding box...\n");
myOut1=open("filteredData.pts","w")
myOut2=open("filteredHarmAmp.pts","w"); myOut3=open("filteredHarmPha.pts","w");
cnt=0
for j in tqdm(range(numNodes)):
    if nx[j][0]>bbox[0] and nx[j][0]<bbox[2]:
        if ny[j][0]>bbox[1] and ny[j][0]<bbox[3]:
            cnt=cnt+1
            myOut1.write(str(cnt)+" "); myOut1.write(str(j+1)+" ");
            myOut1.write(str(nx[j][0])+" "); myOut1.write(str(ny[j][0])+" ");
            myOut1.write(str(nxcpp[j][0])+" "); myOut1.write(str(nycpp[j][0])+" ");
            myOut1.write(str(nz[j][0])+" "); myOut1.write(str(mann[j][0])+" ");
            myOut1.write(str(ed[j][0])+" "); myOut1.write(str(it[j][0])+" ");
            myOut1.write(str(mi[j][0])+"\n")
            for jj in range(numHarm-1):
                myOut2.write(str(harmAMP[j][jj])+" "); myOut3.write(str(harmPHA[j][jj])+" ");
            myOut2.write(str(harmAMP[j][jj+1])+"\n"); myOut3.write(str(harmPHA[j][jj+1])+"\n");
myOut1.close(); myOut2.close(); myOut3.close();

#--- Exit script ---
print("\n"); print("EXIT: Existing script!\n");
end=time.time(); print ("Time elapsed (seconds):",end-start);
print("\n"); sys.exit();


