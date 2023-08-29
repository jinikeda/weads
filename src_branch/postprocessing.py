#!/usr/bin/python3
# File: postprocessing.py
# Name: Peter Bacopoulos
# Date: August 23, 2023

#--- Load modules ---
import numpy as np
from tqdm import tqdm
import time

#--- Initialize code ---
start=time.time()
print("\n"); print("LAUNCH: Launching script!\n");

#--- Read data as scatter points ---
print("   Processing scatter points...\n")
myFile=open("filteredBiomassAccretion.pts","r")
myLines=myFile.readlines(); numPoints=len(myLines);
pts=np.ones((numPoints,21),dtype=float); myFile.seek(0);
for j in range(numPoints):
    myLine=myFile.readline(); myRow=myLine.split();
    pts[j][0]=int(myRow[0]); pts[j][1]=int(myRow[1]);
    pts[j][2]=float(myRow[2]); pts[j][3]=float(myRow[3]);
    pts[j][4]=float(myRow[4]); pts[j][5]=float(myRow[5]);
    pts[j][6]=float(myRow[6]); pts[j][7]=float(myRow[7]);
    pts[j][8]=float(myRow[8]); pts[j][9]=float(myRow[9]);
    pts[j][10]=float(myRow[10]); pts[j][11]=float(myRow[11]);
    pts[j][12]=float(myRow[12]); pts[j][13]=float(myRow[13]);
    pts[j][14]=float(myRow[14]); pts[j][15]=float(myRow[15]);
    pts[j][16]=float(myRow[16]); pts[j][17]=float(myRow[17]);
    pts[j][18]=float(myRow[18]); pts[j][19]=float(myRow[19]);
    pts[j][20]=float(myRow[20])
myFile.close()

#--- Read and write mesh ---
print("   Processing mesh...\n")
myFile=open("fort.14","r"); myOut=open("fort.dt.14","w");
myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split(); numNodes=int(myRow[1]);
nn=np.ones((numNodes,1),dtype=int);
nx=np.ones((numNodes,1),dtype=float);
ny=np.ones((numNodes,1),dtype=float);
nz=np.ones((numNodes,1),dtype=float);
print("      Reading mesh...\n")
for j in tqdm(range(numNodes)):
    myLine=myFile.readline(); myRow=myLine.split();
    N=int(myRow[0]); X=float(myRow[1]); Y=float(myRow[2]); Z=float(myRow[3]);
    nn[j][0]=N; nx[j][0]=X; ny[j][0]=Y; nz[j][0]=Z;
for j in range(numPoints):
    idx=int(pts[j][1]); nz[idx-1][0]=-pts[j][19];
print("\n"); print("      Writing mesh: Nodal table...\n");
for j in tqdm(range(numNodes)):
    myOut.write(str(nn[j][0])+" "); myOut.write(str(nx[j][0])+" ");
    myOut.write(str(ny[j][0])+" "); myOut.write(str(nz[j][0])+"\n");
print("\n"); print("      Writing mesh: Element connectivity and nodestring information...\n");
while True:
    myLine=myFile.readline()
    if not myLine:
        break
    myOut.write(myLine.strip()+"\n")
myFile.close(); myOut.close();

#--- Read everdried ---
print("   Processing everdried...\n")
myFile=open("everdried.63","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline()
ed=np.ones((numNodes,1),dtype=float)
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split(); ed[j][0]=float(myRow[1]);
myFile.close()

#--- Read and write attributes ---
print("   Processing attributes...\n")
nOW=0.022; # Manning's roughness for open-water conversion
myFile=open("fort.13","r"); myOut=open("fort.dt.13","w");
myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split(); numAttributes=int(myRow[0]);
for j in range(numAttributes):
    myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split();
    if myRow[0]=="mannings_n_at_sea_floor":
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split(); mann0=float(myRow[0]);
    else:
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
mann=np.ones((numNodes,2),dtype=float); mann[:,1]=mann[:,1]*mann0;
for j in range(numAttributes):
    myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split();
    if myRow[0]=="mannings_n_at_sea_floor":
        print("      Working bottom roughness...\n")
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split(); numMann=int(myRow[0]); cnt=np.ones((numMann,1),dtype=float);
        print("         Reading values...\n");
        for jj in tqdm(range(numMann)):
            myLine=myFile.readline(); myRow=myLine.split(); NN=int(myRow[0]); mann[NN-1][0]=float(NN); mann[NN-1][1]=float(myRow[1]); cnt[jj][0]=mann[NN-1][0];
        print("\n"); print("         Processing values...\n");
        for jj in range(numPoints):
            idx=int(pts[jj][1]); mann[idx-1][0]=pts[jj][20];
        print("         Writing values...\n")
        for jj in tqdm(range(numMann)):
            NN=int(cnt[jj][0])
            if ed[NN-1][0]>0.0 and mann[NN-1][1]>0.025:
                mann[NN-1][1]=nOW
            myOut.write(str(int(mann[NN-1][0]))+" "); myOut.write(str(mann[NN-1][1])+"\n");
    else:
        myLine=myFile.readline(); myOut.write(myLine.strip()+"\n"); myRow=myLine.split(); numAttr=int(myRow[0]);
        for jj in range(numAttr):
            myLine=myFile.readline(); myOut.write(myLine.strip()+"\n");
myFile.close(); myOut.close();

#--- Exit script ---
print("\n"); print("EXIT: Existing script!\n");
end=time.time(); print ("Time elapsed (seconds):",end-start);


