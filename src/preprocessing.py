#!/usr/bin/python3
# File: preprocessing.py
# Name: Peter Bacopoulos
# Date: August 10, 2023

#--- Load modules ---
import shapefile
from pyproj import Transformer
import numpy as np

#--- Read bounding box ---
print("\n"); print("Processing bounding box...\n");
inEPSG=4269; outEPSG=26919;
transformer=Transformer.from_crs(inEPSG,outEPSG)
sf=shapefile.Reader("GEO_bbox.shp")
bboxGeo=sf.shape(0).bbox; bbox=[None]*4;
bbox[0],bbox[1]=transformer.transform(bboxGeo[1],bboxGeo[0])
bbox[2],bbox[3]=transformer.transform(bboxGeo[3],bboxGeo[2])

#--- Read mesh ---
print("Processing mesh...\n")
myFile=open("fort.14","r")
myLine=myFile.readline()
myLine=myFile.readline(); myRow=myLine.split(); numNodes=int(myRow[1]);
nn=np.ones((numNodes,1),dtype=int);
nx=np.ones((numNodes,1),dtype=float);
ny=np.ones((numNodes,1),dtype=float);
nz=np.ones((numNodes,1),dtype=float);
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split();
    N=int(myRow[0]); X=float(myRow[1]); Y=float(myRow[2]); Z=-float(myRow[3]);
    X,Y=transformer.transform(Y,X); nn[j][0]=N; nx[j][0]=X; ny[j][0]=Y; nz[j][0]=Z;
myFile.close()

#--- Read attributes ---
print("Processing attributes...\n")
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
for j in range(numAttributes):
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
print("Processing everdried...\n")
myFile=open("everdried.63","r")
myLine=myFile.readline()
myLine=myFile.readline()
myLine=myFile.readline()
ed=np.ones((numNodes,1),dtype=float)
for j in range(numNodes):
    myLine=myFile.readline(); myRow=myLine.split(); ed[j][0]=float(myRow[1]);
myFile.close()

# --- Read harmonics ---
print("Processing harmonics...\n")
myFile=open("fort.53","r")
myLine=myFile.readline(); myRow=myLine.split(); numHarm=int(myRow[0]);
harmAMP=np.ones((numNodes,numHarm),dtype=float); harmPHA=np.ones((numNodes,numHarm),dtype=float);
for j in range(numHarm):
    myLine=myFile.readline()
myLine=myFile.readline()
for j in range(numNodes):
    myLine=myFile.readline()
    for jj in range(numHarm):
        myLine=myFile.readline(); myRow=myLine.split();
        harmAMP[j][jj]=float(myRow[0]); harmPHA[j][jj]=float(myRow[1]);
myFile.close()

#--- Filter for bounding box ---
print("Filtering data for bounding box...\n")
myOut1=open("filteredData.pts","w")
myOut2=open("filteredHarmAmp.pts","w"); myOut3=open("filteredHarmPha.pts","w");
cnt=0
for j in range(numNodes):
    if nx[j][0]>bbox[0] and nx[j][0]<bbox[2]:
        if ny[j][0]>bbox[1] and ny[j][0]<bbox[3]:
            cnt=cnt+1
            myOut1.write(str(cnt)+" "); myOut1.write(str(j)+" ");
            myOut1.write(str(nx[j][0])+" "); myOut1.write(str(ny[j][0])+" ");
            myOut1.write(str(nz[j][0])+" "); myOut1.write(str(mann[j][0])+" ");
            myOut1.write(str(ed[j][0])+"\n");
            for jj in range(numHarm-1):
                myOut2.write(str(harmAMP[j][jj])+" "); myOut3.write(str(harmPHA[j][jj])+" ");
            myOut2.write(str(harmAMP[j][jj+1])+"\n"); myOut3.write(str(harmPHA[j][jj+1])+"\n");
myOut1.close(); myOut2.close(); myOut3.close();


