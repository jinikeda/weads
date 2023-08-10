# pyHydroMEM
Python Hydro-MEM code

## Python Version
This code is being developed using Python 3.7.5

## Getting Started

```python hydromem.py --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputEverdriedFile everdried.63 --inputShapeFile PIE_bbox.shp --inEPSG 4269 --outEPSG 26919 --gridSize 200 --outputMEMRasterFile mem --all```

```--all``` will run all (--rasterize, --td, and --mem)

```--rasterize``` will only create the raster files from ADCIRC inputs

```--td``` will only run the tidal datums calculation

```--mem``` will only run mem

## Prerequisites

What things you need to install the software and how to install them

* [ADCIRC Modules](https://github.com/zcobell/ADCIRCModules)
* See requirements.txt

## Docker Image
To enhance the east of use, a docker image for all pre-requisites is available
https://hub.docker.com/repository/docker/mbilskie/pyhydromem/

## Contributors
* [Karim Alizad (USGS)](https://www.usgs.gov/centers/spcmsc)
* [Peter Bacopoulos (LSU)](https://www.lsu.edu/ccr/)
* [Matthew V. Bilskie (UGA)](https://coast.engr.uga.edu/)
* [Davina L. Passeri (USGS)](https://www.usgs.gov/staff-profiles/davina-l-passeri?qt-staff_profile_science_products=0#qt-staff_profile_science_products)
* Eric Swanson (USGS)
* [Scott C. Hagen (LSU)](https://www.lsu.edu/ccr/)

## Acknolwedgements
* [Zachary Cobell](https://thewaterinstitute.org/our-team/zachary-cobell)
