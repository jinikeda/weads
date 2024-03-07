# pyHydroMEM
Python Hydro-MEM code. This code needs a private module pyadcircmodules -> singularity image. Ask Jin

## Python Version
This code is confirmed using Python 3.10 with GDAL3.4.1

## Getting Started
Interactive job on SuperMike3

1. "srun -t 2:00:00 -n8 -N1 -A hpc_ceds2d_1hp -p single --pty /bin/bash" (singularity command is only available on work nodes)
2. "singularity exec -B /home /project/jinikeda/adcircmodules_docker_updated.sif pip install --user -r requirement.txt" (Need to be confirmed)
3. Copy ADCIRCModules folder on your work directory (e.g., /work/jinikeda/ETC/TCB/WEAD/ADCIRCModules)
4. cd /work/jinikeda/ETC/TCB/WEAD/ADCIRCModules/testing/python_tests
5. Copy fort.13, fort.14, fort.53, everdried.63, domain shapfile(aaa.shp etc)
6. "singularity exec -B /work /project/jinikeda/adcircmodules_docker_updated.sif python3 hydromem.py --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputEverdriedFile everdried.63 --inputShapeFile Cut_grd.shp --inEPSG 4269 --outEPSG 26919 --gridSize 200 --outputMEMRasterFile MEM --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --all"

Example of command
inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputEverdriedFile everdried.63 --inputShapeFile PIE_bbox.shp --inEPSG 4269 --outEPSG 26919 --gridSize 200 --outputMEMRasterFile MEM --all```

```--all``` will run all (--rasterize, --td, and --mem)

```--rasterize``` will only create the raster files from ADCIRC inputs

```--td``` will only run the tidal datums calculation

```--mem``` will only run mem

## Prerequisites

What things you need to install the software and how to install them

* [ADCIRC Modules](https://github.com/zcobell/ADCIRCModules) -> **singularity image. chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/singularity_image_sm3/adcircmodules_docker_updated.sif**
* See requirements.txt

## Contributors
Christopher E Kees
Peter Bacopoulos
Jin Ikeda

## Acknolwedgements
* [Karim Alizad (USGS)](https://www.usgs.gov/centers/spcmsc)
* [Peter Bacopoulos (LSU)](https://www.lsu.edu/ccr/)
* [Matthew V. Bilskie (UGA)](https://coast.engr.uga.edu/)
* [Davina L. Passeri (USGS)](https://www.usgs.gov/staff-profiles/davina-l-passeri?qt-staff_profile_science_products=0#qt-staff_profile_science_products)
* Eric Swanson (USGS)
* [Scott C. Hagen (LSU)](https://www.lsu.edu/ccr/)
* [Zachary Cobell](https://thewaterinstitute.org/our-team/zachary-cobell)
