# WEADS: Wetland Ecosystem and Accretion Dynamics Simulator
There are two versions: **raster-based** ("src_raster") and **point-based** ("src_point") **WEADS**. \
This code needs **ADCIRC input/output files** such as fort.13, fort.14, fort.53, and inundationtime.63.
Currently, raster-based WEADS need a private module, pyadcircmodules (see prerequisites for raster version). 

### Python Version
This code is confirmed above Python 3.9 using CI/CD pipeline (Github Action)

## Getting Started

### §1. Install virtual conda env and activation
Type: ***conda env create -f env.yml*** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; autopep8 is used for Python code formatter, flake8 is a linter to check code, and pytest and mock are code testing tools. These packages are not mandatory for running WEADS.  

### §2. Activate virtual env
Type: ***conda activate WEADS_env***

### §3 Create a package: WEADS
Type: ***pip install -e .*** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -e: editable mode (Preferred for package developers, unless using ***pip install .***)

## Contents of the package
**WEADS_Point** for point-based WEADS (Recommended) \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;source codes are located on "src_point" folder \
**WEADS_Raster** for raster-based WEADS (Optional) \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;source codes are located on "src_raster" folder
## Running module in the WEADS package
Show available commands: 
* ***WEADS_Point --help*** 


usage: WEADS_Point [-h] [--all] [--preprocessing] [--td] [--mem] [--postprocessing] [--inputMeshFile INPUTMESHFILE] [--inputAttrFile INPUTATTRFILE]
                   [--inputHarmonicsFile INPUTHARMONICSFILE] [--inputShapeFile INPUTSHAPEFILE] [--outputMEMFile OUTPUTMEMFILE] [--outputMeshFile OUTPUTMESHFILE]
                   [--outputAttrFile OUTPUTATTRFILE] [--inEPSG INEPSG] [--outEPSG OUTEPSG] [--deltaT DELTAT] [--slr SLR]
                   [--inputInundationtimeFile INPUTINUNDATIONTIMEFILE] [--inputvegetationFile INPUTVEGETATIONFILE] [--skipresample]

  **Options:** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     -h, --help            show this help message and exit \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --all                 Run all processing steps \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --preprocessing       Run preprocessing step \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --td                  Run tidal datum step \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --mem                 Run mem step \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --postprocessing      Run postprocessing step \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     --inputMeshFile INPUTMESHFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                         Input mesh <fort.14 file> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputAttrFile INPUTATTRFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Input attribute <fort.13 file> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputHarmonicsFile INPUTHARMONICSFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Input harmonics <fort.53 file> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputShapeFile INPUTSHAPEFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Input domain shape file <xxx.shp> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --outputMEMFile OUTPUTMEMFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Output MEM file:ecology.csv \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --outputMeshFile OUTPUTMESHFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Output mesh file <fort_new.14> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --outputAttrFile OUTPUTATTRFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Output attribute file \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inEPSG INEPSG       Input EPSG code <inEPSGCode> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --outEPSG OUTEPSG     Output EPSG code <outEPSGCode> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --deltaT DELTAT       time step for WEADS simuiation in years \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --slr SLR             Sea level rise \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputInundationtimeFile INPUTINUNDATIONTIMEFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Use inundationtime file for running inunT <inundationtime.63> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputvegetationFile INPUTVEGETATIONFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Path to vegetation file <*.tif> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --skipresample        Skip reprojection and resample to raster domain 

# Example of a command for point-based WEADS 
***WEADS_Point --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputShapeFile Cut_domain.shp --inEPSG 4269 --outEPSG 26914 --deltaT 25 --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --inputInundationtimeFile inundationtime.63 --inputvegetationFile NWI_TX_wetlands.tif --all***

## Prerequisites for WEADS_Point
The user needs to download ***singularity_image*** 

* [ADCIRC Modules](https://github.com/zcobell/ADCIRCModules) -> **singularity image. chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/singularity_image_sm3/adcircmodules_docker_updated.sif**

## Example of a command for raster-based WEADS
Interactive job on SuperMike3

1. "srun -t 2:00:00 -n8 -N1 -A hpc_ceds2d_1hp -p single --pty /bin/bash" (singularity command is only available on work nodes)
2. ""singularity exec -B /home /project/jinikeda/adcircmodules_docker_updated.sif conda env create -f env.yml"" (Need to be confirmed)
3. Copy ADCIRCModules folder on your work directory (e.g., /work/jinikeda/ETC/TCB/WEAD/ADCIRCModules)
4. cd /work/jinikeda/ETC/TCB/WEAD/ADCIRCModules/testing/python_tests
5. Copy fort.13, fort.14, fort.53, everdried.63, domain shapfile(aaa.shp etc)
6. "singularity exec -B /work /project/jinikeda/adcircmodules_docker_updated.sif python3 hydromem.py --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputEverdriedFile everdried.63 --inputShapeFile Cut_grd.shp --inEPSG 4269 --outEPSG 26919 --gridSize 200 --outputMEMRasterFile MEM --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --all"

### Example: 100 years ADCIRC-WEADS simulation on Texas Coastal Bend 

https://github.com/user-attachments/assets/8a179637-1f4e-40b4-a53c-53641dece49e
<p style="text-align: left;"><strong>Figure.1</strong> 100 years ADCIRC-WEADS simulation under intermediate sea level rise scenario on the Coastal Texas Bend</p>

## Contributors
* Jin Ikeda (LSU|Center for Computation and Technology)
* Peter Bacopoulos (LSU|Coastal Ecosystem Design Studio)
* [Christopher E Kees (LSU|Coastal Ecosystem Design Studio Director)](https://www.lsu.edu/ceds/) 

## Acknolwedgements
* [Matthew V. Bilskie (UGA)](https://coast.engr.uga.edu/)
* Scott C. Hagen (emeritus LSU)
* Karim Alizad
* [Zachary Cobell](https://thewaterinstitute.org/our-team/zachary-cobell)
