# WEADS: Wetland Ecosystem and Accretion Dynamics Simulator
There are two versions: **raster-based** ("src_raster") and **point-based** ("src_point") **WEADS**. \
This software will be required to run on High-performance computers or high-spec local machines (more than 8 cores CPU or equivalent with RAM = 16 GB).
This software also needs **ADCIRC input/output files** such as fort.13 (attributes), fort.14 (mesh), fort.53 (elevation harmonic constituents), and inundationtime.63 on the root directory (the same directory as setup.py).
Currently, raster-based WEADS need a dependency module, pyadcircmodules (https://github.com/zcobell/ADCIRCModules). For further details, see *prerequisites* for the raster version. 

### Python Version
The source codes are confirmed above Python 3.9 using CI/CD pipeline (Github Action)

## Getting Started

### §1. Install virtual conda env and activation
Type: ***conda env create -f env.yml*** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; autopep8 is used for Python code formatter, flake8 is a linter to check code, and pytest and mock are code testing tools. These packages are not mandatory for running WEADS.  

### §2. Activate virtual env
Type: ***conda activate WEADS_env***

### §3 Create a package: WEADS
Type: ***pip install -e .*** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -e: editable mode (Preferred for package developers, unless using ***pip install .***) \
Also, the user may need to set ***export PYTHONPATH=$PYTHONPATH:/path_dir_WEADS_xx.py*** and ***source ~/.bashrc** # or .bash_profile, .zshrc, etc* for the Path setting.
(We may move the argument files into the source folders to be determined)

### §4 Copy necessary files to run WEADS
Mandatory: fort.13, fort.14, fort.53, inundationtime.63, domain shapefile \
Optional: vegetation map file

## Contents of the package
**WEADS_Point** for point-based WEADS (** Highly Recommended**) \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;source codes are located on "src_point" folder \
**WEADS_Raster** for raster-based WEADS (Optional) \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;source codes are located on "src_raster" folder
## Running module in the WEADS package

* ***WEADS_Point --help***  Show available commands: 


usage: WEADS_Point [-h] [--all] [--preprocessing] [--td] [--mem] [--postprocessing] [--inputMeshFile INPUTMESHFILE] [--inputAttrFile INPUTATTRFILE]
                   [--inputHarmonicsFile INPUTHARMONICSFILE] [--inputShapeFile INPUTSHAPEFILE] [--outputMEMFile OUTPUTMEMFILE] [--outputMeshFile OUTPUTMESHFILE]
                   [--outputAttrFile OUTPUTATTRFILE] [--inEPSG INEPSG] [--outEPSG OUTEPSG] [--deltaT DELTAT] [--slr SLR]
                   [--inputInundationtimeFile INPUTINUNDATIONTIMEFILE] [--inputvegetationFile INPUTVEGETATIONFILE] [--skip_extracting_raster] [--no_spread_flag]

  **Options:** 

### key commands

```--all``` Run all procedure (--preprocessing, --td, --mem, and --postprocessing)

```--preprocessing``` Create csv files from ADCIRC inputs

```--td``` Calculation tidal datums calculation

```--mem``` Run ecological equiliubum model  

```--postprocessing``` Update ADCIRC input files (fort.13 and fort.14) from ecological productions and rasterization using the nearest neighbor interpolation (only for visualization purposes, not impact entire calculations)

### all commands 
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --deltaT DELTAT       time step for WEADS simulation in years \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --slr SLR             Sea level rise \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputInundationtimeFile INPUTINUNDATIONTIMEFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Use inundationtime file for running inunT <inundationtime.63> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --inputvegetationFile INPUTVEGETATIONFILE \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Path to vegetation file <*.tif> \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --skip_extracting_raster \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Skip extract original raster values \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  --no_spread_flag \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                        Turn off spread flag for vegetation mapping  

# Example of a command for point-based WEADS 
***WEADS_Point --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputShapeFile Study_domain.shp --inEPSG 4269 --outEPSG 26914 --deltaT 25 --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --inputInundationtimeFile inundationtime.63 --inputvegetationFile NWI_resample100.tif --all***

## Prerequisites for WEADS_Raster
The user needs to download ***singularity_image.sif*** to your machine \
GDAL 3.4.1, released 2021/12/27 is already installed in the singularity_image.sif

## Example of a command for raster-based WEADS
Interactive job on SuperMike3 (LSU | HPC cluster) \
**Cautions:** The work procedure strongly depends on your cluster environment. Please communicate with the cluster administrator.

1. *srun -t 2:00:00 -n8 -N1 -A your_allocation -p single --pty /bin/bash* (singularity command is only available on work nodes)
2. Install other Python packages *singularity exec -B /target-path /singularity-image-path/adcircmodules_docker_updated.sif pip install --user -r requirements.txt* 
3. Copy fort.13, fort.14, fort.53, inundationtime.63, domain shapefile(aaa.shp etc)
4. *singularity exec -B /work /singularity-image-path/adcircmodules_docker_updated.sif python3 WEADS_Raster.py --inputMeshFile fort.14 --inputAttrFile fort.13 --inputHarmonicsFile fort.53 --inputShapeFile Study_domain.shp --inEPSG 4269 --outEPSG 26914 --gridSize 100 --deltaT 25 --outputMEMRasterFile MEM --slr 0.0 --outputAttrFile fort_new.13 --outputMeshFile fort_new.14 --inputInundationtimeFile inundationtime.63 --inputvegetationFile NWI_resample100.tif --all*

For the details of the available command, type ***WEADS_Raster --help***

### Example: 100 years ADCIRC-WEADS simulation on Texas Coastal Bend 
https://github.com/user-attachments/assets/8a179637-1f4e-40b4-a53c-53641dece49e
<p style="text-align: left;"><strong>Figure.1</strong> 100 years ADCIRC-WEADS simulation under intermediate sea level rise scenario on the Coastal Texas Bend</p>

## Contributors
* Jin Ikeda (LSU|Center for Computation and Technology)
* Peter Bacopoulos (LSU|Coastal Ecosystem Design Studio)
* [Christopher E Kees (LSU|Coastal Ecosystem Design Studio Director)](https://www.lsu.edu/ceds/) 

## Acknowledgments
* [Matthew V. Bilskie (UGA)](https://coast.engr.uga.edu/)
* Scott C. Hagen (emeritus LSU)
* Karim Alizad
* [Zachary Cobell](https://thewaterinstitute.org/our-team/zachary-cobell)
