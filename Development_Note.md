# Development Note (Aug 22, 24 version)

-  **Discussed the development direction on June 12th**
---

## Input/Output data storage 
- Folder Path: chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO (need to update)


## Running file
- WEADS_Raster or WEADS_Point.py (Argument file to run WEAD)
  
## Main Python module/function files (src folders)
- grd2dem.py (Convert ADCIRC unstructured mesh file to tiff file)
- hyconn.py (Calculate hydro connectivity and classify ocean, land, and pond using ADCIRC output file: everdried.63)
- tidaldatums.py (Calculate tidal datums in each pixel on ocean region. **This part is still computationally expensive on Sep/12**)
- tidaldatumsidw.py (Interpolate the tidal datum on land and partially wet regions)
- <span style="color:blue"> **mem.py** </span> (Calculate vegetation productivities)
  - read_NWI_file.py (provide a multi-species domain with expansion: Jin completed on Sep 3)
- rast2adc.py (renew an ADCIRC mesh file: fort.14 based on mem simulation)
- update_nodal_attributes.py (renew an ADCIRC nodal attribute file: fort.13 based on mem simulation)
  - manning.py (Specify Manning's roughnesses based on biomass productivities but probably doesn't use this file)
  
## Supplemental files 
- raster.py (rasterization functions)
- basics.py (check existing a file)

# To-do-list

## Development (mandatory)
**src_point will be main source code and eliminate src_raster approach**
  - rasterize (elimiate sigularity image file)
**src/mem.py**
  - The only problem with reading vegetation classification is that a segmentation fault happens when reading a large tiff file with a docker.<br>input file: Folder Path/inputs/Wetlands_NWI_forMeshRefinement/Region3_NWI_LC_Reclassify_wetlandsOnly.tif,<br> output file: Folder Path/outputs/Domain_classification_distribution_resample100.tif NWI.py is a currently not elegant.
  

## Development (desirable)
- Compare Morris's Excel VBA sheet (Pete: https://github.com/cekees/weads/tree/main/sandbox)
- Modify average_horizontal_eddy_viscosity_in_sea_water_wrt_depth based on update fort.14
- Consider the climate and catastrophic aspects of WEADS development (Shabnam starts climate analysis and literature review).
- Pete and Shabnam work on soil cohorts using the NUMAR model approach.

## Modifications (desirable)

## Modifications (optional)
- Read Netcdf outputs

# NOTE

## Books
- A Blue carbon primer: the state of coastal wetland carbon science, practice and policy
- Another one is not sure

  
NWI classification
- 8 = salt marsh (regularly flooded) -> follow with tidal cycle
- 9 = mangrove
- 20 = irregularly flooded marsh

WATTE classification (hydroMEM will also incorporate with WATTE)
- Input_Water_Class = "40" #Don't use 0 (no data)
- Input_Other_Class= "55" land
- Input_Marsh_Class = "16,23,32" #Classification(s), Marsh, low, medium, high productivity

**Jin’s priority list**

1. Nanobind and parallelization


**Jin’s request order**


## Pete's development
## Soil Cohorts ##
## Excluded **pyadcircmodule** dependency (src_branch folder: This folder will be deleted/merged later)

## Development, in no particular order


