
# Development Note (Aug 28, 24 version)

-  **Discussed the development direction on Aug 28th**
---

## Input/Output data storage 
- Folder Path: chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO -> publically available path (need to update)

## Running file
- WEADS_Raster.py or WEADS_Point.py (Argument file to run WEAD)

# To-do-list

## Development (mandatory)
**src_point will be the main source code and eliminate the dependency of the src_raster approach**
  - rasterize (eliminate singularity image file)
**src_raster/mem.py**
  - The only problem with reading vegetation classification is that a segmentation fault happens when reading a large tiff file with a docker.<br>input file: Folder Path/inputs/Wetlands_NWI_forMeshRefinement/Region3_NWI_LC_Reclassify_wetlandsOnly.tif,<br> output file: Folder Path/outputs/Domain_classification_distribution_resample100.tif NWI.py is a currently not elegant.
  
## Development (desirable)
- Compare Morris's Excel VBA sheet
- Modify average_horizontal_eddy_viscosity_in_sea_water_wrt_depth based on update fort.14
- Consider the climate and catastrophic aspects of WEADS development (Shabnam starts climate analysis and literature review).
- Pete and Shabnam will work on soil cohorts using the NUMAR model approach.

## Modifications (desirable)

## Modifications (optional)
- Read Netcdf outputs

**Jin’s priority list**

1. Excluded **pyadcircmodule** dependency (src_raster folder: This folder will be deleted/merged later)
2. Nanobind and parallelization

## Main Python module/function files (src raster folders)
- tidaldatums.py (Calculate tidal datums in each pixel on ocean region. **This part is still computationally expensive on Sep 12**)
  

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
