# Development Note (Aug 28 version)

---

## Running file
- hydromem.py (Argument file for running WAED)
  
## Procedure files (src folders)
- grd2dem.py (Convert ADCIRC unstructured mesh file to tiff file)
- hyconn.py (Calculate hydro connectivity and classify ocean, land, and pond)
- tidaldatums.py (Calculate tidal datums in each pixel on water region **this part is computationally expensive)
- tidaldatumsidw.py (Interpolate the tidal datum on land regions)
- mem.py (Calculate vegetation productivities)
- rast2adc.py (renew a mesh file: fort.14 based on mem simulation)
- update_nodal_attributes.py (renew nodal attribute file: fort.13 based on mem simulation)
  - manning.py (Specify Manning's roughnesses based on biomass productivities but probably doesn't use this file)
  
## Supplemental files 
- raster.py (rasterization functions)
- basics.py (check existing a file)


# To-do-list

## Development (mandatory)
**src/mem.py**
- Read current wetland marsh/mangrove distributions (Jin).
  - Jin made a distribution mapping using read_NWI_file.py, which still needs to be incorporated into the package sequence (16 cores take 30 mins).<br> input file: chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO/inputs/Wetlands_NWI_forMeshRefinement/Region3_NWI_LC_Reclassify_wetlandsOnly.tif,<br> output file: chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO/outputs/Domain_classification_distribution_resample100.tif
  - Jin added a multi-species option in hydromem.py. src/mem.py needs to be modified (ongoing)
  
- Examine and calculate each species' (ecological) response (Pete).

## Development (desirable)
- Replace private module: **pyadcircmodules**
- Read maxele.63 and create a max inundation map (Pete creates a Python file). Jin modifies hydromem.py (ongoing)
~~- MEM 5 classifications in src/mem.py (Jin done on Aug 22)~~
- Pete modifies src/tidaldatums.py (avoid double for loops) -> Jin and Chris will work on Cython with parallelization.

## Modifications (desirable)
~~- src/hyconn.py (currently used for loops and time-consuming + not sure about pond classification; Jin slightly modified the code Aug 19)~~

## Modifications (optional)
- src/tidaldatumsidw.py: GDAL IDW may not perform well (Jin and Linoj used another approach in CRMS2MAP, which needs to be considered, pending Jin Aug 22)

# NOTE
## Input/Output data storage 
chenier.cct.lsu.edu:/data/CCR_data/ACTIVE/ESLR2021TCB/WEAD/IO
---

NWI classification
- 8 = salt marsh (regularly flooded) -> follow with tidal cycle
- 9 = mangrove
- 20 = irregularly flooded marsh

WATTE classification (hydroMEM also incorporate with WATTE)
- Input_Water_Class = "40" #Don't use 0 (no data)
- Input_Other_Class= "55" land
- Input_Marsh_Class = "16,23,32" #Classification(s), Marsh, String

**Jin’s priority list**

1. Read current wetland marsh/mangrove distributions in src/mem.py
2. Inundation level part (hydromem.py)
3. WATTE modification + Evaluate productivity and Inundation level map
4. Cython and parallelization
5. Check fort.13

**Jin’s request order**

1. maxinundepth.63
2. Inundationdepth.tif (maxinundationdepth.63 -> tiff file)


## Pete's note for ecology.py (this code will be deleted later)
## Development, in no particular order

- Adjustment of Manning’s n for open-water conversion – TBD
- NWI classification is supportive for model initialization; however, how to evolve into the future with multi-type distribution? (Jin will start considering)

## Variables, listed in Pythonic order
0. idx, local (local index)
1. idx, global (global index)
2. x, geo (longitude)
3. y, geo (latitude)
4. x, cpp (easting)
5. y, cpp (northing)
6. z (elevation)
7. n (manning)
8. ed (everdried)
9. it (inundation time)
10. mi (maximum inundation depth)
11. mlw (mean low water)
12. msl (mean sea level)
13. mhw (mean high water)
14. mlw, interp (mean low water, interpolated)
15. mhw, interp (mean high water, interpolated)
16. nwi (integer-valued National Wetlands Inventory)
17. B (biomass)
18. A (accretion)
19. z-adj (elevation, adjusted)
20. n-adj (Manning, adjusted)

## Modules, listed in no particular order
- numpy, math, sys, tqdm, time

## Scripts, listed in order of execution with respective runtime (mm:ss) for TCBa
1. preprocessing.py (ADCIRC or another model) - 04:26
2. hydrodynamics.py (agnostic—points) - 09:56
3. interpdatums.py (agnostic—points) - 35:42
4. ecology.py(d) (agnostic—points) - 00:19
5. postprocessing.py (ADCIRC or another model) - XX:XX
   Total runtime: 50 minutes, 23 seconds + XX:XX (postprocessing)

a TCB: 3,527,549 mesh nodes; 8 attributes; 23 tidal constituents; 3 static global outputs
b Bounding box set externally via lower-left corner (lon/lat) and upper-right corner (lon/lat)
c Three wetland types: salt marsh (regularly flooded); mangrove; irregularly flooded marsh
d Accounts for multi-type calculation of biomass and accretion with updating of elevation and Manning’s n

## Files, listed according to type
### Input (6, derived from present simulation)
- fort.14
- fort.NWI.14
- fort.13
- fort.53
- everdried.63
- inundationtime.63
- maxinundepth.63 -> Need this file output (Jin needs to confirm.)

### Input (1, assigns bounding box for ROI—region of interest)
- ROI.pts

### Intermediate (6)
- filteredData.pts (clipped to bounding box)
- filteredHarmAmp.pts (clipped to bounding box)
- filteredHarmPha.pts (clipped to bounding box)
- filteredTidalDatums.pts (clipped to bounding box)
- filteredTidalDatumsIDW.pts (clipped to bounding box, intertidal points only)
- filteredBiomassAccretion.pts (clipped to bounding box, intertidal points only)

### Output (2, configured for subsequent simulation)
- fort.dt.14
- fort.dt.13


