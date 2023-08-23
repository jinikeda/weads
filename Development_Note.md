# Development note (Aug 22 version)

## Variables, listed in pythonic order
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
- maxinundepth.63 -> 

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
