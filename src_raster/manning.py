BL = 370.0
BU = 750.0              # Lower and upper bounds of medium productivity
nL = 0.035
nM = 0.05
nH = 0.07      # Manning's roughness based on biomass density
no = 0.022                         # Manning's roughness for open-water reclassification

band = rasterMAN.GetRasterBand(1)
man = band.ReadAsArray()
manA = man  # Initialize new Manning's field with old Manning's field

manA[B > 0.001] = nL
manA[B > BL] = nM
manA[B > BU] = nH
manA[w < 0.5] = no  # New assignments of Manning's roughness
