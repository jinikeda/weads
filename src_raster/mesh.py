

# ----------------------------------------------------------
# F U N C T I O N   F I L E E X I S T S
# ----------------------------------------------------------
# Determine if a file exists
# result = function(file)
# ----------------------------------------------------------
def isInMesh(x, y, myMesh):
    if not os.path.exists(file):
        print(file + ' NOT FOUND.\tPROGRAM EXIT.')
        exit()
# ----------------------------------------------------------

# ----------------------------------------------------------
# F U N C T I O N    C O M P U T E C E N T R O I D
# ----------------------------------------------------------
#
# Computes the centroid of an ADCIRC element
# result = function(elementID,myMesh,meshConn,transformer)
# ----------------------------------------------------------


def computeCentroid(elementID, myMesh, meshConn, transformer):
    # First, find the centroid of the element
    x1 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][0])).x()
    x2 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][1])).x()
    x3 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][2])).x()
    y1 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][0])).y()
    y2 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][1])).y()
    y3 = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][2])).y()
    xC, yC = transformer.transform(
        ((y1 + y2 + y3) / 3.0),
        ((x1 + x2 + x3) / 3.0))
    return xC, yC
# ----------------------------------------------------------

# ----------------------------------------------------------
# F U N C T I O N    E L E M E N T B O U N D I N G B O X
# ----------------------------------------------------------
#
# Computes the bounding box of an ADCIRC element
# result = function(elementID,myMesh,meshConn,transformer,convert)
# ----------------------------------------------------------


def elementBoundingBox(elementID, myMesh, meshConn, transformer, convert):
    # Compute the bounding box of element
    # i/Users/mvb49270/Documents/repo/pyHydroMEM/src/interpolationMethods.py
    elemX = [None] * 3
    elemY = [None] * 3
    elemX[0] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][0])).x()
    elemX[1] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][1])).x()
    elemX[2] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][2])).x()
    elemY[0] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][0])).y()
    elemY[1] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][1])).y()
    elemY[2] = myMesh.node(myMesh.nodeIndexById(meshConn[elementID][2])).y()
    if convert == 0:
        elem_xMin = min(elemX)
        elem_yMin = min(elemY)
        elem_xMax = max(elemX)
        elem_yMax = max(elemY)
    elif convert == 1:
        elem_xMin, elem_yMin = transformer.transform(
            min(elemY), min(elemX))
        elem_xMax, elem_yMax = transformer.transform(
            max(elemY), max(elemX))

    return elem_xMin, elem_yMin, elem_xMax, elem_yMax, elemX, elemY

# ----------------------------------------------------------

# ----------------------------------------------------------
# F U N C T I O N    F I N D R A S T E R S T E N C I L
# ----------------------------------------------------------
#
# Computes the row/col stencil for a given element
# Find the raster cell of element i's bbox
# (xMin, yMin & xMax, yMax)
# result = function(bbox, element x-y min/max, gridSize)
# ----------------------------------------------------------


def findRasterStencil(bbox, numRows, numCols,
                      elem_yMin, elem_yMax,
                      elem_xMin, elem_xMax, gridSize):

    minRow = math.ceil((bbox[3] - elem_yMin) / gridSize)
    maxRow = math.ceil((bbox[3] - elem_yMax) / gridSize)
    minCol = math.ceil((elem_xMin - bbox[0]) / gridSize)
    maxCol = math.ceil((elem_xMax - bbox[0]) / gridSize)
    if minRow < 0:
        minRow = 0
    if maxRow > numRows:
        maxRow = numRows
    if minCol < 0:
        minCol = 0
    if maxCol > numCols:
        maxCol = numCols

    return minRow, maxRow, minCol, maxCol

# ----------------------------------------------------------

# ----------------------------------------------------------
# F U N C T I O N    C O M P U T E E L E M E N T A R E A
# ----------------------------------------------------------
#
# Compute the element area using linear basis functions
# result = function(bbox, element x-y min/max, gridSize)
# ----------------------------------------------------------


def computeElementArea(elemX, elemY, transformer):
    x1, y1 = transformer.transform(elemY[0], elemX[0])
    x2, y2 = transformer.transform(elemY[1], elemX[1])
    x3, y3 = transformer.transform(elemY[2], elemX[2])

    a1 = x3 - x2
    a2 = x1 - x3
    a3 = x2 - x1
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2

    elementArea = (b1 * a2 - b2 * a1) / 2.0

    return elementArea, a1, a2, a3, b1, b2, b3,\
        x1, x2, x3, y1, y2, y3
# ----------------------------------------------------------
