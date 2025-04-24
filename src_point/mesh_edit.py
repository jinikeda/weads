#!/usr/bin/env python
# File: mesh_edit.py
# Developer: Jin Ikeda
# Last modified: April 22, 2025

# Development Notes:
""" This script is used to preprocess the input files for the WEAD project. The input file is the ADCIRC mesh file (fort.14)
For example, thin_layer.geojson, and the script reads the input file and adjusts the z values within the domain shapefile (Thin_layer.geojson).
The script outputs the modified nodes and attributes in the model input.
For thin layer replacement, we do not expect to change the attributes significantly, so we only modified the z values.
The script is used in the WEADS project."""

##########################################################################
# --- Load internal modules ---
from . import general_functions as gfa
import pandas as pd
import geopandas as gpd
import numpy as np
import shutil

def mesh_edit(inputMeshFile, inputadjustFile, inEPSG, z_adjust=0.05):

    # --- GLOBAL PARAMETERS ---
    ndv = -99999.0
    inputEPSG = inEPSG  # Input EPSG code, e.g., 4269 for NAD83
    outputMeshFile = "fort_updated.14"
    backupMeshFile = inputMeshFile + ".bk"

    # === Load Vector File ===
    if inputadjustFile.endswith((".shp", ".geojson")):
        gdf_polygon = gpd.read_file(inputadjustFile)
    else:
        raise ValueError("Unsupported file format. Please use .shp or .geojson.")

    print("Original CRS:", gdf_polygon.crs)
    print("Number of features:", len(gdf_polygon))

    # === Convert CRS to Match Input Mesh ===
    gdf_polygon = gdf_polygon.to_crs(epsg=inEPSG)

    # === Load and Create GeoDataFrame from fort.14 ===
    ADCIRC_nodes, numNodes, numElements = gfa.read_fort14(inputMeshFile, output_Flag=False)
    gdf_ADCIRC = gfa.create_gdf(ADCIRC_nodes, ['x', 'y'], [], f"EPSG:{inEPSG}", None)

    # === Check crs ===
    assert gdf_ADCIRC.crs == gdf_polygon.crs, "CRS of ADCIRC nodes and polygon must match."

    # === Modify Elevation Within Polygon ===
    within_mask = gdf_ADCIRC.within(gdf_polygon.unary_union)
    gdf_ADCIRC.loc[within_mask, "z"] += z_adjust  # Increase Z by 0.05 meters

    print(f" Modified {within_mask.sum()} points inside the polygon(s).")

    # === Save Modified Nodes (CSV Format, Optional) ===
    print(gdf_ADCIRC.head(20))
    gdf_ADCIRC.loc[within_mask, ["node_id", "x", "y", "z"]].to_csv("ADCIRC_nodes_domain.csv", index=False)
    print(" Saved modified nodes to 'ADCIRC_nodes_domain.csv'.")

    # === Update ADCIRC Mesh ===
    shutil.copy(inputMeshFile, outputMeshFile)
    gfa.update_ADCIRC_mesh(
        outputMeshFile,
        node_id=gdf_ADCIRC["node_id"].values,
        x=gdf_ADCIRC["x"].values,
        y=gdf_ADCIRC["y"].values,
        new_z = gdf_ADCIRC["z"].values  # Inside of the function, apply negative to match ADCIRC elevation convention
    )

    print(f" Updated mesh written to '{outputMeshFile}'")

    # === Backup Original Mesh File ===
    shutil.copy(inputMeshFile, backupMeshFile)
    print(f"Backup created: {backupMeshFile}")

    # Overwrite the original file with the updated file
    shutil.move(outputMeshFile, inputMeshFile)
    print(f"{outputMeshFile} has been moved to {inputMeshFile}")



