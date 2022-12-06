#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 12:54:57 2022

@author: ghiggi
"""

# TODO:
# - add_HealpixMesh,
# - add_EquiangularMesh,
# - add_GaussianLegendreMesh,
# - add_CubedMesh
# --> pgysp graph?

# - reshape equiangular to lat, lon core dim
# - reshape_equiangular_to_unstructured

# - Check_mesh()  --> Check Polygon mpatches

# - 3D plots

## Add from shapefile
# - da.sphere.add_mesh_from_shp(fpath)  # poly.shp
# - da.sphere.add_nodes_from_shp(fpath) # point.shp
# - da.sphere.save_mesh_to_shp()
# - da.sphere.save_nodes_to_shp()

## Spherical polygons area computations (now planar assumption)
# - https://github.com/anutkk/sphericalgeometry
# - https://stackoverflow.com/questions/4681737/how-to-calculate-the-area-of-a-polygon-on-the-earths-surface-using-python
