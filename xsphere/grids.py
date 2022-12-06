#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:55:53 2021

@author: ghiggi
"""
import numpy as np
import matplotlib.patches as mpatches
from scipy.spatial import SphericalVoronoi
 
#-----------------------------------------------------------------------------.
# ##############################
#### Coordinates conversion ####
# ##############################
# radius = 6371.0e6
def lonlat2xyz(longitude, latitude, radius=1):
    """From 2D geographic coordinates to cartesian geocentric coordinates."""
    lon, lat = np.deg2rad(longitude), np.deg2rad(latitude)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return x, y, z


def xyz2lonlat(x,y,z, radius=1):
    """From cartesian geocentric coordinates to 2D geographic coordinates."""
    latitude = np.rad2deg(np.arcsin(z / radius))  
    longitude = np.rad2deg(np.arctan2(y, x)) 
    return longitude, latitude 

## Testing 
# x, y, z = xsphere.lonlat2xyz(lon, lat)
# lon1, lat1 = xsphere.xyz2lonlat(x,y,z)   
# np.testing.assert_allclose(lon, lon1) 
# np.testing.assert_allclose(lat, lat1) 

##----------------------------------------------------------------------------.
## Conversion to spherical coordinate is buggy 
# def xyz2sph(x,y,z):
#     """From cartesian geocentric coordinates to spherical polar coordinates."""
#     r = np.sqrt(x**2 + y**2 + z**2)  
#     theta = np.arccos(z/r) 
#     phi = np.arctan(y, x)
#     return theta, phi, r

# def sph2xyz(theta, phi, radius=1):
#     """From spherical polar coordinates to cartesian geocentric coordinates."""
#     x = radius * np.sin(theta) * np.cos(phi)
#     y = radius * np.sin(theta) * np.sin(phi)
#     z = radius * np.cos(theta)
#     return x, y, z

# def lonlat2sph(longitude, latitude, radius=1):
#     """From 2D geographic coordinates to spherical polar coordinates."""
#     x, y, z = lonlat2xyz(longitude=longitude, latitude=latitude, radius=radius)
#     return xyz2sph(x,y,z)
    
# def sph2lonlat(theta, phi, radius=1):
#     """From spherical polar coordinates to 2D geographic coordinates."""
#     x, y, z = sph2xyz(theta=theta, phi=phi, radius=1)
#     return xyz2lonlat(x,y,z)

## Testing 
# x, y, z = xsphere.lonlat2xyz(lon, lat)
# theta, phi, r = xsphere.xyz2sph(x,y,z)
# x1, y1, z1 = xsphere.sph2xyz(theta, phi, r)
# np.testing.assert_allclose(x, x1)
# np.testing.assert_allclose(y, y1)
# np.testing.assert_allclose(z, z1)

#-----------------------------------------------------------------------------.


def get_polygons_2D_coords(lon_bnds, lat_bnds):
    """Create a list of numpy [x y] array polygon vertex coordinates from CDO lon_bnds and lat_bnds matrices.""" 
    # Input: n_polygons x n_vertices 
    # Output: list (for each polygon) of numpy_array [x, y] polygon vertex coordinates
    n_polygons = lon_bnds.shape[0]
    n_vertices = lon_bnds.shape[1]
    list_polygons_xy = list()
    for i in range(n_polygons):
        poly_corners = np.zeros((n_vertices, 2), np.float64)
        poly_corners[:,0] = lon_bnds[i,:]
        poly_corners[:,1] = lat_bnds[i,:]
        list_polygons_xy.append(poly_corners)
    return(list_polygons_xy)


def get_PolygonPatchesList_from_latlon_bnds(lon_bnds, lat_bnds, fill=True):
    """Create a list of Polygon mpatches from CDO lon_bnds and lat_bnds."""
    # Construct list of polygons 
    l_polygons_xy = get_polygons_2D_coords(lon_bnds=lon_bnds, lat_bnds=lat_bnds)
    l_Polygon_patch = [mpatches.Polygon(xy=p, closed=False, fill=fill) for p in l_polygons_xy]
    return l_Polygon_patch   


def get_PolygonPatchesList(l_polygons_xy, fill=True):
    """Create Polygon mpatches from a numpy [x y] array with polygon vertex coordinates."""
    # Construct list of mpatches.Polygon
    l_Polygon_patch = [mpatches.Polygon(xy=p, closed=False, fill=fill) for p in l_polygons_xy]
    return l_Polygon_patch


def SphericalVoronoiMesh(lon, lat):
    """
    Infer the mesh of a spherical sampling from the mesh node centers provided in 2D geographic coordinates.
    
    Parameters
    ----------
    lon : numpy.ndarray
        Array of longitude coordinates (in degree units).
    lat : numpy.ndarray
        Array of latitude coordinates (in degree units).

    Returns
    -------
    list_polygons_lonlat : list
        List of numpy.ndarray with the polygon mesh vertices for each graph node.

    """
    # Convert to geocentric coordinates 
    radius = 6371.0e6 # radius = 1 can also be used
    x, y, z = lonlat2xyz(lon, lat, radius=radius)
    coords = np.column_stack((x,y,z))
    # Apply Spherical Voronoi tesselation
    sv = SphericalVoronoi(coords,
                          radius=radius, 
                          center=[0, 0, 0])
    ##-------------------------------------------------------------------------.
    # SphericalVoronoi object methods 
    # - sv.vertices : Vertex coords
    # - sv.regions : Vertex ID of each polygon
    # - sv.sort_vertices_of_regions() : sort indices of vertices to be clockwise ordered
    # - sv.calculate_areas() : compute the area of the spherical polygons 
    ##-------------------------------------------------------------------------.
    # Sort vertices indexes to be clockwise ordered
    sv.sort_vertices_of_regions()
    # Retrieve area 
    area = sv.calculate_areas() 
    ##-------------------------------------------------------------------------.
    # Retrieve list of polygons coordinates arrays
    list_polygons_lonlat = []
    for region in sv.regions:  
        tmp_xyz = sv.vertices[region]    
        tmp_lon, tmp_lat = xyz2lonlat(tmp_xyz[:,0],tmp_xyz[:,1],tmp_xyz[:,2], radius=radius)    
        list_polygons_lonlat.append(np.column_stack((tmp_lon, tmp_lat)))
    ##-------------------------------------------------------------------------.
    return list_polygons_lonlat, area
    
 



