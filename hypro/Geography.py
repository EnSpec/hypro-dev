#!/usr/bin/env python3
#
# Copyright (C) 2021 by Nanfeng Liu <nliu58@wisc.edu>
# Licensed GNU GPL v3; see `../LICENSE` for complete terms.
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison

""" Functions to process map coordinate systems. """

import os

import numpy as np
import osgeo
from osgeo import gdal, osr


def set_crs_mapping(crs):
    """ Set CRS axis mapping to (x,y) format. """
    if int(osgeo.__version__[0]) >= 3:
        # By default, GDAL 3 respects axis ordering specified by the CRS
        crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)


def get_utm_zone(lon):
    """ Calculate UTM zone.

    Parameters
    ----------
    lon: float
        Longitude, in degrees. West: negative, East: positive.

    Returns
    -------
    zone: int
        UTM zone number.
    """

    zone = int(1+(lon+180.0)/6.0)

    return zone


def is_northern(lat):
    """ Determine if it is northern hemisphere.

    Parameters
    ----------
    lat: float
        Latitude, in degrees. Northern: positive, Southern: negative.

    Returns
    -------
    int
        1: northern, 0: southern.
    """

    if lat < 0.0:
        return 0
    else:
        return 1


def define_utm_crs(lon, lat):
    """ Define a UTM map coordinate system.

    Parameters
    ----------
    lon: float
        Longitude, in degrees. West: negative, East: positive.
    lat: float
        Latitude, in degrees. Northern: positive, Southern: negative.

    Returns
    -------
    crs: osr.SpatialReference object
        UTM map coordinate system.
    """

    crs = osr.SpatialReference()
    crs.SetWellKnownGeogCS("WGS84")
    zone = get_utm_zone(lon)
    crs.SetUTM(int(zone), int(is_northern(lat)))
    
    set_crs_mapping(crs)

    return crs


def define_wgs84_crs():
    """ Define a WGS84 map coordinate system.

    Returns
    -------
    crs: osr.SpatialReference object
        WGS84 map coordinate system.
    """

    crs = osr.SpatialReference()
    crs.SetWellKnownGeogCS("WGS84")
    
    set_crs_mapping(crs)

    return crs


def get_raster_crs(file):
    """ Get the map coordinate system of a raster image.

    Parameters
    ----------
    file: str
        Georeferenced image filename.

    Returns
    -------
    crs: osr object
        Map coordinate system.
    """

    ds = gdal.Open(file, gdal.GA_ReadOnly)
    prj = ds.GetProjection()
    ds = None

    crs = osr.SpatialReference()
    crs.ImportFromWkt(prj)
    
    set_crs_mapping(crs)

    return crs


def get_grid_convergence(lon, lat, map_crs):
    """ Get grid convergence angles.

    Parameters
    ----------
    lon: list of floats
        Longitude. West: negative; East: positive.
    lat: list of floats
        Latitde. North: positive; South: negative.

    Reuturns
    --------
    grid_convergence: array of floats
        Grid convergence in degrees.
    """

    lon, lat = np.array(lon), np.array(lat)
    if map_crs.GetAttrValue('PROJECTION').lower() == 'transverse_mercator':
        lon0 = map_crs.GetProjParm('central_meridian')
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
        lon0 = np.deg2rad(lon0)
        grid_convergence = np.arctan(np.tan(lon-lon0)*np.sin(lat))
        grid_convergence = np.rad2deg(grid_convergence)
    else:
        delta_lat = 1e-4
        lon_lat_0 = np.array([lon, lat]).T
        lon_lat_1 = np.array([lon, lat+delta_lat]).T
        wgs84_crs = define_wgs84_crs()
        transform = osr.CoordinateTransformation(wgs84_crs, map_crs)
        xy0 = np.array(transform.TransformPoints(lon_lat_0))
        xy1 = np.array(transform.TransformPoints(lon_lat_1))
        dx = xy1[:,0]-xy0[:,0]
        dy = xy1[:,1]-xy0[:,1]
        grid_convergence = np.abs(np.rad2deg(np.arcsin(dx/np.sqrt(dx**2+dy**2))))
        index = dx*dy>0
        grid_convergence[index] = -grid_convergence[index]

    return grid_convergence


def get_map_crs(dem, longitude, latitude):
    """ Get map coordinate system.

    Notes
    -----
    If `dem` is a file, the map coordinate system should be
    the same as that of the dem file; otherwise define a UTM coordinate system
    based on the longitude and latitude.

    Parameters
    ----------
    dem: str or float
        DEM image filename, or user-specified DEM value.
    longitude, latitude: float
        Longitude and latitude.

    Returns
    -------
    map_crs: osr object
        Map coordinate system.
    """

    if os.path.isfile(dem):
        map_crs = get_raster_crs(dem)
    else:
        map_crs = define_utm_crs(longitude, latitude)

    return map_crs
