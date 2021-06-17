#!/usr/bin/env python3
#
# Copyright (C) 2021 Brendan Heberlein <bheberlein@wisc.edu>
# Licensed GNU GPL v3; see `../LICENSE` for complete terms.
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison

""" Spatial reference systems. """

import osgeo
from osgeo import osr

def set_axis_mapping(crs):
    # Ensure (x, y) coordinate ordering
    if int(osgeo.__version__[0]) >= 3:
        # By default, GDAL 3 respects the axis ordering specified by the SRS
        # See: https://github.com/OSGeo/gdal/issues/1546
        crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)


def get_crs(**kwargs):

    crs = osr.SpatialReference()
    key, val = [*kwargs.items()].pop()
    key = key.lower()

    if key == 'epsg':
        crs.ImportFromEPSG(val)
    elif key == 'esri':
        crs.ImportFromESRI(val)
    elif key == 'usgs':
        crs.ImportFromUSGS(val)
    elif key == 'proj4':
        crs.ImportFromProj4(val)
    elif key == 'wkt':
        crs.ImportFromWkt(val)
    else:
        raise KeyError(f'Invalid key "{key}".')

    set_axis_mapping(crs)

    return crs


def utm_crs(zone, north=True, datum='WGS84'):
    crs = osr.SpatialReference()
    crs.SetWellKnownGeogCS(datum)
    crs.SetUTM(zone, north)
    return crs

# Define coordinate systems
WGS84 = get_crs(epsg=4326)
WGS84_UTM15N = get_crs(epsg=32615)
WGS84_UTM16N = get_crs(epsg=32616)