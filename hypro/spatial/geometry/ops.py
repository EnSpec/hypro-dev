#!/usr/bin/env python3
#
# Copyright (C) 2021 Brendan Heberlein <bheberlein@wisc.edu>
# Licensed GNU GPL v3; see `../LICENSE` for complete terms.
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison

""" Geometric operations. """

from osgeo import ogr
from shapely import wkb

def reproject(geom, transform):
    ''' Reproject a `shapely` geometry using an instance of
         `osr.CoordinateTransformation`.
        
        ARGUMENTS:
        ---------
         geom (shapely.geometry.base.BaseGeometry): The geometry to be
            reprojected; an instance of some `shapely` geometry type.
            
         transform (osr.CoordinateTransformation): A coordinate
            transformation. NOTE: the source CRS must match that
            of the input geometry.
        
        RETURNS:
        -------
         shapely.geometry.base.BaseGeometry: The reprojected geometry.
    '''
    # Convert input geometry to an `ogr.Geometry` type via WKB
    ogr_geom = ogr.CreateGeometryFromWkb(wkb.dumps(geom))
    # Do the coordinate transformation using OGR
    ogr_geom.Transform(transform)
    # Convert back to `shapely` geometry type
    reprojected = wkb.loads(ogr_geom.ExportToWkb())
    
    return reprojected