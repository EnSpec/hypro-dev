#!/usr/bin/env python3
#
# Copyright (C) 2021 Brendan Heberlein <bheberlein@wisc.edu>
# Licensed GNU GPL v3; see `../LICENSE` for complete terms.
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison

""" ... """

from pathlib import Path

from osgeo import ogr, osr
from shapely import geometry, wkb

# import numpy as np

# from hypro.spatial.crs import WGS84
from spatial.crs import WGS84

GEOM_MAP = {
    geometry.Point: ogr.wkbPoint,
    geometry.LineString: ogr.wkbLineString,
    geometry.Polygon: ogr.wkbPolygon,
    geometry.LinearRing: ogr.wkbLinearRing
}

def common_geometry_type(geometries):
    """ Find a common geometry type among a list of `shapely` geometries. """

    types = set(type(g) for g in geometries)

    if len(types) > 1:
        raise TypeError('All geometries must share a common type!')

    return types.pop()

def geometry_to_shp(geoms, fn, layername=None, crs=WGS84):
    """ Write `shapely` geometries to an ESRI shapefile.
    
    By default, geometries are assumed to be referenced to WGS84.
    
    NOTES
    -----
    - `srs` must be the SRS of all input geometries; no
      reprojection will occur if the SRS is wrong!

    TODO
    ----
    - Allow for writing of attribute fields
    """

    fn = Path(fn)

    # Initialize the shapefile datasource
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(str(fn))

    # Set layer name
    if not layername:
        layername = fn.stem
    # Create layer
    layer = geometry_to_layer(geoms, crs, ds, name=layername)

    # Write projection to file
    with fn.with_suffix('.prj').open(mode='w') as f:
        f.write(crs.ExportToWkt())

    layer = None
    ds = None


def geometry_from_shp(shp, to_srs=None):
    """ Load an ESRI shapefile to `shapely` geometries.
    
    Each feature is treated as a single, separate geometry.
    Optionally, geometries will be reprojected into a
    given target SRS, if different from the native SRS.
    
    PARAMETERS
    ----------
    shp : str or Path
        Filepath pointing to the input shapefile.
    
    to_srs : osr.SpatialReference
        Desired target SRS.
    
    RETURNS
    -------
    list
        Sequence of `shapely` geometries.
    
    TODO
    ----
    - Allow for retrieval of feature attribute, e.g. by
      optionally applying a user-specified function to
      each feature which would control the output type
    """
    
    # Path objects are coerced to sring types
    if isinstance(shp, Path): shp = str(shp)

    handle = ogr.Open(shp)
    layer = handle.GetLayer(0)
    
    transform = None
    
    # If a target SRS is given...
    if to_srs is not None:
        from_srs = layer.GetSpatialRef()
        # ...and it's different from the native SRS;
        if not from_srs.IsSame(to_srs):
            # Define a coordinate transformation to be applied
            transform = osr.CoordinateTransformation(from_srs, to_srs)
    
    geoms = [geometry_from_feature(f, transform=transform) for f in layer]
    
    # Close the shapefile
    handle = None
    
    return geoms


def geometry_from_feature(feature, transform=None):
    """ Load an OGR feature to a `shapely` geometry.
    
    An optional coordinate transformation can be applied.
    """
    # Get the geometry
    geom_ogr = feature.GetGeometryRef()
    
    # Transform the geometry, if needed
    if transform: geom_ogr.Transform(transform)
    
    # Convert to `shapely` geometry via WKB representation
    geom_wkb = geom_ogr.ExportToWkb()
    geom = wkb.loads(geom_wkb)
    
    return geom


def geometry_to_feature(geom, layer_defn, transform=None):
    """ Coerce a `shapely` geometry to an OGR feature.
    
    An optional coordinate transformation can be applied.
    """
    
    # Convert to OGR geometry via WKB representation
    geom_ogr = ogr.CreateGeometryFromWkb(geom.wkb)
    
    # Transform the geometry, if needed
    if transform: geom_ogr.Transform(transform)
    
    # Create a new feature
    feature = ogr.Feature(layer_defn)
    feature.SetGeometry(geom_ogr)
    
    return feature


def geometry_to_layer(geoms, crs, dataset, name=''):
    """ Create an OGR layer containing the input `shapely` geometries. """
    
    geometry_type = GEOM_MAP[common_geometry_type(geoms)]
    
    # Create layer
    layer = dataset.CreateLayer(name, srs=crs, geom_type=geometry_type)
    
    ### TODO: Add attribute fields here
    
    # Get layer definition
    layer_defn = layer.GetLayerDefn()
    
    # Add each geometry as a separate feature
    for geom in geoms:
        assert GEOM_MAP[type(geom)] is geometry_type
        feat = geometry_to_feature(geom, layer_defn, transform=None)
        layer.CreateFeature(feat)
    
    dataset.FlushCache()
        
    return layer