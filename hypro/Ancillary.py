#!/usr/bin/env python3
#
# Copyright (C) 2021
# Brendan Heberlein <bheberlein@wisc.edu>
# Nanfeng Liu <nliu58@wisc.edu>
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison
#
# HyPro is free software, made available
# for use & modification according to the
# GNU General Public License 3.0 (GPL v3).
# See `../LICENSE` for complete terms.

import os
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

from ENVI import read_envi_header, empty_envi_header, write_to_envi


def make_ancillary_header(shape, srs, gt, *, descr='', nodata=-9999.0):
    
    hdr = empty_envi_header()
    hdr['file type'] = 'ENVI Standard'
    hdr['description'] = descr
    
    hdr['samples'] = shape[1]
    hdr['lines'] = shape[0]
    try:
        hdr['bands'] = shape[2]
    except IndexError:
        hdr['bands'] = 1
    
    hdr['byte order'] = 0
    hdr['header offset'] = 0
    hdr['interleave'] = 'bsq'
    hdr['data type'] = 4
    
    hdr['data ignore value'] = nodata
    
    hdr['coordinate system string'] = srs.ExportToWkt()
    
    hdr['map info'] = [
        srs.GetAttrValue('projcs').replace(',', ''),
        1, 1, gt[0], gt[3], gt[1], gt[1], ' ',' ',
        srs.GetAttrValue('datum').replace(',', ''),
        srs.GetAttrValue('unit')
    ]
    
    
    if srs.GetAttrValue('PROJECTION').lower() == 'transverse_mercator':
        hdr['map info'][7] = srs.GetUTMZone()
        if gt[3]>0.0:
            hdr['map info'][8] = 'North'
        else:
            hdr['map info'][8] = 'South'
    
    return hdr


def calculate_ancillary_products(flight_dict, nodata=-9999.0):
    
    msk_fn = Path(flight_dict['background_mask_file'])
    dem_fn = Path(flight_dict['merged_dem_file'])
    sca_fn = Path(flight_dict['merged_sca_file'])
    
    # Get flightline basename
    basename = '_'.join(os.path.basename(msk_fn).split('_')[:2])
    
    # Read the background mask
    msk_ds = gdal.Open(str(msk_fn))
    msk = msk_ds.GetRasterBand(1).ReadAsArray().astype(bool)
    
    # While the dataset is open, grab the geotransform array
    gt = msk_ds.GetGeoTransform()
    
    msk_ds = None
    
    # Calculate DEM aspect
    aspect_option = gdal.DEMProcessingOptions(trigonometric=False, computeEdges=True)
    aspect_ds = gdal.DEMProcessing(f'/vsimem/_tmp_{basename}_Aspect.tif',
                                   str(dem_fn), 'aspect', options=aspect_option)
    
    # Calculate DEM slope
    slope_option = gdal.DEMProcessingOptions(computeEdges=True)
    slope_ds = gdal.DEMProcessing(f'/vsimem/_tmp_{basename}_Slope.tif',
                                  str(dem_fn), 'slope', options=slope_option)

    # Convert slope & aspect to radians
    dem_aspect = np.radians(aspect_ds.GetRasterBand(1).ReadAsArray())
    dem_slope = np.radians(slope_ds.GetRasterBand(1).ReadAsArray())

    # Close temporary files
    aspect_ds = None
    slope_ds = None

    # Load SCA image header
    sca_hdr = read_envi_header(sca_fn.with_suffix('.hdr'))

    # Read SCA image
    sca_shape = (sca_hdr['bands'], sca_hdr['lines'], sca_hdr['samples'])
    sca_image = np.memmap(sca_fn, dtype='float32', mode='r', shape=sca_shape)
    
    # Get sensor viewing angles from SCA image
    view_zenith =  np.radians(sca_image[0])
    view_azimuth = np.radians(sca_image[1])
    
    # Calculate solar angles
    solar_zenith = np.radians(float(sca_hdr['sun zenith']))
    solar_azimuth = np.radians(float(sca_hdr['sun azimuth']))
    
    relative_azimuth = dem_aspect-solar_azimuth
    
    # Generate flat SZA image
    sza_image = np.full_like(view_zenith, nodata)
    sza_image[~msk] = solar_zenith

    # Generate flat SAA image
    saa_image = np.full_like(view_azimuth, nodata)
    saa_image[~msk] = solar_azimuth
    
    # Cosine of I
    cos_i = (np.cos(solar_zenith)*np.cos(dem_slope) + 
             np.sin(solar_zenith)*np.sin(dem_slope) * np.cos(relative_azimuth))

    # Phase
    # Wanner et al. JGRA 1995, eq. 51
    # Schlapfer et al. IEEE-TGARS 2015, eq. 2
    cos_phase = np.cos(solar_zenith)*np.cos(view_zenith) + np.sin(solar_zenith)*np.sin(view_zenith)*np.cos(relative_azimuth)
    phase = np.arccos(cos_phase)
    
    # Set no data values
    dem_aspect[msk] = nodata
    dem_slope[msk] = nodata
    view_zenith[msk] = nodata
    view_azimuth[msk] = nodata
    relative_azimuth[msk] = nodata
    cos_i[msk] = nodata

    # Get SRS
    srs = osr.SpatialReference()
    srs.ImportFromWkt(sca_hdr['coordinate system string'])
    
    anc_dir = Path(flight_dict['merge_dir'])/'ancillary'
    anc_dir.mkdir(exist_ok=True)
    
	# Write DEM slope
    slope_hdr = make_ancillary_header(dem_slope.shape, srs, gt, descr='DEM slope, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_DSMSlope', dem_slope, slope_hdr)
    
    # Write DSM aspect
    aspect_hdr = make_ancillary_header(dem_aspect.shape, srs, gt, descr='DEM aspect, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_DSMAspect', dem_aspect, aspect_hdr)
    
    # Write sensor zenith angles
    vza_hdr = make_ancillary_header(view_zenith.shape, srs, gt, descr='Sensor viewing zenith angle, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_VZA', view_zenith, vza_hdr)
    
    # Write sensor azimuth angles
    vaa_hdr = make_ancillary_header(view_azimuth.shape, srs, gt, descr='Sensor viewing azimuth angle, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_VAA', view_azimuth, vaa_hdr)
    
    # Write solar zenith angle
    sza_hdr = make_ancillary_header(sza_image.shape, srs, gt, descr='Solar zenith angle, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_SZA', sza_image, sza_hdr)
    
    # Write solar azimuth angle
    saa_hdr = make_ancillary_header(saa_image.shape, srs, gt, descr='Solar azimuth angle, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_SAA', saa_image, saa_hdr)
	
	# Write relative azimuth 
    raa_hdr = make_ancillary_header(relative_azimuth.shape, srs, gt, descr='Relative viewing azimuth angle, in [rad]', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_RAA', relative_azimuth, raa_hdr)
    
	# Write cosine I
    cosi_hdr = make_ancillary_header(cos_i.shape, srs, gt, descr='Cosine of I', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_CosI', cos_i, cosi_hdr)
    
	# Write phase
    phase_hdr = make_ancillary_header(phase.shape, srs, gt, descr='Phase', nodata=-9999.0)
    write_to_envi(anc_dir/f'{basename}_Phase', phase, phase_hdr)