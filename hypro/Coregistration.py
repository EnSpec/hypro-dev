import os
import warnings

from pathlib import Path

import numpy as np
from osgeo import gdal

from geoarray import GeoArray
from arosics.CoReg_local import COREG_LOCAL
from arosics.DeShifter import DESHIFTER

import logging

logger = logging.getLogger(__name__)

def coerce_to_geoarray(dataset):
    """ Attempt to coerce inputs of various types to a GeoArray object. """
    
    if isinstance(dataset, GeoArray):
        geo_array = dataset
    
    elif isinstance(dataset, gdal.Dataset):
        gt = dataset.GetGeoTransform()
        crs = dataset.GetProjection()
        crs_wkt = crs.ExportToWkt()
        geo_array = GeoArray(dataset, geotransform=gt, projection=crs_wkt)
    
    elif isinstance(dataset, np.memmap):
        ### Need to pass gt & crs_wkt as arguments...
        raise NotImplementedError
    
    # elif isinstance(dataset, hypro.spatial.RasterMap):
    #     raise NotImplementedError
    
    elif isinstance(dataset, Path):
        geo_array = GeoArray(str(dataset))
    
    elif type(dataset) is str:
        geo_array = GeoArray(dataset)
    
    else:
        raise TypeError(f'Cannot coerce type "{type(dataset)}" to GeoArray.')
    
    return geo_array

_default_kwargs_CRL = {
    'grid_res': 50, # Tie point spacing in x-direction, in terms of target image pixels
    'window_size': (164, 164), # Size of matching window
    'nodata': (-9999.0,)*2, # No-data values for reference & target images
    'max_iter': 8,
    'r_b4match': 1, # Reference band for matching
    's_b4match': 1, # Target band for matching
    'CPUs': 1,
    'max_shift': 50,
}

def coregister(shifted_image, reference_image, **kwargs):
    """ Get `COREG_LOCAL` instance to align two images by warping to local shifts. """
    
    ### What does it look like when coregistration fails?
    
    # Ignore annoying warnings about lossy PROJ strings (AROSICS should fix)
    warnings.filterwarnings('ignore', category=UserWarning, module='pyproj')
    
    # Load input images as `GeoArray` instances
    shifted_image = coerce_to_geoarray(shifted_image)
    reference_image = coerce_to_geoarray(reference_image)
    
    # Update default coregistration arguments with keyword arguments
    coreg_kwargs = {**_default_kwargs_CRL, **kwargs}
    
    # Do coregistration
    CRL = COREG_LOCAL(reference_image, shifted_image, **coreg_kwargs)
    # result = CRL.correct_shifts()
    
    # return result
    return CRL


def do_coregistration(flight_dict, config):
    """ Do sensor-to-sensor coregistration for multisensor acquisitions.
    
    NOTES
    -----
    Expects the following information to be present in the configuration file:
    
    {
      "do_coregistration": true,
      "reference_sensor": "vnir",
      "reference_band": 185,
      "target_bands": {"swir": 8}
    }
    
    """
    
    if len(flight_dict['sensors']) > 1:
        
        coreg_dict = config.get('Coregistration', {})
        
        if coreg_dict.get('do_coregistration', False) is True:
            
            reference_sensor = coreg_dict.get('reference_sensor', [*flight_dict['sensors'].keys()].pop())
            reference_image = flight_dict['sensors'][reference_sensor]['ortho_rdn_image_file']
            reference_band = coreg_dict.get('reference_band', 1)
            
            logger.info(f'Co-register radiance images to reference sensor: {reference_sensor}.')
            logger.info(f'Using reference band: {reference_band}')
            
            for sensor_index, sensor_dict in flight_dict['sensors'].items():
                
                # Don't co-register the reference image
                if sensor_index == reference_sensor:
                    continue
                
                # Get file names of input & output images
                shifted_image = sensor_dict['ortho_rdn_image_file']
                deshifted_image = shifted_image+'CoReg'
                
                if os.path.exists(deshifted_image):
                    logger.info('Write the coregistered image to %s.' %deshifted_image)
                    continue
                
                logger.info('Co-registering sensor: %s' %sensor_index)
                
                # Get matching band for target image
                target_band = coreg_dict['target_bands'].get(sensor_index, 1)
                logger.info('Using matching band: %s' %target_band)
                
                logger.info('Preparing for co-registration...')
                # Initialize coregistration result
                CRL = coregister(shifted_image, reference_image,
                                    r_b4match=reference_band,
                                    s_b4match=target_band,
                                    grid_res=100)
                logger.info('Calculating local shifts...')
                # Calculate shifts
                result = CRL.correct_shifts()
                logger.info('Co-registration complete.')
                
                logger.info('Write the co-registered image to %s.' %deshifted_image)
                # Save deshifted image
                result['GeoArray_shifted'].save(deshifted_image)
                # Update sensor dictionary
                sensor_dict['ortho_rdn_image_file'] = deshifted_image
                
                logger.info('Apply co-registration shifts to orthorectified DEM.')
                # Apply shift correction to ortho DEM
                ortho_dem = sensor_dict['ortho_dem_image_file']
                deshifted_dem = ortho_dem+'CoReg'
                dem_result = DESHIFTER(ortho_dem, CRL.coreg_info).correct_shifts()
                # Save & update sensor dictionary
                dem_result['GeoArray_shifted'].save(deshifted_dem)
                logger.info('Write the co-registered DEM to %s.' %deshifted_dem)
                sensor_dict['ortho_dem_image_file'] = deshifted_dem
                
                logger.info('Apply co-registration shifts to orthorectified SCA.')
                # Apply shift correction to ortho SCA
                ortho_sca = sensor_dict['ortho_sca_image_file']
                deshifted_sca = ortho_sca+'CoReg'
                sca_result = DESHIFTER(ortho_sca, CRL.coreg_info).correct_shifts()
                # Save & update sensor dictionary
                sca_result['GeoArray_shifted'].save(deshifted_sca)
                logger.info('Write the co-registered SCA to %s.' %deshifted_sca)
                sensor_dict['ortho_sca_image_file'] = deshifted_sca
                
        else:
            logger.info('Image coregistration disabled; moving on.')
    else:
        logger.info('Not a multisensor dataset; moving on.')