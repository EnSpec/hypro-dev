#!/usr/bin/env python3
#
# Copyright (C) 2021 by Nanfeng Liu <nliu58@wisc.edu>
# Licensed GNU GPL v3; see `../LICENSE` for complete terms.
#
# Environmental Spectroscopy Laboratory
# Dept. of Forest and Wildlife Ecology
# University of Wisconsin - Madison

""" Functions to do atmospheric corrections. """

import os
import logging

import numpy as np

logger = logging.getLogger(__name__)


def atm_corr_band(atm_lut_WVC, atm_lut_VIS, atm_lut_VZA, atm_lut_RAA, atm_lut,
                  wvc_image, vis_image, vza_image, raa_image, rdn_image,
                  bg_mask):
    """ Do atmospheric correction for one band.

    Parameters
    ----------
    atm_lut_WVC, atm_lut_VIS, atm_lut_VZA, atm_lut_RAA: list of floats
        Atmospheric lookup table water vapor column, visibility, view zenith and relative azimuth angle grids.
    atm_lut: ndarray
        Atmospheric lookup table, shape = (RHO, WVC, VIS, VZA, RAA).
    wvc_image, vis_image, vza_image, raa_image: 2D arrays.
        Water vapor column, visibility, view zenith and relative azimuth angle images.
    rdn_image: 2D array
        Radiance image.
    bg_mask: 2D bool array
        Background mask.

    Returns
    -------
    rho: 2D array
        Surface reflectance.
    """

    from scipy.interpolate import RegularGridInterpolator

    # Interpolate the lookup table
    idx = ~bg_mask
    pts = np.array([wvc_image[idx], vis_image[idx], vza_image[idx], raa_image[idx]]).T
    
    interp_params = (atm_lut_WVC, atm_lut_VIS, atm_lut_VZA, atm_lut_RAA)
    
    def alt_interpolator(alt):
        return RegularGridInterpolator(interp_params, alt)
    
    # Calculate interpolated radiance at albedo [0.0, 0.5, 1.0]
    rdn_000 = alt_interpolator(atm_lut[0])(pts)
    rdn_050 = alt_interpolator(atm_lut[1])(pts)-rdn_000
    rdn_100 = alt_interpolator(atm_lut[2])(pts)-rdn_000

    del alt_interpolator, pts

    # Do atmospheric corrections
    S = (rdn_100-2*rdn_050)/(rdn_100-rdn_050+1e-10)
    F = rdn_100*(1-S)
    A = rdn_image[idx]-rdn_000
    rho = np.zeros(rdn_image.shape)
    rho[idx] = A/(F+S*A)

    # Clear data
    del idx
    del A, S, F
    del rdn_000, rdn_050, rdn_100

    return rho


def atm_corr_image(flight_dict):
    """ Whole-image atmospheric correction.

    Parameters
    ----------
    flight_dict: dict
        Flight dictionary.
    """

    if os.path.exists(flight_dict['refl_file']):
        logger.info('Write the reflectance image to %s.' %flight_dict['refl_file'])
        return

    from ENVI    import read_envi_header, write_envi_header
    from AtmLUT  import read_binary_metadata

    # Read radiance image header
    rdn_header = read_envi_header(os.path.splitext(flight_dict['merged_rdn_file'])[0]+'.hdr')

    # Read atmospheric lookup table
    atm_lut_metadata = read_binary_metadata(flight_dict['resampled_atm_lut_file']+'.meta')
    atm_lut_metadata['shape'] = tuple([int(v) for v in atm_lut_metadata['shape']])
    
    atm_lut_WVC = np.array([float(v) for v in atm_lut_metadata['WVC']])
    atm_lut_VIS = np.array([float(v) for v in atm_lut_metadata['VIS']])
    atm_lut_VZA = np.array([float(v) for v in atm_lut_metadata['VZA']])
    atm_lut_RAA = np.array([float(v) for v in atm_lut_metadata['RAA']])

    atm_lut = np.memmap(flight_dict['resampled_atm_lut_file'],
                        dtype=atm_lut_metadata['dtype'],
                        mode='r',
                        shape=atm_lut_metadata['shape']) # shape=(RHO, WVC, VIS, VZA, RAA, WAVE)

    # Read scan angles (SCA) image
    sca_header = read_envi_header(os.path.splitext(flight_dict['merged_sca_file'])[0]+'.hdr')
    saa = float(sca_header['sun azimuth'])
    sca_image = np.memmap(flight_dict['merged_sca_file'],
                          dtype='float32',
                          shape=(sca_header['bands'],
                                 sca_header['lines'],
                                 sca_header['samples']))
    # View zenith angle (VZA)
    vza_image = np.copy(sca_image[0,:,:])
    # Relative azimuth angle (RAA)
    raa_image = saa-sca_image[1,:,:]
    raa_image[raa_image<0] += 360.0
    raa_image[raa_image>180] = 360.0-raa_image[raa_image>180]
    
    # Clear data
    sca_image.flush()
    del sca_header, saa
    del sca_image

    # Read WVC & VIS images
    wvc_header = read_envi_header(os.path.splitext(flight_dict['wvc_file'])[0]+'.hdr')
    tmp_wvc_image = np.memmap(flight_dict['wvc_file'],
                          mode='r',
                          dtype='float32',
                          shape=(wvc_header['lines'],
                                 wvc_header['samples']))
    wvc_image = np.copy(tmp_wvc_image)

    vis_header = read_envi_header(os.path.splitext(flight_dict['vis_file'])[0]+'.hdr')
    tmp_vis_image = np.memmap(flight_dict['vis_file'],
                              dtype='float32',
                              mode='r',
                              shape=(vis_header['lines'],
                                     vis_header['samples']))
    vis_image = np.copy(tmp_vis_image)
    tmp_wvc_image.flush()
    tmp_vis_image.flush()
    del wvc_header, vis_header
    del tmp_vis_image, tmp_wvc_image

    # Read background mask
    bg_header = read_envi_header(os.path.splitext(flight_dict['background_mask_file'])[0]+'.hdr')
    bg_mask = np.memmap(flight_dict['background_mask_file'],
                        dtype='bool',
                        mode='r',
                        shape=(bg_header['lines'],
                               bg_header['samples']))

    idx = ~bg_mask
    
    max_WVC = atm_lut_WVC.max()
    max_VIS = atm_lut_VIS.max()
    max_VZA = atm_lut_VZA.max()
    max_RAA = atm_lut_RAA.max()
    
    # Enforce cutoffs in ALT parameters
    wvc_image[wvc_image>=max_WVC] = max_WVC-0.1
    vis_image[vis_image>=max_VIS] = max_VIS-0.1
    vza_image[vza_image>=max_VZA] = max_VZA-0.1
    raa_image[raa_image>=max_RAA] = max_RAA-0.1

    del max_WVC, max_VIS, max_VZA, max_RAA
    
    # Remove outliers in WVC and VIS
    wvc = wvc_image[idx]
    avg_wvc = wvc.mean()
    std_wvc = wvc.std()
    wvc_image[(np.abs(wvc_image-avg_wvc) > 2*std_wvc)&(~bg_mask)] = avg_wvc
    del wvc
    
    vis = vis_image[idx]
    avg_vis = vis.mean()
    std_vis = vis.std()
    vis_image[(np.abs(vis_image-avg_vis)>2*std_vis)&(~bg_mask)] = avg_vis
    del vis
    
    del idx

    fid = open(flight_dict['refl_file'], 'wb')
    # Do atmospheric correction
    info = 'Bands = '
    for band in range(rdn_header['bands']):
        if band%20==0:
            info += '%d, ' %(band+1)
        if (rdn_header['wavelength'][band]>=1340.0 and rdn_header['wavelength'][band]<=1440.0) or (rdn_header['wavelength'][band]>=1800.0 and rdn_header['wavelength'][band]<=1980.0) or rdn_header['wavelength'][band]>=2460.0:
            fid.write(np.zeros((rdn_header['lines'], rdn_header['samples'])).astype('float32').tostring())
        else:
            offset = rdn_header['header offset']+4*band*rdn_header['lines']*rdn_header['samples']# in bytes
            rdn_image = np.memmap(flight_dict['merged_rdn_file'],
                                  dtype='float32',
                                  mode='r',
                                  offset=offset,
                                  shape=(rdn_header['lines'],
                                         rdn_header['samples']))
            refl = atm_corr_band(atm_lut_WVC, atm_lut_VIS, atm_lut_VZA, atm_lut_RAA, np.copy(atm_lut[...,band]),
                                 wvc_image, vis_image, vza_image, raa_image, rdn_image,
                                 bg_mask)
            fid.write(refl.astype('float32').tostring())
            rdn_image.flush()
            del refl, rdn_image

    fid.close()
    info += '%d, Done!' %band
    logger.info(info)

    # Clear data
    del wvc_image, vis_image, vza_image, raa_image
    atm_lut.flush()
    bg_mask.flush()
    del atm_lut, bg_mask

    rdn_header['description'] = 'Reflectance [0-1]'
    write_envi_header(os.path.splitext(flight_dict['refl_file'])[0]+'.hdr', rdn_header)
    logger.info('Write the reflectance image to %s.' %flight_dict['refl_file'])
