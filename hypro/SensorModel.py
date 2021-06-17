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

""" Functions to create a sensor model. """

import os
import logging

import numpy as np

logger = logging.getLogger(__name__)


def make_sensor_model(sensor_model_file, fov, ifov, samples, if_rotated):
    """ Generate a sensor model.

    Parameters
    ----------
    sensor_model_file: str
        The sensor model filename.
    fov: float
        Sensor fov [deg].
    ifov: float
        Sensor instaneous fov [mrad].
    samples: int
        Image columns.
    if_rotated: bool
        If the sensor is 180 degree rotated.
    """

    if os.path.exists(sensor_model_file):
        logger.info('Write the sensor model to %s.' %sensor_model_file)
        return

    sensor_model = np.zeros((samples, 3))
    sensor_model[:,0] = np.arange(samples)

    fov = np.deg2rad(fov)
    if if_rotated:
        xs = np.linspace(np.tan(fov/2), -np.tan(fov/2), num=samples)
    else:
        xs = np.linspace(-np.tan(fov/2), np.tan(fov/2), num=samples)

    sensor_model[:,1] = np.arctan(xs)
    sensor_model[:,2] = ifov/1000
    np.savetxt(sensor_model_file,
               sensor_model,
               header='pixel    vinkelx    vinkely',
               fmt='%d    %.10f    %.10f')
    del sensor_model, xs

    logger.info('Write the sensor model to %s.' %sensor_model_file)


def determine_if_rotated(imu_gps_file):
    """ Determine if the sensor is 180 degree rotated.

    Parameters
    ----------
    imu_gps_file: str
        IMUGPS file

    Returns
    -------
    True or False: bool
        Whether the sensor is 180 degree rotated.
    """

    imugps = np.loadtxt(imu_gps_file)
    flag = np.sign((imugps[-1,3] - imugps[0,3])*imugps[0,6])

    if flag == -1:
        logger.info('The sensor is 180 degree rotated.')
        return True
    else:
        logger.info('The sensor is not 180 degree rotated.')
        return False
