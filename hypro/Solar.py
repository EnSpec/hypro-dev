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

import numpy as np

def get_sun_earth_distance(when):
    """ Get sun-earth distance of a day.

    Parameters
    ----------
    when: datetime object
        Date and time.

    Returns
    -------
    d: float
        Sun-Earth distance.
    """

    import numpy as np
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    sun_earh_distance_file = os.path.abspath(os.path.join(cur_dir, os.pardir, 'data/solar-data/sun_earth_distance.dat'))
    d = np.loadtxt(sun_earh_distance_file)
    doy = when.timetuple().tm_yday
    d = d[doy-1]

    return d

def get_sun_angles(longitude, latitude, utc_time):
    """ Calculate the Sun's position.

    References
    ----------
    (1) Manuel Blanco-Muriel, et al. (2001). Computing the solar vector. Solar Energy, 70(5), 431-441.
    (2) The C code is available from: http://www.psa.es/sdg/sunpos.htm

    Parameters
    ----------
    longitude: float
        Longitude, in degrees. West: negative, East: positive.
    latitude: float
        Latitude, in degrees. Northern: positive, Southern: negative.
    utc_time: datetime object
        UTC time.

    Returns
    -------
    ZenithAngle: float
        Sun zenith angle, in degrees.
    AzimuthAngle: float
        Sun azimuth angle, in degrees.
    """

    rad = np.pi/180
    EarthMeanRadius = 6371.01 # in [km]
    AstronomicalUnit = 149597890 # in [km]

    DecimalHours = utc_time.hour+(utc_time.minute+utc_time.second/60.0)/60.0
    Aux1 = int((utc_time.month-14)/12)
    Aux2 = int(1461*(utc_time.year+4800+Aux1)/4 +
                 367*(utc_time.month-2-12*Aux1)/12 -
                 3*(utc_time.year + 4900+ Aux1)/100/4 +
                 utc_time.day-32075)
    JulianDate = Aux2-0.5+DecimalHours/24.0
    ElapsedJulianDays = JulianDate-2451545.0

    Omega = 2.1429-0.0010394594*ElapsedJulianDays
    MeanLongitude = 4.8950630+ 0.017202791698*ElapsedJulianDays
    MeanAnomaly = 6.2400600+ 0.0172019699*ElapsedJulianDays
    EclipticLongitude = MeanLongitude + 0.03341607*np.sin(MeanAnomaly)+ 0.00034894*np.sin(2*MeanAnomaly)-0.0001134-0.0000203*np.sin(Omega)
    EclipticObliquity = 0.4090928-6.2140e-9*ElapsedJulianDays+0.0000396*np.cos(Omega)
    Sin_EclipticLongitude = np.sin(EclipticLongitude)
    Y = np.cos(EclipticObliquity)*Sin_EclipticLongitude
    X = np.cos(EclipticLongitude)
    RightAscension = np.arctan2(Y, X)
    if RightAscension < 0.0:
        RightAscension = RightAscension + np.pi*2
    Declination = np.arcsin(np.sin(EclipticObliquity)*Sin_EclipticLongitude)

    GreenwichMeanSiderealTime = 6.6974243242+0.0657098283*ElapsedJulianDays+DecimalHours
    LocalMeanSiderealTime = (GreenwichMeanSiderealTime*15+longitude)*rad
    HourAngle = LocalMeanSiderealTime-RightAscension
    LatitudeInRadians = latitude*rad
    Cos_Latitude = np.cos(LatitudeInRadians)
    Sin_Latitude = np.sin(LatitudeInRadians)
    Cos_HourAngle= np.cos(HourAngle)
    ZenithAngle = np.arccos(Cos_Latitude*Cos_HourAngle*np.cos(Declination)+np.sin(Declination)*Sin_Latitude)
    Y = -np.sin(HourAngle)
    X = np.tan(Declination)*Cos_Latitude-Sin_Latitude*Cos_HourAngle
    Azimuth = np.arctan2(Y,X)
    if Azimuth < 0.0:
        Azimuth = Azimuth + np.pi*2
    Parallax = (EarthMeanRadius/AstronomicalUnit)*np.sin(ZenithAngle)
    ZenithAngle = (ZenithAngle+Parallax)/rad
    AzimuthAngle = Azimuth/rad

    return [ZenithAngle, AzimuthAngle]