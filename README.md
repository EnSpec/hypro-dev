<p align="center" text-align="center">
<br>
<img width="200" alt="HyPro" src="https://user-images.githubusercontent.com/18175461/113648777-cca85280-9652-11eb-9438-5b4d324cfd9b.png">
<br>
<br>
HyPro v1.0.0
<br>
© 2021 Nanfeng Liu
<br>
Licensed under GNU GPLv3
</p>

<br>

# Introduction

## Description

This package provides methods for geometric, radiometric & atmospheric corrections on hyperspectral data cubes to convert raw sensor DN values to at-surface reflectance.

The package main script is `HyspexPro.py`, which draws on methods defined throughout the module.

(1) HyspexPro aims to do geometric and radiometric corrections on HySpex images. Our imaging system consists of two HySpex cameras (VNIR-1800 and SWIR-384) and one GPS-Inertial Measurement Unit (IMU) sensor (iMAR iTraceRT-F400). The VNIR-1800 sensor has 182 spectral bands within 400-1000 nm (spectral sampling = 3.26 nm). The SWIR-384 sensor has 288 bands within 960-2500 nm (spectral sampling = 5.45 nm). For more details about these cameras, please refer to the [NEO website](https://www.hyspex.no/products/). The iTraceRT-F400 sensor records the attitudes (roll, pitch and heading angles) and positions (longitude, latitude and altitude) of the imaging platform.

(2) Geometric correction focuses on removing image geometric distortions mainly caused by platform motions, which can be characterized by roll, pitch and heading angles. We adapted the ray-tracing method proposed by Meyer<sup>[[1]](#Meyer1994)</sup> to calculate the map coordinates (easting, northing and elevation) of each raw image pixel (a process called geo-referencing). Then, the `warp` tool provided by [GDAL](https://gdal.org) (Geospatial Data Abstraction Library) was used to resample the raw image to regular map grids (a process called geo-rectification).

(3) Radiometric correction focuses on reducing radiometric distortions mainly caused by sensor smile effects (shifts in sensor center wavelengths and full-widths-at-half-maximum (FWHMs)) and atmospheric effects (absorptions and scatterings). The spectral matching method proposed by Gao<sup>[[2]](#Gao2004)</sup> was adapted to detect sensor smile effects. Then, raw radiance spectra were resampled to common center wavelengths using a cubic interpolation. In the atmospheric correction step, two atmosphere parameters (water vapor column and visibility) were retrieved from image spectra. The water vapor column was estimated via the Continuum Interpolated Band Ratio (CIBR) technique proposed by Kaufman.<sup>[[3]](#Kaufman1992)</sup> The estimation of visibility was based on the Dense Dark Vegetation (DDV) method proposed by Kaufman.<sup>[[4]](#Kaufman1997)</sup> All radiometric correction steps require an atmosphere lookup table which consists of the at-senor radiance simulated under different atmospheric and sun-target-view geometric conditions. An open-source radiative transfer model, LibRadTran, was used for this purpose. For more details about LibRadTran, please refer to the [website](http://www.libradtran.org/doku.php).

<br>

Below is the standard workflow of HyspexPro:

<br>

<p align="center" text-align="center">
<img width="800" alt="DN to Reflectance Pipeline" src="https://user-images.githubusercontent.com/18175461/113648788-d16d0680-9652-11eb-8e62-5cfdadf43e5b.png">
</p>

<br>

(4) The input dataset of HyspexPro includes:
  * HySpex raw digital number (DN) images;
  * HySpex IMU & GPS data extracted from HySpexNAV;
  * HySpex sensor model files provided by NEO;
  * Surface model (DSM, DEM, DTM, etc.);
  * HySpex lab radiometric re-calibration data provided by NEO (optional);
  * Ground control points data (optional).

<br>

(5) The output dataset of HyspexPro includes:
  * HySpex ground surface reflectance images;
  * View angles (zenith and azimuth);
  * Surface model.

<br>


## Permissions & Disclaimer

HyPro is Python package that provides tools for processing raw hyperspectral
image data from pushbroom sensors.

HyPro v1.0.0 is authored by Nanfeng Liu (nliu58@wisc.edu), with contributions from Adam Chlus (chlus@wisc.edu) & Brendan Heberlein (bheberlein@wisc.edu). HyPro v1.0.0 is licensed for public use under [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html). For details of the licensing agreement, refer to [`LICENSE`](https://github.com/EnSpec/hypro/blob/main/LICENSE). In brief:

```text
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

## References

<a name="Meyer1994"><sup><b>[1]</b></sup></a> Meyer P. (1994). Remote Sens. Environ. 49(2): 118-130. _doi:[10.1016/0034-4257(94)90048-5](https://doi.org/10.1016/0034-4257(94)90048-5)_

<a name="Gao2004"><sup><b>[2]</b></sup></a> Gao B., et al. (2004). Remote Sens. Environ. 90(4): 424-433. _doi:[10.1016/j.rse.2003.09.002](https://doi.org/10.1016/j.rse.2003.09.002)_

<a name="Kaufman1992"><sup><b>[3]</b></sup></a> Kaufman Y., et al. (1992). IEEE Trans. Geosci. Remote Sens. 30(5): 871–884. _doi:[10.1109/36.175321](https://doi.org/10.1109/36.175321)_

<a name="Kaufman1997"><sup><b>[4]</b></sup></a> Kaufman Y., et al. (1997). IEEE Trans. Geosci. Remote Sens. 35(5): 1286-1298. _doi:[10.1109/36.628795](https://doi.org/10.1109/36.628795)_
