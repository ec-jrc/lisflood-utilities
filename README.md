# Lisflood Utilities

This repository hosts source code of LISFLOOD utilities.
Go to [Lisflood OS page](https://ec-jrc.github.io/lisflood/) for more information.

Other useful resources

| **Project**         | **Documentation**                                         | **Source code**                                                 |
| ------------------- | --------------------------------------------------------- | --------------------------------------------------------------- |
| Lisflood            | [Model docs](https://ec-jrc.github.io/lisflood-model/)    | https://github.com/ec-jrc/lisflood-code                         |
|                     | [User guide](https://ec-jrc.github.io/lisflood-code/)     |                                                                 |
| Lisvap              | [Docs](https://ec-jrc.github.io/lisflood-lisvap/)         | https://github.com/ec-jrc/lisflood-lisvap                       |
| Calibration tool    | [Docs](https://ec-jrc.github.io/lisflood-calibration/)    | https://github.com/ec-jrc/lisflood-calibration                  |
| Lisflood Utilities  |                                                           | https://github.com/ec-jrc/lisflood-utilities (this repository)  |
| Lisflood Usecases   |                                                           | https://github.com/ec-jrc/lisflood-usecases                     |


## Intro and documentation

Lisflood Utilities is a set of tools to help LISFLOOD users (or any users of PCRaster/NetCDF files)
to execute some mundane tasks that are necessary to operate LISFLOOD.

The documentation about the available tools and how to use them can be found in the [Wiki](https://github.com/ec-jrc/lisflood-utilities/wiki) of this repository.

## Installation

### Requisites

The easy way is to use conda environment as they incapsulate C dependencies as well, so you wouldn't need to install libraries.

Otherwise, ensure you have properly installed the following software:

- Python 3.5+
- GDAL C library and software
- NetCDF4 C library

### Install

If you use conda, create a new env (or use an existing one) and install gdal and lisflood-utilities:

```bash
conda create --name myenv python=3.7 -c conda-forge
conda activate myenv
conda install -c conda-forge pcraster eccodes gdal
pip install lisflood-utilities
```

If you don't use conda but a straight Python3 virtualenv:

```bash
source /path/myenv/bin/activate
pip install lisflood-utilities
```

If GDAL library fails to install, ensure to install the same package version of the C library you have on your system. You may also need to setup paths to GDAL headers.

To check which version of GDAL libraries you have installed on your computer, use `gdal-config`

```bash
sudo apt-get install libgdal-dev libgdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
gdal-config --version  # 3.0.1
pip install GDAL==3.0.1
```

> **Note:** if you previously installed an older version of `lisflood-utilities`, it is highly recommended to remove it before installing the newest version:

```bash
pip uninstall lisflood-utilities
pip install -e./
```

## Using `lisfloodutilities` programmatically 

You can use lisflood utilities in your Python programs. As an example, the script below creates the mask map for a set of stations (_stations.txt_). The mask map is a boolean map with 1 and 0. 1 is used for all (and only) the pixels hydrologically connected to one of the stations. The resulting mask map is in PCRaster format.

```python
from lisfloodutilities.cutmaps.cutlib import mask_from_ldd
from lisfloodutilities.nc2pcr import convert
from lisfloodutilities.readers import PCRasterMap

ldd = 'tests/data/cutmaps/ldd_eu.nc'
clonemap = 'tests/data/cutmaps/area_eu.map'
stations = 'tests/data/cutmaps/stations.txt'

ldd_pcr = convert(ldd, clonemap, 'tests/data/cutmaps/ldd_eu_test.map', is_ldd=True)[0]
mask, outlets_nc, maskmap_nc = mask_from_ldd(ldd_pcr, stations)
mask_map = PCRasterMap(mask)
print(mask_map.data)
```