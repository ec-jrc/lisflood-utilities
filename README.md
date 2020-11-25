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


## Intro

Lisflood Utilities is a set of tools to help LISFLOOD users (or any users of PCRaster/netCDF files)
to execute some mundane tasks that are necessary to operate lisflood.
Here's a list of utilities you can find in lisflood-utilities package.

* __pcr2nc__ is a tool to convert PCRaster maps to netCDF4 files.
  - convert a single map into a NetCDF4 file
  - convert a time series of maps into a NetCDF4 mapstack
  - support for WGS84 and ETRS89 (LAEA) reference systems
  - fine tuning of output files (compression, significant digits, etc.)
 
* __nc2pcr__ is a tool to convert a netCDF file into PCRaster maps.
  - convert 2D variables in single PCRaster maps
  - netCDF4 mapstacks are not supported yet

* __cutmaps__ is a tool to cut netcdf files in order to reduce size, using either
  - a bounding box of coordinates
  - a bounding box of matrix indices
  - an existing boolean area mask
  - a list of stations and a LDD (in netCDF or PCRaster format) **Note: PCRaster must be installed in the conda env**
 
* __compare__ is a package containing a set of simple Python classes that helps to compare 
netCDF, PCRaster and TSS files.

The package contains convenient classes for reading/writing:

* PCRasterMap
* PCRasterReader
* NetCDFMap
* NetCDFWriter

### Installation

#### Requisites
The easy way is to use conda environment as they incapsulate C dependencies as well, so you wouldn't need to install libraries.

Otherwise, ensure you have properly installed the following software:

- Python 3.5+
- GDAL C library and software
- netCDF4 C library

#### Install
If you use conda, create a new env (or use an existing one) and install gdal and lisflood-utilities:

```bash
conda create --name myenv python=3.7 -c conda-forge
conda activate myenv
conda install -c conda-forge pcraster gdal
pip install lisflood-utilities
```

If you don't use conda but a straight python3 virtualenv:

```bash
source /path/myenv/bin/activate
pip install lisflood-utilities
```

If GDAL library fails to install, ensure to install the same package version of the
C library you have on your system. You may also need to setup paths to gdal headers.

To check which version of GDAL libraries you have installed on your computer, use gdal-config

```bash
sudo apt-get install libgdal-dev libgdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
gdal-config --version  # 3.0.1
pip install GDAL==3.0.1
```

## pcr2nc

### Usage

> __Note:__ This guide assumes you have installed the program with pip tool.
> If you cloned the source code instead, just substitute the executable `pcr2nc` with `python pcr2nc_script.py` that is in the root folder of the cloned project.

The tool takes three command line input arguments:

- -i, --input: It can be a path to a single file, a folder or a unix-like widlcard expression like _/path/to/files/dis00*_
- -o, --output_file: Path to the output nc file
- -m, --metadata: Path to a yaml or json file containing configuration for the NetCDF4 output file.

Unless the input is a single file, the resulting NetCDF4 file will be a mapstack according to a time dimension.

Input as a folder containing PCRaster maps. In this case, the folder must contain ONLY PCraster files and the output will be a mapstack.

```bash
pcr2nc -i /path/to/input/ -o /path/to/output/out.nc -m ./nc_metadata.yaml
```

Input as a path to a single map. In this case, the output won't be a mapstack.

```bash
pcr2nc -i /path/to/input/pcr.map -o /path/to/output/out.nc -m ./nc_metadata.yaml
```

Input as a _Unix style pathname pattern expansion_. The output will be a mapstack. __Note that in this case the input argument must be contained in double quotes!__

```bash
pcr2nc -i "/path/to/input/pcr00*" -o /path/to/output/out.nc -m ./nc_metadata.json
```

#### Writing metadata configuration file

Format of resulting NetCDF4 file is configured into a metadata configuration file. This file can be written in YAML or JSON format.

An example of a metadata configuration file is the following

```yaml
variable:
  shortname: dis
  description: Discharge
  longname: discharge
  units: m3/s
  compression: 9
  least_significant_digit: 2
source: JRC Space, Security, Migration
reference: JRC Space, Security, Migration
geographical:
  datum: WGS84
time:
  calendar: proleptic_gregorian
  units: days since 1996-01-01

```

#### Variable section

In `variable` section you can configure metadata for the main variable:

- `shortname`: A short name for the variable
- `longname`: The long name version
- `description`: A description for humans
- `units`: The units of the variable
- `compression`: Optional, integer number between 1 and 9, default 0 (no compression). If present the output nc file will be compressed at this level.
- `least_significant_digit`: Optional, integer number, default 2. From NetCDF4 documentation:

> If your data only has a certain number of digits of precision
(say for example, it is temperature data that was measured with a precision
of 0.1 degrees), you can dramatically improve zlib compression by quantizing
(or truncating) the data using the least_significant_digit keyword argument
to createVariable. The least significant digit is the power of ten of the
smallest decimal place in the data that is a reliable value.
For example if the data has a precision of 0.1,
then setting least_significant_digit=1 will cause data the data to be
quantized using `numpy.around(scale*data)/scale`, where `scale = 2**bits`,
and bits is determined so that a precision of 0.1 is retained
(in this case bits=4). Effectively, this makes the compression 'lossy'
instead of 'lossless', that is some precision in the data is sacrificed for the sake of disk space.

#### Source and reference

`source` and `reference` add information for the institution that is providing the NetCDF4 file.

#### Geographical section

In `geographical` section the only setting to configure is `datum`. 
Currently, pcr2nc supports the following list:

  * `WGS84`
  * `ETRS89`
  * `GISCO`

#### Time section

This section is optional and is only required if the output file is a mapstack (a timeseries of georeferenced 2D arrays)
In this section you have to configure `units` and `calendar`.

- `units`: Can be one of the following strings (replacing placeholders with the actual date):
    - `hours since YYYY-MM-DD HH:MM:SS`
    - `days since YYYY-MM-DD`
- `calendar`: A recognized calendar identifier, like `proleptic_gregorian`, `gregorian` etc.

## nc2pcr

This tool converts single maps netCDF (time dimension is not supported yet) to PCRaster format.

### Usage

```bash
nc2pcr -i /path/to/input/file.nc -o /path/to/output/out.map -c /path/to/clone.map
```

If input file is a LDD map, you must add the `-l` flag:

```bash
nc2pcr -i /path/to/input/ldd.nc -o /path/to/output/ldd.map -c /path/to/clone.map -l
```

## Cutmaps: a NetCDF files cookie-cutter

This tool cut netcdf files, using a mask, a bounding box or a list of stations along with a LDD map.  

### Usage:
The tool accepts as input:

* a mask map (either PCRaster or netCDF format) or 
  - alternatively, matrix indices in the form xmini_xmaxi:ymini_ymaxi or
  - alternatively, coordinates bounding box in the form xmin_xmax:ymin_ymax
  - alternatively, list of stations with coordinates and a LDD map.
* a path to a folder containing netCDF files to cut or a static dataset path like LISFLOOD static files. 
* a path to a folder where to write cut files.

The following command will cut all netcdf files inside _/workarea/Madeira/lai/_ folder 
and produced files will be writte in current folder. 
The cookie-cutter that will be used is _/workarea/Madeira/maps/MaskMap/Bacia_madeira.nc_. 
This file is a mask (boolean map with 1 only in the area of interest) where cutmaps takes the bounding box from.
The mask can also be in PCRaster format.

```bash
cutmaps -m /workarea/Madeira/maps/MaskMap/Bacia_madeira.nc -f /workarea/Madeira/lai/ -o ./
```

**Indices can also be passed as an argument (using -c argument instead of -m). Knowing your area of interest from your netCDF files, 
you can determine indices of the array and you can pass in the form `imin_i_max:jmin_jmax`.**

```bash
cutmaps -c 150_350:80_180 -f /workarea/Madeira/lai/ -o ./
```

**Example with coordinates and path to EFAS/GloFAS static data (-S option), with -W to allow overwriting existing files in output directory:**

```bash
cutmaps -S /home/projects/lisflood-eu -c 4078546.12_4463723.85:811206.57_1587655.50 -o /Work/Tunisia/cutmaps -W
```

**Example with stations.txt and LDD**

Given a LDD map and a list of stations in a text file, each row having coordinates X/Y or lat/lon and an index, separated by tabs:

```text
4297500	1572500 1
4292500	1557500 2
4237500	1537500 3
4312500	1482500 4
4187500	1492500 5
```

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.map -N stations.txt -o /Work/Tunisia/cutmaps
```

If ldd is in netCDF format, LDD will be converted to PCRaster format, first.

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.nc -N stations.txt -o /Work/Tunisia/cutmaps
``` 

If you experience problems, you can try to pass a path to a PCRaster clone map.

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.nc -C area.map -N stations.txt -o /Work/Tunisia/cutmaps
```
You will find the produced mask.map and mask.nc for your area in the same folder of ldd map; you will need it for lisflood/lisvap executions.
You will also have outlets.map/outlets.nc based on stations.txt, which let you produce gauges TSS if configured in LISFLOOD.

## compare utility

This tool let you compare two netcdf datasets. You can configure it with tolerances (atol, rtol, thresholds for percentage of tolerated different values).
You can also set the option to write diff files, so that you can inspect maps and differences with a tool like Panoply

```text
usage: compare [-h] -a DATASET_A -b DATASET_B -m MASKAREA [-M SUBMASKAREA]
               [-e] [-s] [-D] [-r RTOL] [-t ATOL] [-p MAX_DIFF_PERCENTAGE]
               [-l MAX_LARGEDIFF_PERCENTAGE]

Compare netCDF outputs: 0.12.12

optional arguments:
  -h, --help            show this help message and exit
  -a DATASET_A, --dataset_a DATASET_A
                        path to dataset version A
  -b DATASET_B, --dataset_b DATASET_B
                        path to dataset version B
  -m MASKAREA, --maskarea MASKAREA
                        path to mask
  -e, --array-equal     flag to compare files to be identical
  -s, --skip-missing    flag to skip missing files in comparison
  -D, --save-diffs      flag to save diffs in netcdf files for visual
                        comparisons. Files are saved in ./diffs folder of
                        current directory.For each file presenting
                        differences, you will find files diffs, original A and
                        B (only for timesteps where differences are found).
  -r RTOL, --rtol RTOL  rtol
  -t ATOL, --atol ATOL  atol
  -p MAX_DIFF_PERCENTAGE, --max-diff-percentage MAX_DIFF_PERCENTAGE
                        threshold for diffs percentage
  -l MAX_LARGEDIFF_PERCENTAGE, --max-largediff-percentage MAX_LARGEDIFF_PERCENTAGE
                        threshold for large diffs percentage
```

## gfit

### Introduction

The tool Gfit2 was created to extract flood warning thresholds from a map stack of daily discharge for several years, given in NetCDF format.
It is designed to work with any other variable as well as with different sampling (sub- or super- daily)

The tool is made of three files:

1. Gfit2.sh: the main file, used to run the tool (e.g., ./Gfit2.sh ./settingFile_test.sh)
2. settingFile_test.sh: a setting file, where paths and values of the input variables are defined. It can be renamed, as long as the correct name is called in the script (e.g.: set1.sh --> ./Gfit2.sh ./set1.sh)
3. gfit2.r: the r script that  performs the extreme value analysis

### Requirements

You need to have installed the following:

- CDO
- R

### How to run
The tool operates in a Linux environment (also as a job with qsub). To run the tool you'll need:
- R (the script was tested with R version 3.5.0)
- the following R packages: ncdf4, lmomco, ismev
- CDO (the script was tested with CDO version 1.6.5.1)

### What does the tool do

The Gfit2 tool performs the following steps:
- It takes the input file in NetCDF and extract the series of annual maxima. an optional number of warm-up years ($warmup_yrs) are excluded at the beginning of the data, to remove potential spin-up effect. Also the map of average value is computed
- An extreme value fitting is performed on each pixel of the map, using L-moments and a 2-parameter Gumbel distribution. As option, the user can limit the analysis to a specific number of years. Also, one can use the option to remove from the fitting all values smaller than the long-term average. The fitting is performed only where at least 5 annual maxima are available, otherwise a NA is returned on the output map
- Output return level maps (corresponding to user-selected years of recurrence interval) are saved in a netcdf file (return_levels.nc) and as ascii files. If a clone map in PCRaster format is provided, maps are also saved in PCRaster (return level maps and parameters of the Gumbel distribution)
- After the EV fitting, the tool estimates some further statistics from the long-term input file, including minimum, maximum, and different percentile maps. This part can be commented out if not of interest

### Reference
[L-Moments: Analysis and Estimation of Distributions Using Linear Combinations of Order Statistics](https://www.jstor.org/stable/2345653)

Hosking, J.R.M., 1990. L-Moments: Analysis and Estimation of Distributions Using Linear Combinations of Order Statistics. J. R. Stat. Soc. Ser. B Methodol. 52, 105124.

## Using lisfloodutilities programmatically TODO

You can use lisflood utilities in your python programs.

### lisfloodutilities.readers

### lisfloodutilities.writers

### lisfloodutilities.compare

### lisfloodutilities.cutmaps

```python
from lisfloodutilities.cutmaps.cutlib import mask_from_ldd
from lisfloodutilities.nc2pcr import convert
from lisfloodutilities.readers import PCRasterMap

ldd = 'tests/data/cutmaps/ldd_eu.nc'
clonemap = 'tests/data/cutmaps/area_eu.map'
stations = 'tests/data/cutmaps/stations.txt'

ldd_pcr = convert(ldd, clonemap, 'tests/data/cutmaps/ldd_eu_test.map', is_ldd=True)[0]
mask = mask_from_ldd(ldd_pcr, stations)
mask_map = PCRasterMap(mask)
print(mask_map.data)
```

### lisfloodutilities.nc2pcr

### lisfloodutilities.pcr2nc
