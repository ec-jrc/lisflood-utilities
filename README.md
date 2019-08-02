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


## lisflood-utilities

User guide of each tool is placed under the relative directory in a README.md file.

### pcr2nc

__pcr2nc__ is a tool to convert PCRaster maps to netCDF4 files.
It's developed at JRC E1 directorate, as a Floods group initiative.

- convert a single map into a NetCDF4 file
- convert a time series of maps into a NetCDF4 mapstack
- support for WGS84 and ETRS89 (LAEA) reference systems
- fine tuning of output files (compression, significant digits, etc.)

#### Installation

There are two ways to install the software: one via pip tool and one by cloning the repository.

##### Requisites

Ensure you have properly installed the following software:

- Python 3.5+
- GDAL C library and software
- netCDF4 C library

Create a python3 virtualenv for using the software and activate it.

If you have virtualenvwrapper:
```bash
$ workon pcr2nc
```

Otherwise just execute the activate script
```bash
$ source /path/to/virtualenvs/pcr2nc/bin/activate
```

##### Install via pip tool

Activate the virtualenv and then:

```bash
$ pip install pcr2nc
```

After the install was complete, you still have to install the proper GDAL package,
according to the version of gdal library that is installed on your machine.
You will also need C GDAL headers to properly install python GDAL wrapper.

E.g.

```bash
sudo apt-get install libgdal-dev libgdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
gdal-config --version  # 2.2.3
pip install GDAL==2.2.3
```

##### Install by cloning the repository

Assuming your Python 3 virtualenv is called pcr2nc and you have virtualenvwrapper installed:

```bash
git clone https://github.com/ec-jrc/lisflood-utilities
cd lisflood-utilities/pcr2nc
workon pcr2nc
```

Install requirements

```bash
$ pip install -r requirements.txt
```

If GDAL library fails to install, ensure to install the same package version of the
library you have on your system.
To check which version of GDAL libraries you have installed on your computer, use gdal-config

```bash
gdal-config --version
2.1
```

Example: you have installed gdal 2.1, then:

```bash
$ pip install GDAL==2.1
$ pip install -r requirements.txt
```

#### Usage

> __Note:__ This guide assumes you have installed the program with pip tool.
> If you cloned the source code instead, just substitute the executable `pcr2nc` with `python pcr2nc_script.py` that is in the root folder of the cloned project.

The tool takes three command line input arguments:

- -i, --input: It can be a path to a single file, a folder or a unix-like widlcard expression like _/path/to/files/dis00*_
- -o, --output_file: Path to the output nc file
- -m, --metadata: Path to a yaml or json file containing configuration for the NetCDF4 output file.

Unless the input is a single file, the resulting NetCDF4 file will be a mapstack according to a time dimension.

##### Example of usages:

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

##### Writing metadata configuration file

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

##### Variable section

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

##### Source and reference

`source` and `reference` add information for the institution that is providing the NetCDF4 file.

##### Geographical section

In `geographical` section the only setting to configure is `datum`. 
Currently, pcr2nc supports the following list:

  * `WGS84`
  * `ETRS89`
  * `GISCO`

##### Time section

This section is optional and is only required if the output file is a mapstack (a timeseries of georeferenced 2D arrays)
In this section you have to configure `units` and `calendar`.

- `units`: Can be one of the following strings (replacing placeholders with the actual date):
    - `hours since YYYY-MM-DD HH:MM:SS`
    - `days since YYYY-MM-DD`
- `calendar`: A recognized calendar identifier, like `proleptic_gregorian`, `gregorian` etc.

## Cutmaps: a NetCDF files cookie-cutter


### Usage:
TODO

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