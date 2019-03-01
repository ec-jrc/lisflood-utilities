# pcr2nc

__pcr2nc__ is a tool to convert PCRaster maps to netCDF4 files.
It's developed at JRC E1 directorate, as a Floods group initiative.

- convert a single map into a NetCDF4 file
- convert a time series of maps into a NetCDF4 mapstack
- support for WGS84 and ETRS89 (LAEA) reference systems
- fine tuning of output files (compression, significant digits, etc.)

## Installation

There are two ways to install the software: one via pip tool and one by cloning the repository.

### Requisites

Ensure you have properly installed the following software:

- Python 3.5+
- GDAL C library and software
- netCDF4 C library

Create a python3.5 virtualenv for using the software and activate it.

If you have virtualenvwrapper:
```bash
$ workon pcr2nc
```

Otherwise just execute the activate script
```bash
$ source /path/to/virtualenvs/pcr2nc/bin/activate
```


### Install by cloning the repository
Assuming your virtualenv is called pcr2nc and you have virtualenvwrapper installed:

```bash
$ git clone https://github.com/domeniconappo/pcr2nc.git
$ cd pcr2nc
```

Install requirements
```bash
$ pip install -r requirements.txt
```

If GDAL library fails to install, ensure to install the same package version of the
library you have on your system.
Example: you have installed gdal 2.1, then:

```bash
$ pip install GDAL==2.1
$ pip install -r requirements.txt
```

### Install via pip tool

Activate the virtualenv and then:

```bash
$ pip install git+https://github.com/domeniconappo/pcr2nc.git
```

After the install was complete, you still have to install the proper GDAL package,
according to the version of gdal library that is installed on your machine.

E.g.
```bash
$ pip install GDAL==2.1
```

## Usage

> __Note:__ This guide assumes you have installed the program with pip tool.
> If you cloned the source code instead, just substitute the executable `pcr2nc` with `python pcr2nc_script.py` that is in the root folder of the cloned project.

The tool takes three command line input arguments:

- -i, --input: It can be a path to a single file, a folder or a unix-like widlcard expression like _/path/to/files/dis00*_
- -o, --output_file: Path to the output nc file
- -m, --metadata: Path to a yaml or json file containing configuration for the NetCDF4 output file.

Unless the input is a single file, the resulting NetCDF4 file will be a mapstack according to a time dimension.

#### Example of usages:

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

### Writing metadata configuration file

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

In `geographical` section the only setting to configure is `datum`. Currently, `WGS84`, `ETRS89` and `GISCO` are the reference systems that are supported by pcr2nc.

#### Time section

This section is optional and is only required if the output file is a mapstack.
In this section you have to configure `units` and `calendar`.

- `units`: Can be one of the following strings (replacing placeholders with the actual date):
    - `hours since YYYY-MM-DD HH:MM:SS`
    - `days since YYYY-MM-DD`
- `calendar`: A recognized calendar identifier, like `proleptic_gregorian`, `gregorian` etc.
