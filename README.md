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

Lisflood Utilities is a set of tools to help LISFLOOD users (or any users of PCRaster/NetCDF files)
to execute some mundane tasks that are necessary to operate lisflood.
Here's a list of utilities you can find in lisflood-utilities package.

* __[catchstats](#catchstats)__ calculates catchment statistics (mean, sum, std, min, max...) from NetCDF4 files given masks created with [`cutmaps`](#cutmaps:-a-NetCDF-files-cookie-cutter).

* __[cddmap](#cddmap)__ is a tool to generate correlation decay distance (CDD) maps starting from station timeseries

* __[compare](#compare)__ is a package containing a set of simple Python classes that helps to compare 
NetCDF, PCRaster and TSS files.

* __[cutmaps](#cutmaps)__ is a tool to cut NetCDF files in order to reduce size, using either
  - a bounding box of coordinates
  - a bounding box of matrix indices
  - an existing boolean area mask
  - a list of stations and a LDD ("local drain direction" in NetCDF or PCRaster format)

* __[gridding](#gridding)__ is a tool to interpolate meteo variables observations stored in text files containing (lon, lat, value) into grids.
  - uses inverse distance interpolation
  - input file names must use format: \<var\>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt
  - option to store all interpolated grids in a single NetCDF4 file
  - option to store each interpolated grid in a GeoTIFF file
  - output files are compressed
  - grids are setup in the configuration folder and are defined by a dem.nc file
  - meteo variables parameters are defined in the same configuration folder

* __[nc2pcr](#nc2pcr)__ is a tool to convert a NetCDF file into PCRaster maps.
  - convert 2D variables in single PCRaster maps
  - NetCDF4 mapstacks are not supported yet

* __[ncextract](#ncextract)__ is a tool to extract values from NetCDF4 (or GRIB) file(s) at specific coordinates.

* __[pcr2nc](#pcr2nc)__ is a tool to convert PCRaster maps to NetCDF4 files.
  - convert a single map into a NetCDF4 file
  - convert a time series of maps into a NetCDF4 mapstack
  - support for WGS84 and ETRS89 (LAEA) reference systems
  - fine tuning of output files (compression, significant digits, etc.)
  
> **Note**: PCRaster must be installed in the Conda environment.

* __[thresholds](#thresholds)__ is a tool to compute the discharge return period thresholds from NetCDF4 file containing a discharge time series.

* __[water-demand-historic](#water-demand-historic)__ is a package allowing to generate sectoral (domestic, livestock, industry, and thermoelectric) water demand maps with monthly to yearly temporal steps for a range of past years, at the users’ defined spatial resolution and geographical extent. These maps are required by the LISFLOOD OS [water use module](https://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use)

* __[waterregions](#waterregions)__ is a package containing two scripts that allow to create and verify a water regions map, respectively.

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
- NetCDF4 C library

#### Install
If you use conda, create a new env (or use an existing one) and install gdal and lisflood-utilities:

```bash
conda create --name myenv python=3.7 -c conda-forge
conda activate myenv
conda install -c conda-forge pcraster eccodes gdal
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

Note: if you previously installed an older version of the lisflood-utilities, it is highly recommended to remove it before installing the newest version:

```bash
pip uninstall lisflood-utilities
pip install -e./
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
  variable_x_name: lon
  variable_y_name: lat
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

In the `geographical` section you can configure `datum` and name of the x and y variables. As `variable_x_name` and `variable_y_name` you should use 'lon' and 'lat' for geographical coordinates (e.g. WGS84) and 'x' and 'y' for projected coordinates (e.g. ETRS89).
 
Currently, pcr2nc supports the following list for `datum`:

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

This tool converts single maps NetCDF (time dimension is not supported yet) to PCRaster format.

### Usage

```bash
nc2pcr -i /path/to/input/file.nc -o /path/to/output/out.map [-c /path/to/clone.map optional]
```

If input file is a LDD map, you must add the `-l` flag:

```bash
nc2pcr -i /path/to/input/ldd.nc -o /path/to/output/ldd.map  -l [-c /path/to/clone.map optional]
```

## cutmaps

This tool cuts NetCDF files using either a mask, a bounding box, or a list of stations along with a LDD (local drain direction) map.  

### Usage

#### In the command prompt

The tool requires a series of arguments:

* The area to be extracted can be defined in one of the following ways:
    - `-m`, `--mask`: a mask map (either PCRaster or NetCDF format).
    - `-i`, `--cuts_indices`: a bounding box defined by matrix indices in the form `-i imin imax jmin jmax` (the indices must be integers).
    - `-c`, `--cuts`: a bounding box defined by coordinates in the form `-c xmin xmax ymin ymax` (the coordinates can be integer or floating point numbers; x = longitude, y = latitude).
    - `-N`, `-stations`: a list of stations included in a tab separated text file. This approach requires a LDD (local drain direction) map as an extra input, defined with the argument `-l` (`-ldd`).
* The files to be cut may be defined in one of the following ways:
    - `-f`, `--folder`: a folder containing NetCDF files.
    - `-F`, `--file`: a single netCDF file to be cut.
    - `-S`, `--subdir`: a directory containing a number of folders
* The resulting files will be saved in the folder defined by the argument `-o` ( `--outpath`).

There are additional optional arguments

* `-W`, `--overwrite`: it allows to overwrite results.
* `-C`, `--clonemap`: it can be used to define a clone map when the LDD input map (argument `-l`) is in NetCDF format.

##### Examples 

**Using a mask**

The following command will cut all NetCDF files inside a specific folder (argument `-f`) using a mask (argument `-m`). The mask is a boolean map (1 only in the area of interes) that `cutmaps` uses to create a bounding box. The resulting files will be written in the current folder (argument `-o`). 

```bash
cutmaps -m /workarea/Madeira/maps/MaskMap/Bacia_madeira.nc -f /workarea/Madeira/lai/ -o ./
```

**Using indices**

The following command cuts all the maps in an input folder (argument `-f`) using a bounding box defined by matrix indices (argument `-i`). Knowing your area of interest from your NetCDF files, you can determine indices of the array and pass them in the form `-i imin imax jmin jmax` (integer numbers).

```bash
cutmaps -i "150 350 80 180" -f /workarea/Madeira/lai/ -o ./
```

**Using coordinates**

The following command cuts all the maps in an input directory containing several folders (argument `-S`) using a bounding box defined by coordinates (argument `-c`). The argument `-W` allows to overwrite pre-existing files in the output directory (argument `-o`):

```bash
cutmaps -S /home/projects/lisflood-eu -c "4078546.12 4463723.85 811206.57 1587655.50" -o /Work/Tunisia/cutmaps -W
```

**Using station coordinates and a local drain direction map**

The TXT file with stations must have a specific format as in the example below. Each row represents a stations, and it contains three columns separated by tabs that indicated the X and Y coordinates (or lon and lat) and the ID of the station.

```text
4297500	1572500	1
4292500	1557500	2
4237500	1537500	3
4312500	1482500	4
4187500	1492500	5
```

The following command will cut all the maps in a specific folder (`-f` argument) given a LDD map (`-l` argument) and the previous text file (`-N` argument), and save the results in a folder defined by the argument `-o`.

```bash
cutmaps -f /home/projects/lisflood-eu -l ldd.map -N stations.txt -o /Work/Tunisia/cutmaps
```

If the LDD is in NetCDF format, it will be first converted into PCRaster format.

```bash
cutmaps -f /home/projects/lisflood-eu -l ldd.nc -N stations.txt -o /Work/Tunisia/cutmaps
``` 

If you experience problems, you can try to pass a path to a PCRaster clone map using the `-C` argument.

```bash
cutmaps -f /home/projects/lisflood-eu -l ldd.nc -C area.map -N stations.txt -o /Work/Tunisia/cutmaps
```

#### In a Python script

You can use the `cutmaps` tool programmatically from a Python code. You just need to import the tool and then run the following command:

```Python
from lisfloodutilities.cutmaps.cutlib import cutmap

cutmap(file_to_cut, file_out, x_min, x_max, y_min, y_max, use_coords=True)
```

Where `use_coords=True` means that it will cut the map using coordinates, while `use_coords=False` will use indices.

### Output

Apart from the cut files in the output folder specified in the command prompt, `cutmaps` produces other outputs in the folder where the LDD map is stored:

* _mask.map_ and _mask.nc_ for your area of interest, which may be needed in subsequent LISFLOOD/LISVAP executions.
* _outlets.map_ and _outlets.nc_ based on _stations.txt_, which will let you produce gauges TSS if configured in LISFLOOD.


## compare

This tool compares two NetCDF datasets. You can configure it with tolerances (absolute `--atol`, relative `--rtol`, thresholds for percentage of tolerated different values `--max-diff-percentage`). You can also set the option `--save-diffs` to write files with the diffences, so that you can inspect maps and differences with tools like [Panoply](https://www.giss.nasa.gov/tools/panoply/).

```text
usage: compare [-h] -a DATASET_A -b DATASET_B -m MASKAREA [-M SUBMASKAREA]
               [-e] [-s] [-D] [-r RTOL] [-t ATOL] [-p MAX_DIFF_PERCENTAGE]
               [-l MAX_LARGEDIFF_PERCENTAGE]

Compare NetCDF outputs: 0.12.12

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
  -D, --save-diffs      flag to save diffs in NetCDF files for visual
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

## thresholds

The thresholds tool computes the discharge return period thresholds using the method of L-moments.
It is used to post-process the discharge from the LISFLOOD long term run.
The resulting thresholds can be used in a flood forecasting system to define the flood warning levels.

### Usage:
The tool takes as input a Netcdf file containing the annual maxima of the discharge signal. LISFLOOD computes time series of discharge values (average value over the selected computational time step). The users are therefore required to compute the annual maxima. As an example, this step can be achieved by using CDO (cdo yearmax), for all the details please refer to [https://code.mpimet.mpg.de/projects/cdo/embedded/index.html#x1-190001.2.5](https://code.mpimet.mpg.de/projects/cdo/embedded/index.html#x1-190001.2.5)

The output NetCDF file contains the following return period thresholds [1.5, 2, 5, 10, 20, 50, 100, 200, 500], together with the Gumbel parameters (sigma and mu).

```text
usage: thresholds [-h] [-d DISCHARGE] [-o OUTPUT]

Utility to compute the discharge return period thresholds using the method of L-moments.
Thresholds computed: [1.5, 2, 5, 10, 20, 50, 100, 200, 500]

options:
  -h, --help            show this help message and exit
  -d DISCHARGE, --discharge DISCHARGE
                        Input discharge files (annual maxima)
  -o OUTPUT, --output OUTPUT
                        Output thresholds file
```

## water-demand-historic

This utility allows to create water demand maps at the desired resolution and for the desired geographical areas. The maps indicate, for each pixel, the time-varying water demand map to supply for domestic, livestock, industrial, and thermoelectric water consumption. The temporal discretization is monthly for domestic and energy demand, yearly for industrial and livestock demand. The maps of sectoral water demand are required by the LISFLOOD OS [water use module](https://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use/). Clearly, the sectoral water demand maps and the scripts of this utility can be used also for other applications, as well as for stand-alone analysis of historical water demand for anthropogenic use.

#### Input
The creation of the sectoral water demand maps requires a template map that defines the desired geographical area and spatial resolution. The generation of the maps relies on a number of external datasets (examples are the [Global Human Settlement - Datasets - European Commission (europa.eu)](https://ghsl.jrc.ec.europa.eu/datasets.php) and [FAO AQUASTAT Dissemination System](https://data.apps.fao.org/aquastat/?lang=en)). The locations of the template map, of the input datasets and files, of the output folder, and other users’ choices (e.g. start year and end year) are specified in a configuration file. The syntax of the configuration file is pre-defined and an example is provided to the users. The complete list of external datasets, the instructions on how to prepare (i) the external dataset, (ii) the template map, (iii) the input folder, (iv) the output folder, and (v) the configuration file are explained into details [here](src/lisfloodutilities/water-demand-historic/README.md)

#### Output
Four sectoral water demand maps in netCDF-4 format. The geographical extent and the spatial resolution are defined by the template map (users-defined input file). Each netCDF-4 file has 12 months, for each year included in the temporal time span identified by the user. Sectoral water demand data with lower (yearly) temporal resolution are repeated 12 times. 

#### Usage
The methodology includes five main steps. The instructions on how to retrieve the scrips, create the environment including all the required packages, and use the utility are provided [here](src/lisfloodutilities/water-demand-historic/README.md)

#### Important notes on documentation and data availability
The complete list of external datasets, the instructions on how to retrieve the external datasets, the methodology, and the usage of the scripts are explained into details [here](src/lisfloodutilities/water-demand-historic/README.md). The README file provides detailed technical information about the input datasets and the usage of this utility. The methodology is explained in the manuscript: Choulga, M., Moschini, F., Mazzetti, C., Grimaldi, S., Disperati, J., Beck, H., Salamon, P., and Prudhomme, C.: Technical note: Surface fields for global environmental modelling, EGUsphere, 2023 ([preprint](https://doi.org/10.5194/egusphere-2023-1306)).

The global sectoral water demand maps at 3 arcmin (or 0.05 degrees) resolution, 1979-2019, produced using the scripts of this utility can be downloaded from [Joint Research Centre Data Catalogue - LISFLOOD static and parameter maps for GloFAS - European Commission (europa.eu)](https://data.jrc.ec.europa.eu/dataset/68050d73-9c06-499c-a441-dc5053cb0c86)



## waterregions

The modelling of water abstraction for domestic, industrial, energetic, agricoltural and livestock use   can require a map of the water regions. The concept of water regions and information for their definition are explained [here](htpst://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use/). 
Since groundwater and surface water resources demand and abstraction are spatially distributed inside each water region, each model set-up must include all the pixels of the water region. This requirement is crucial for the succes of the calibration of the model. This utility allows the user to meet this requirement.
More specifically, this utility can be used to:
1. create a water region map which is consistent with a set of calibration points: this purpose is achieved by using the script define_waterregions.
2. verify the consistency between an existing water region map and an exixting map of calibration catchments: this purpose is achieved by using the script verify_waterregions
It is here reminded that when calibrating a catchment which is a subset of a larger computational domain, and the option wateruse is switched on, then the option groudwatersmooth must be switched off. The explanation of this requirement is provided in the chapter [Water use](https://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use/) of the LISFLOOD documentation. 

#### Requirements
python3, pcraster 4.3. The protocol was tested on Linux.

### define_waterregions
This utility allows to create a  water region map which is consistent with a set of calibration points. The protocol was created by Ad De Roo (Unit D2, Joint Research Centre).

#### Input 
- List of the coordinates of the calibration points. This list must be provided in a .txt file with three columns: LONGITUDE(or x), LATITUDE(or y), point ID.
- LDD map can be in NetCDF format or pcraster format. When using pcraster format, the following condition must be satisfied: *PCRASTER_VALUESCALE=VS_LDD*. 
- Countries map in NetCDF format or pcraster format. When using pcraster format, the following condition must be satisfied: *PCRASTER_VALUESCALE=VS_NOMINAL*. This map shows the political boundaries of the Countries, each Coutry is identified by using a unique ID. This map is used to ensure that the water regions are not split accross different Countries.
- Map of the initial definition of the water regions in NetCDF format or pcraster format. When using pcraster format, the following condition must be satisfied: *PCRASTER_VALUESCALE=VS_NOMINAL*. This map is used to attribute a water region to areas not included in the calibration catchments. In order to create this map, the user can follow the guidelines provided [here](https://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use/).
- file *.yaml* or *.json* to define the metadata of the output water regions map in NetCDF format. An example of the structure of these files is provided [here](tests/data/waterregions)

##### Input data provided by this utility:
This utility provides three maps of [Countries IDs](tests/data/waterregions): 1arcmin map of Europe (EFAS computational domain), 0.1 degree and 3arcmin maps of the Globe. ACKNOWLEDGEMENTS: both the rasters were retrieved by upsampling the original of the World Borders Datase provided by  http://thematicmapping.org/ (the dataset is available under a Creative Commons Attribution-Share Alike License).

#### Output
Map of the water regions which is consistent with the calibration catchments. In other words, each water region is entirely included in one calibration catchment.  The test to check the consistency between the newly created water regions map and the calibration catchments is implemented internally by the code and the outcome of the test is printed on the screen. 
In the output map, each water region is identified by a unique ID. The format of the output map can be NetCDF or pcraster.

#### Usage
The following command lines allow to produce a water region map which is consistent with the calibration points (only one commad line is required: each one of the command lines below shows a different combination of input files format):

*python define_waterregions.py -p calib_points_test.txt -l ldd_test.map -C countries_id_test.map -w waterregions_initial_test.map -o my_new_waterregions.map* <br>

*python define_waterregions.py -p calib_points_test.txt -l ldd_test.nc -C countries_id_test.nc -w waterregions_initial_test.nc -o my_new_waterregions.nc -m metadata.test.json* <br>

*python define_waterregions.py -p calib_points_test.txt -l ldd_test.map -C countries_id_test.nc -w waterregions_initial_test.map -o my_new_waterregions.nc -m metadata.test.yaml* <br>


The input maps can be in nectdf format or pcraster format (the same command line can accept a mix of pcraster and NetCDF formats).It is imperative to write the file name in full, that is including the extension (which can be either ".nc" or ".map").<br>
The utility can return either a pcraster file or a NetCDF file. The users select their preferred format by specifying the extension of the file in the output option (i.e. either ".nc" or ".map"). <br>
The metadata file in .yaml format must be provided only if the output file is in NetCDF format.<br>

The code internally verifies that the each one of the newly created water regions is entirely included  within one calibration catchments. If this condition is satisfied, the follwing message in printed out: *“OK! Each water region is completely included inside one calibration catchment”*. If the condition is not satisfied, the error message is *“ERROR: The  water regions WR are included in more than one calibration catchment”*. Moreover, the code provides the list of the water regions WR and the calibration catchments that do not meet the requirment. This error highlight a problem in the input data: the user is recommended to check (and correct) the list of calibration points and the input maps.

The input and output arguments are listed below. 


```text
usage: define_waterregions.py [-h] -p CALIB_POINTS -l LDD -C COUNTRIES_ID -w
                              WATERREGIONS_INITIAL -o OUTPUT_WR

Define Water Regions consistent with calibration points: {}

optional arguments:
  -h, --help            show this help message and exit
  -p CALIB_POINTS, --calib_points CALIB_POINTS
                        list of calibration points: lon or x, lat or y, point id. File extension: .txt,
  -l LDD, --ldd LDD     LDD map, file extension: .nc or .map
  -C COUNTRIES_ID, --countries_id COUNTRIES_ID
                        map of Countries ID, fike extension .nc or .map 
  -w WATERREGIONS_INITIAL, --waterregions_initial WATERREGIONS_INITIAL
                        initial map of water regions, file extension: .nc or .map
  -o OUTPUT_WR, --output_wr OUTPUT_WR
                        output map of water regions, file extension: .nc or .map 
  -m METADATA, --metadata_file METADATA
                        Path to metadata file for NetCDF, .yaml or .json format                     
```



### verify_waterregions

This function allows to verify the consistency between a  water region map and a map of calibration catchments. This function must be used when the water region map and the map of calibration catchments have been defined in an independent manner (i.e. not using the utility **define_waterregions**). The function verify_waterregions verifies that each water region map is entirely included in one calibration catchment. If this condition is not satisfied, an error message is printed on the screen. 

#### Input
- Map of calibration catchments in NetCDF format.
- Water regions map in NetCDF format.

#### Output
The output is a message on the screen. There are two options:
- 'OK! Each water region is completely included inside one calibration catchment.'
- 'ERROR: The  water regions WR are included in more than one calibration catchment’: this message is followed by the list of the water regions and of the catchment that raised the isuue.
In case of error message, the user can implement the function **define_waterregions**.

#### Usage
The following command line allows to produce a water region map which is consistent with the calibration points:

*python verify_waterregions.py -cc calib_catchments_test.nc -wr waterregions_test.nc*

The input and output arguments are listed below. All the inputs are required. 

```text
usage: verify_waterregions.py [-h] -cc CALIB_CATCHMENTS -wr WATERREGIONS

Verify that the Water Regions map is consistent with the map of the
calibration catchments

optional arguments:
  -h, --help            show this help message and exit
  -cc CALIB_CATCHMENTS, --calib_catchments CALIB_CATCHMENTS
                        map of calibration catchments, NetCDF format
  -wr WATERREGIONS, --waterregions WATERREGIONS
                        map of water regions, NetCDF format
```

NOTE:
The utility **pcr2nc** can be used to convert a map in pcraster format into NetCDF format.


## gridding

This tool is used to interpolate meteo variables observations stored in text files containing (lon, lat, value) into grids.
It uses inverse distance interpolation method from pyg2p.

#### Requirements
python3, pyg2p

### Usage

> __Note:__ This guide assumes you have installed the program with pip tool.
> If you cloned the source code instead, just substitute the executable `gridding` with `python bin/gridding` that is in the root folder of the cloned project.

The tool requires four mandatory command line input arguments:

- -i, --in: Set input folder path with kiwis/point files
- -o, --out: Set output folder base path for the tiff files or the NetCDF file path.
- -c, --conf: Set the grid configuration type to use. Right now only 5x5km, 1arcmin are available.
- -v, --var: Set the variable to be processed. Right now only variables pr,pd,tn,tx,ws,rg,pr6,ta6 are available.

The input folder must contain the meteo observation in text files with file name format:  \<var\>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt
The files must contain the columns longitude, latitude, observation_value is separated by TAB and without the header.
Not mandatory but could help to store the files in a folder structure like: ./YYYY/MM/DD/\<var\>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt

Example of command that will generate a NetCDF file containing the precipitation (pr) grids for March 2023:

```bash
gridding -i /meteo/pr/2023/ -o /meteo/pr/pr_MARCH_2023.nc -c 1arcmin -v pr -s 202303010600 -e 202304010600
```

The input and output arguments are listed below and can be seen by using the help flag.

```bash
gridding --help
```

```text
usage: gridding [-h] -i input_folder -o {output_folder, NetCDF_file} -c
                {5x5km, 1arcmin,...} -v {pr,pd,tn,tx,ws,rg,...}
                [-d files2process.txt] [-s YYYYMMDDHHMISS] [-e YYYYMMDDHHMISS]
                [-q] [-t] [-f]

version v0.1 ($Mar 28, 2023 16:01:00$) This script interpolates meteo input
variables data into either a single NETCDF4 file or one GEOTIFF file per
timestep. The resulting NetCDF is CF-1.6 compliant.

optional arguments:
  -h, --help            show this help message and exit
  -i input_folder, --in input_folder
                        Set input folder path with kiwis/point files
  -o {output_folder, NetCDF_file}, --out {output_folder, NetCDF_file}
                        Set output folder base path for the tiff files or the
                        NetCDF file path.
  -c {5x5km, 1arcmin,...}, --conf {5x5km, 1arcmin,...}
                        Set the grid configuration type to use.
  -v {pr,pd,tn,tx,ws,rg,...}, --var {pr,pd,tn,tx,ws,rg,...}
                        Set the variable to be processed.
  -d files2process.txt, --dates files2process.txt
                        Set file containing a list of filenames to be
                        processed in the form of
                        <var>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt
  -s YYYYMMDDHHMISS, --start YYYYMMDDHHMISS
                        Set the start date and time from which data is
                        imported [default: date defining the time units inside
                        the config file]
  -e YYYYMMDDHHMISS, --end YYYYMMDDHHMISS
                        Set the end date and time until which data is imported
                        [default: 20230421060000]
  -q, --quiet           Set script output into quiet mode [default: False]
  -t, --tiff            Outputs a tiff file per timestep instead of the
                        default single NetCDF [default: False]
  -f, --force           Force write to existing file. TIFF files will be
                        overwritten and NetCDF file will be appended.
                        [default: False]
```


## cddmap

This tool is used to generate correlation decay distance (CDD) maps starting from station timeseries

#### Requirements
python3, pyg2p

### Usage

cddmap [directory]/[--analyze]/[--merge-and-filter-jsons]/--generatemap] [--start first_station] [--end last_station] [--parallel] [--only-extract-timeseries timeseries_keys_file] [--maxdistance max_distance_in_km]

The tool requires an input argument indicating the station timeseries main folder, and calculates the CDD for each stations as well as correlations and distances files. Outputs the results in a txt file containing station coordinates and CDD values.
After creating the CDD txt file, it can be used with one of the following commands:

- --analyze: read cdd file previously created for postprocessing  
- --merge-and-filter-jsons: merge all cdd files in a folder and filters out a list of stations.
- --generatemap: generate a NetCDF CDD map file using CDD txt file and angular distance weighted interpolation between station points
- --start and --end arguments are used to split the task in many sub tesks, evaluating only the stations between "start" and "end", since the CDD evaluation can be very time-demanding. 
- --only-extract-timeseries: in combination with path of the station's main folder, extracts the timeseries specified in the timeseries_keys_file txt list of keys
- --parallel: enable CDD evaluation in parallel on multiple cores. It will require more memory
- --maxdistance: evaluates only station that are clores then maxdistance in km

The input folder must contain the meteo observation in text files

Example of command that will generate txt files for the CDD of precipitation (pr), in parallel mode, for station that are closer then 500 kms:

```bash
cddmap /meteo/pr --parallel --maxdistance 500
```


## ncextract

The `ncextract` tool extracts time series from (multiple) NetCDF or GRIB file(s) at user defined coordinates.

### Usage

The tool takes as input a CSV file containing point coordinates and a directory containing one or more NetCDF or GRIB files. The CSV files must contain only three columns: point identifier, and its two coordinates. The name of the coordinate fields must match those in the NetCDF or GRIB files. For example:

```text
ID,lat,lon
0010,40.6083,-4.2250
0025,37.5250,-6.2750
0033,40.5257,-6.4753
```

The output is a file containing the time series at the pixels corresponding to the provided coordinates, in chronological order. The function supports two otput formats: CSV or NetCDF.

```text
usage: ncextract.py [-h] -i INPUT -d DIRECTORY -o OUTPUT [-nc]

Utility to extract time series of values from (multiple) NetCDF files at specific coordinates.
Coordinates of points of interest must be included in a CSV file with at least 3 columns named id,
lat, lon.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV file (id, lat, lon)
  -d DIRECTORY, --directory DIRECTORY
                        Input directory with .nc files
  -o OUTPUT, --output OUTPUT
                        Output file. Two extensions are supported: .csv or .nc
```

#### Use in the command prompt

The following command extracts the discharge time series from EFAS simulations (NetCDF files in the directory _EFAS5/discharge/maps_) in a series of points where gauging stations are located (file _gauging_stations.csv_), and saves the extraction as a CSV file.

```bash
ncextract -i ./gauging_stations.csv -d ./EFAS5/discharge/maps/ -o ./EFAS5/discharge/timeseries/results_ncextract.csv
```

#### Use programmatically

The function can be imported in a Python script. It takes as inputs two `xarray.Dataset`: one defining the input maps and the other the points of interest. The result of the extraction can be another `xarray.Dataset`, or saved as a file either in CSV or NetCDF format.

```Python
from lisfloodutilities.ncextract import extract_timeseries

# load desired input maps and points
# maps: xarray.Dataset
# points: xarray.Dataset

# extract time series and save in a xarray.Dataset
ds = extract_timeseries(maps, points, output=None)
```


## catchstats

The `catchstats` tool calculates catchment statistics given a set of input NetCDF files and a set of mask NetCDF files.

### Usage

#### In the command prompt

The tool takes as input a directory containing the NetCDF files from which the statistics will be computed, and another directory containing the NetCDF files that define the catchment boundaries, which can be any of the outputs of  `cutmaps` (not necessarily the file _my_mask.nc_). The input files can be the LISFLOOD static maps (no temporal dimension) or stacks of maps with a temporal dimension. The mask NetCDF files must be named after the catchment ID, as this name will be used to identify the catchment in the output NetCDF. For instance, _0142.nc_ would correspond to the mask of catchment 142. Optionally, an extra NetCDF file can be passed to the tool to account for different pixel area; in this case, the statistics will be weighted by this pixel area map.

Only some statistics are currently available: mean, sum, std (standard deviation), var (variance), min, max, median and count. The weighing based on pixel area does not affect the statistics min, max, median nor count.

The output are NetCDF files (as many as catchments in the mask directory) containing the resulting statistics.

```text
usage: catchstats.py [-h] -i INPUT -m MASK -s STATISTIC -o OUTPUT -a AREA [-W]

Utility to compute catchment statistics from (multiple) NetCDF files.
The mask map is a NetCDF file with values in the area of interest and NaN elsewhere.
The area map is optional and accounts for varying pixel area with latitude.

options:
  -h, --help
                          show this help message and exit
  -i INPUT, --input INPUT
                          directory containing the input NetCDF files
  -m MASK, --mask MASK
                          directory containing the mask NetCDF files
  -s STATISTIC, --statistic STATISTIC
                          list of statistics to be computed. Possible values: mean, sum, std, var, min, max, median, count
  -o OUTPUT, --output OUTPUT
                          directory where the output NetCDF files will be saved
  -a AREA, --area AREA
                          NetCDF file of pixel area used to weigh the statistics
  -W, --overwrite
                          overwrite existing output files
```

**Example**

The following command calculates the average and total precipitation for a set of catchemtns from the dataset EMO-1. The static map _pixarea.nc_ is used to account for varying pixel area.

```bash
catchstats -i ./EMO1/pr/ -m ./masks/ -s mean sum -o ./areal_precipitation/ -a ./EFAS5/static_maps/pixarea.nc
```

#### In a Python script

The tool can be imported in a Python script to be able to save in memory the output. This function takes in a `xarray.Dataset` with the input maps from which statistics will be computed, a dictionary of `xarray.DataArray` with the catchment masks, and optionally the weighing map. By default, the result is a `xarray.Dataset`, but NetCDF files could be written, instead, if a directory is provided in the `output` attribute.

```Python
# import function
from lisfloodutilities.catchstats import catchment_statistics

# load desired input maps and catchment masks
# maps: xarray.Dataset
# masks: Dict[int, xarray.DataArray]

# compute statistics and save in a xarray.Dataset
ds = catchment_statistics(maps, masks, statistic=['mean', 'sum'], weight=None, output=None)
```

### Ouput

The structure of the output depends on whether the input files include a temporal dimension or not:

* If the input files DO NOT have a time dimension, the output has a single dimension: the catchment ID. It contains as many variables as the combinations of input variables and statistics. For instance, if the input variables are "elevation" and "gradient" and three statistics are required ("mean", "max", "min"), the output will contain 6 variables: "elevation_mean", "elevation_max", "elevation_min", "gradient_mean", "gradient_max" and "gradient_min".
* If the input files DO have a time dimension, the output has two dimensions: the catchment ID and time. The number of variables follows the same structure explained in the previous point. For instance, if the input files are daily maps of precipitation (variable name "pr") and we calculate the mean and total precipitation over the catchment, the output will contain two dimensions ("ID", "time") and two variables ("pr_mean", "pr_sum").


## lfcoords

This tool finds the appropriate coordinates in the LISFLOOD river network of any point, provided that the catchment area is known. A thourough explanation of the method can be found in [Burek and Smilovic (2023)](https://essd.copernicus.org/articles/15/5617/2023/).

First, it uses the original coordinates and catchment area to find the most accurate pixel in a high-resolution map. [Burek and Smilovic (2023)](https://essd.copernicus.org/articles/15/5617/2023/) use MERIT [(Yamazaki et al., 2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR024873), which has a spatial resolution of 3 arc-seconds. The result of this first step is, for every point, a new value of coordinates and area, and a shapefile of the catchment polygon in high-resolution.

Second, it finds the pixel in the low-resolution grid (LISFLOOD static maps) that better matches the catchment shape derived in the previous step. As a result, for each point we obtained a new value of coordinates and area, and a new shapefile of the catchment polygon in low-resolution.

### Usage

#### In the command prompt

The tool can be executed from the command prompt by indicating a configuration file.

```text
usage: lfcoords.py [-h] -c CONFIG_FILE

Correct the coordinates of a set of stations to match the river network in the
LISFLOOD static map.
First, it uses a reference value of catchment area to find the most accurate
pixel in a high-resolution map.
Second, it finds the pixel in the low-resolution map that better matches the
catchment shape derived from the high-resolution map.

options:
  -h, --help
                          show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG FILE
                          path to the YML configuration file
```

**Example**

```bash
lfcoords --config-file config.yml
```

##### Configuration file

The configuration file defines the input files, the folder where the resulting shapefiles will be saved, and some thresholds used in the process. A template of the configuration file can be found [here](./src/lisfloodutilities/lfcoords/config.yml). Below you find an example:

```yml
input:
    stations: points.csv # ID, lat, lon, area
    ldd_fine: ldd_MERIT.tif
    upstream_fine: uparea_MERIT.tif # km2
    ldd_coarse: ldd_LISFLOOD.tif
    upstream_coarse: uparea_LISFLOOD.nc # m2
            
output_folder: ./shapefiles/

conditions:
    min_area: 25 # km2
    abs_error: 50 # km2
    pct_error: 1 # %
```

##### Inputs

The tool requires 5 inputs:

* A CSV file of the stations to be located in the LISFLOOD grid. This file must contain four columns with four specific names: 'ID', 'area' in km2, 'lat', 'lon'. Below you find an example of the stations CSV file:

```csv
ID,area,lat,lon
429,35399,49.018,12.144
436,26448,48.947,12.015
439,37687,48.88,12.747
```

* A map of the local drainage directions in high-resolution, e.g., MERIT.
* A map of the upstream area in high-resolution, e.g., MERIT. The units of this map must be km2, same units as the _area_ field in the CSV file.
* A map of the local drainage directions in low-resolution, i.e., the LISFLOOD static map.
* A map of the upstream area in low-resolution, i.e., the LISFLOOD static map. The units of this map are m2 (instead of km2), as these are the units used in LISFLOOD; the code converts internally this map into km2.

All maps can be provided either in TIFF or NetCDF format.

##### Outputs

The main output is a new **CSV file** saved in the same directory as the input CSV file and named similarly, but with a suffix indicating the resolution of the LISFLOOD grid. For instance, in the configuration file above the input CSV file is named _stations.csv_ and the resolution of the LISFLOOD grid is 3 arcmin, so the output CSV files will be named _stations_3min.csv_. The CSV contains 6 new columns defining the coordinates and catchment area in both the high-resolution (`3sec` in the example) and low-resolution grids (`3min` in the example). Example:

```csv
ID,area,area_3min,area_3sec,lat,lat_3min,lat_3sec,lon,lon_3min,lon_3sec
429,35399,35216,35344,49.018,49.025,49.022083,12.144,12.125,12.14375
436,26448,26334,26394,48.947,48.925,48.94625,12.015,12.025,12.014583
439,37687,37540,37605,48.88,48.925,48.879583,12.747,12.675,12.74625
```

Besides, the tool creates **shapefiles** of the catchment polygons derived for both the high and low resolution grids. The shapefiles are saved in two subdirectories inside the `output_folder` directory defined in the configuration file. In each of these subdirectories, there will be one file for each station.

#### In a Python script

The functions that compose the `lfcoords` tool can be imported in a Python script to be able to save in memory the output:

```Python
# import function
from lisfloodutilities.lfcoords import Config
from lisfloodutilities.lfcoords.finer_grid import coordinates_fine
from lisfloodutilities.lfcoords.coarser_grid import coordinates_coarse

# read the configuration file
cfg = Config(CONFIG_FILE)

# find coordinates in high-resolution (HR) grid
points_HR = coordinates_fine(cfg, save=False)

# find coordinates in LISFLOOD
points_LF = coordinates_coarse(cfg, points_HR, save=False)
```

## Using `lisfloodutilities` programmatically 

You can use lisflood utilities in your python programs. As an example, the script below creates the mask map for a set of stations (stations.txt). The mask map is a boolean map with 1 and 0. 1 is used for all (and only) the pixels hydrologically connected to one of the stations. The resulting mask map is in pcraster format.

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
