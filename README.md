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

* __[pcr2nc](#pcr2nc)__ is a tool to convert PCRaster maps to NetCDF4 files.
  - convert a single map into a NetCDF4 file
  - convert a time series of maps into a NetCDF4 mapstack
  - support for WGS84 and ETRS89 (LAEA) reference systems
  - fine tuning of output files (compression, significant digits, etc.)
 
* __[nc2pcr](#nc2pcr)__ is a tool to convert a NetCDF file into PCRaster maps.
  - convert 2D variables in single PCRaster maps
  - NetCDF4 mapstacks are not supported yet

* __[cutmaps](#cutmaps)__ is a tool to cut NetCDF files in order to reduce size, using either
  - a bounding box of coordinates
  - a bounding box of matrix indices
  - an existing boolean area mask
  - a list of stations and a LDD ("local drain direction" in NetCDF or PCRaster format)
  
> **Note**: PCRaster must be installed in the Conda environment.
 
* __[compare](#compare)__ is a package containing a set of simple Python classes that helps to compare 
NetCDF, PCRaster and TSS files.

* __[thresholds](#thresholds)__ is a tool to compute the discharge return period thresholds from NetCDF4 file containing a discharge time series.

* __[waterregions](#waterregions)__ is a package containing two scripts that allow to create and verify a water regions map, respectively.

* __[gridding](#gridding)__ is a tool to interpolate meteo variables observations stored in text files containing (lon, lat, value) into grids.
  - uses inverse distance interpolation
  - input file names must use format: \<var\>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt
  - option to store all interpolated grids in a single NetCDF4 file
  - option to store each interpolated grid in a GeoTIFF file
  - output files are compressed
  - grids are setup in the configuration folder and are defined by a dem.nc file
  - meteo variables parameters are defined in the same configuration folder

* __[cddmap](#cddmap)__ is a tool to generate correlation decay distance (CDD) maps starting from station timeseries

* __[ncextract](#ncextract)__ is a tool to extract values from NetCDF4 files at specific coordinates.

* __[catchstats](#catchstats)__ calculates catchment statistics (mean, sum, std, min, max...) from NetCDF4 files given masks created with [`cutmaps`](#cutmaps:-a-NetCDF-files-cookie-cutter).

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

Note: if you previously installed an older version of the lisflood-utilitiies, it is highly recommended to remove it before installing the newest version:

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

The tool requires a series of arguments:

* The area to be extracted can be defined in one of the following ways:
    - `-m`, `--mask`: a mask map (either PCRaster or NetCDF format).
    - `-i`, `--cuts_indices`: a bounding box defined by matrix indices in the form `-i imin imax jmin jmax` (the indices must be integers).
    - `-c`, `--cuts`: a bounding box defined by coordinates in the form `-c xmin xmax ymin ymax` (the coordinates can be integer or floating point numbers; x = longitude, y = latitude).
    - `-N`, `-stations`: a list of stations included in a tab separated text file. This approach requires a LDD (local drain direction) map as an extra input, defined with the argument `-l` (`-ldd`).
* The files to be cut may be defined in one of the following ways:
    - `-f`, `--folder`: a folder containing NetCDF files.
    - `-S`, `--static-data`: a directory containint the LISFLOOD static maps. 
* The resulting files will be saved in the folder defined by the argument `-o` ( `--outpath`).

There are additional optional arguments

* `-W`, `--overwrite`: it allows to overwrite results.
* `-C`, `--clonemap`: it can be used to define a clone map when the LDD input map (argument `-l`) is in NetCDF format.

#### Examples 

**Using a mask**

The following command will cut all NetCDF files inside a specific folder (argument `-f`) using a mask (àrgument `-m`). The mask is a boolean map (1 only in the area of interes) that `cutmaps` uses to create a bounding box. The resulting files will be written in the current folder (argument `-o`). 

```bash
cutmaps -m /workarea/Madeira/maps/MaskMap/Bacia_madeira.nc -f /workarea/Madeira/lai/ -o ./
```

**Using indices**

The following command cuts all the maps in an input folder (argument `-f`) using a bounding box defined by matrix indices (argument `-i`). Knowing your area of interest from your NetCDF files, you can determine indices of the array and pass them in the form `-i imin imax jmin jmax` (integer numbers).

```bash
cutmaps -i "150 350 80 180" -f /workarea/Madeira/lai/ -o ./
```

**Using coordinates**

The following command cuts all the static maps in an input folder (argument `-S`) using a bounding box defined by coordinates (argument `-c`). The argument `-W` allows to overwrite pre-existing files in the output directory (argument `-o`):

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

The following command will cut all the static maps in a specific folder (`-S` argument) given a LDD map (`-l` argument) and the previous text file (`-N` argument), and save the results in a folder defined by the argument `-o`.

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.map -N stations.txt -o /Work/Tunisia/cutmaps
```

If the LDD is in NetCDF format, it will be first converted into PCRaster format.

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.nc -N stations.txt -o /Work/Tunisia/cutmaps
``` 

If you experience problems, you can try to pass a path to a PCRaster clone map using the `-C` argument.

```bash
cutmaps -S /home/projects/lisflood-eu -l ldd.nc -C area.map -N stations.txt -o /Work/Tunisia/cutmaps
```

### Output

Apart from the cut files in the output folder specified in the command prompt, `cutmaps` produces other outputs in the folder where the LDD map is stored:

* _mask.map_ and _mask.nc_ for your area of interest, which may be needed in subsequent LISFLOOD/LISVAP executions
* _outlets.map_ and _outlets.nc_ based on _stations.txt_, which will let you produce gauges TSS if configured in LISFLOOD.

## compare

This tool let you compare two NetCDF datasets. You can configure it with tolerances (`atol`, `rtol`, thresholds for percentage of tolerated different values).
You can also set the option to write diff files, so that you can inspect maps and differences with a tool like Panoply.

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
- file *.yaml* or *.json* to define the metadata of the output water regions map in NetCDF format. An example of the structure of these files is provided [here](https://github.com/ec-jrc/lisflood-utilities/tree/master/tests/data/waterregions)

##### Input data provided by this utility:
This utility provides three maps of [Countries IDs](https://github.com/ec-jrc/lisflood-utilities/tree/master/tests/data): 1arcmin map of Europe (EFAS computational domain), 0.1 degree and 3arcmin maps of the Globe . ACKNOWLEDGEMENTS: both the rasters were retrieved by upsampling the original of the World Borders Datase provided by  http://thematicmapping.org/ (the dataset is available under a Creative Commons Attribution-Share Alike License).

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

The `ncextract` tool extracts the time series of values from (multiple) NetCDF file(s) at user defined coordinates.

### Usage

The tool takes as input a CSV file containing point coordinates (structured in 3 columns: id, lat, lon) and a directory containing one or more NetCDF files.

The output is a CSV file (or optionally a NetCDF file) containing the values at the  points corresponding to the provided coordinates, in chronological order.

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
                        Output file (default is CSV, use -nc for NetCDF)
  -nc, --nc             Output to NetCDF
```


## catchstats

The `catchstats` tool calculates catchment statistics given a set of input NetCDF files and a set of mask NetCDF files.

### Usage

The tool takes as input a directory containing the NetCDF files from which the statistics will be computed, and another directory containing the NetCDF files that define the catchment boundaries, which can be any of the outputs of  `cutmaps` (not necessarily the file _my_mask.nc_). The input files can be the LISFLOOD static maps (no temporal dimension) or stacks of maps with a temporal dimension. The mask NetCDF files must be named after with the catchment ID, as this name will be used to identify the catchment in the output NetCDF; for instance: _0142.nc_ would correspond to the mask of catchment 142. Optionally, an extra NetCDF file can be passed to the tool to account for different pixel area; in this case, the statistics will be weighted by this pixel area map.

Only some statistics are currently available: mean, sum, std (standard deviation), var (variance), min, max, median. The weighing based on pixel area does not affect the statistics min, max and median.

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
                          directory containint the input NetCDF files
  -m MASK, --mask MASK
                          directory containing the mask NetCDF files
  -s STATISTIC, --statistic STATISTIC
                          list of statistics to be computed. Possible values: mean, sum, std, var, min, max, median
  -o OUTPUT, --output OUTPUT
                          directory where the output NetCDF files will be saved
  -a AREA, --area AREA
                          NetCDF file of pixel area used to weight the statistics
  -W, --overwrite
                          overwrite existing output files
```

#### Example

The following command calculates the average and total precipitation for a set of catchemtns from the dataset EMO-1. The static map _pixarea.nc_ is used to account for varying pixel area.

```bash
catchstats -i ./EMO1/pr/ -m ./masks/ -s mean sum -o ./areal_precipitation/ -a ./EFAS5/static_maps/pixarea.nc
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


