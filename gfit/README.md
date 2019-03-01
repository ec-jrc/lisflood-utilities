# GFIT

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
