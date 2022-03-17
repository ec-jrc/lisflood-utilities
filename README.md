# Overview

Module to generate daily meteorological input files for the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model. The following variables are created: 1) precipitation (`pr.nc`), 2) temperature (`ta.nc`), 3) reference potential evaporation (`et.nc`), 4) open water potential evaporation (`ew.nc`), and 5) bare soil potential evaporation (`es.nc`). 

The module consists of several similar scripts tailored to different input meteorological datasets, such as ISIMIP3b, W5E5, and MSWEP. The following should be kept in mind when using the scripts:
1. The scripts calculate potential evaporation directly from the input data following the [LISVAP](https://github.com/ec-jrc/lisflood-lisvap) approach. Separately running the LISVAP module is, therefore, no longer necessary. 
2. The data are resampled and subsetted to the resolution and area of the template map (located at `templatemap_path` specified in the configuration file).
3. Air temperature and air pressure are downscaled to the template map resolution (up to 1 km using a simple delta lapse rate correction). 
4. The output is written to the `scratch_folder` and moved to the `output_folder` once the processing is done. The `scratch_folder` should point to a storage location dedicated to a lot of file accesses.
5. The ISIMIP3b and W5E5 scripts load the input data into memory (using the `diskless=True` argument) to avoid read errors which frequently occurred on the system used for developing the scripts.

# Data

All scripts require GMTED2010 1-km mean surface elevation data from the [EarthEnv website](
https://data.earthenv.org/topography/elevation_1KMmn_GMTEDmn.tif). Put `elevation_1KMmn_GMTEDmn.tif` in the `gmted2010_folder` specified in the configuration file.

## main_ISIMIP3b_projections.py

Requires ISIMIP daily meteorological forcing data for the historical period and all climate change scenarios. Download [this](https://data.isimip.org/api/v1/datasets/filelist/?page=1&climate_scenario=ssp119&climate_scenario=ssp126&climate_scenario=ssp245&climate_scenario=ssp370&climate_scenario=ssp460&climate_scenario=ssp534-over&climate_scenario=ssp585&climate_scenario=historical&query=&ISIMIP3b=time_step&simulation_round=ISIMIP3b&time_step=daily) file list and download the data with `wget -c -i isimip3b.txt`. Put the data in `isimip3b_folder`.


## main_MSWX_MSWEP_reanalysis.py

Requires [MSWX](www.gloh2o.org/mswx) and [MSWEP](www.gloh2o.org/mswep) daily meteorological data. Follow the download instructions on the respective web pages. Use the following filter file for rclone:
```
+ MSWX_V100/Past/**/Daily/*.nc
+ MSWX_V100/Past/**/Monthly/*.nc
+ MSWEP_V280/Past/Daily/*.nc
+ MSWEP_V280/Past/Monthly/*.nc
+ MSWEP_V280/NRT/Daily/*.nc
+ MSWEP_V280/NRT/Monthly/*.nc
- *
```
Put the data in `mswx_folder` and `mswep_folder`, respectively. 

## main_W5E5_reanalysis.py

Requires W5E5

# System requirements

The script can be run on a normal desktop PC (Windows and Linux) with 16 GB or more of physical memory.

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/lisflood-meteo-forcing
cd lisflood-meteo-forcing
```
Produce a configuration file with the correct paths and folders based on the provided template.

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python main_ISIMIP3b_projections.py <config file>
```
If the environment creation step fails, we recommend creating the environment and installing the packages as follows:
```
conda create -n <env> -c conda-forge scipy pandas numpy netcdf4 matplotlib rasterio scikit-image
```
