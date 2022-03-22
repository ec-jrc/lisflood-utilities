# Overview

Module to generate daily meteorological input files for the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model. The module consists of several similar scripts tailored to different input meteorological datasets, such as ISIMIP3b, W5E5, and MSWX/MSWEP. The following variables are created: precipitation (`pr.nc`), 2-m air temperature (`ta.nc`), reference potential evaporation (`et.nc`), open water potential evaporation (`ew.nc`), and bare soil potential evaporation (`es.nc`).

The following should be kept in mind when using the scripts:
1. The scripts calculate potential evaporation directly from the input data following the [LISVAP](https://github.com/ec-jrc/lisflood-lisvap) approach. Separately running the LISVAP module is, therefore, no longer necessary. 
2. The data are resampled and subsetted to the resolution and area of the template map (located at `templatemap_path` specified in the configuration file).
3. Air temperature and air pressure are downscaled to the template map resolution (up to 1 km using a simple delta lapse-rate correction).
4. The output is written to the `scratch_folder` and moved to the `output_folder` once the processing is done. The `scratch_folder` should point to a storage location dedicated to a lot of file accesses.
5. The ISIMIP3b and W5E5 scripts load the input data into memory (using the `diskless=True` argument) to avoid read errors which frequently occurred on the system used for developing the scripts.
6. The ISIMIP3b script can be run simultaneously multiple times to simultaneously process multiple scenarios/models. 

# Data

The following datasets are required depending on the script in question:

1. GMTED2010 1-km mean surface elevation data can be downloaded from the [EarthEnv website](
https://data.earthenv.org/topography/elevation_1KMmn_GMTEDmn.tif) (required for all scripts). Put `elevation_1KMmn_GMTEDmn.tif` in `gmted2010_folder` (specified in the configuration file).

2. ISIMIP3b daily meteorological data (historical and climate projections). Download [this](https://data.isimip.org/api/v1/datasets/filelist/?page=1&climate_scenario=ssp119&climate_scenario=ssp126&climate_scenario=ssp245&climate_scenario=ssp370&climate_scenario=ssp460&climate_scenario=ssp534-over&climate_scenario=ssp585&climate_scenario=historical&query=&ISIMIP3b=time_step&simulation_round=ISIMIP3b&time_step=daily) file list and download the data with `wget -c -i file-list.txt`. Put the data in `isimip3b_folder`.

3. [MSWX](http://www.gloh2o.org/mswx) and [MSWEP](http://www.gloh2o.org/mswep) daily historical meteorological data. Follow the download instructions on the respective web pages. Use [this](rclone_filter_file.txt) filter file for rclone. Put the data in `mswx_folder` and `mswep_folder`, respectively. Retain the folder structure (e.g., `<mswx_folder>/Past/Temp/Daily/2007133.nc`).

4. GSWP3-W5E5 daily historical meteorological data (`obsclim` and `counterclim`). Download [this](https://data.isimip.org/api/v1/datasets/filelist/?page=1&tree=ISIMIP3a&InputData=climate&atmosphere=gswp3-w5e5&climate_scenario=counterclim&climate_scenario=obsclim&time_step=daily&climate_forcing=gswp3-w5e5) file list and download the data with `wget -c -i file-list.txt`. Put the data in `w5e5_folder`.

# System requirements

The scripts can be run on a normal desktop PC (Windows and Linux) with 32 GB or more of physical memory.

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
