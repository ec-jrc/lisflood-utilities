# Overview

Scripts to produce csv files with station and discharge data needed for the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). The scripts require a database containing discharge and catchment boundary data in HDF5 format as input. This database cannot be made publicly available but can be requested from the author of this script.

ISIMIP daily meteo forcing data for all scenarios: Download file list using [this search criteria](https://data.isimip.org/search/climate_scenario/ssp119/climate_scenario/ssp126/climate_scenario/ssp245/climate_scenario/ssp370/climate_scenario/ssp460/climate_scenario/ssp534-over/climate_scenario/ssp585/climate_scenario/historical/query//simulation_round/ISIMIP3b/time_step/daily/) and then download with `wget -c -i isimip3b.txt`
 

Download GMTED2010 elevation data from https://www.earthenv.org/topography (dataset:Elevation, Aggregation: Median

# Instructions

conda create -n <env> -c conda-forge pcraster h5py pandas numpy netcdf4 matplotlib rasterio scikit-image


