# Overview

Scripts to generate daily meteorological input files for the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model. The following variables are created:
1. precipitation (`pr.nc`), 
2. temperature (`ta.nc`), 
3. reference potential evaporation (`et.nc`), 
4. open water potential evaporation (`ew.nc`), and 
5. soil potential evaporation (`es.nc`). 

For efficiency, the scripts calclate potential evaporation independently and thus the extra task of running the [LISVAP](https://github.com/ec-jrc/lisflood-lisvap) module is no longer required. The scripts downscale air temperature and air pressure to the output resolution is selected (up to 1 km using a simple delta lapse rate correction). The maps are resampled and subsetted to match the resolution and area of the template map (located at the `templatemap_path` specified in the configuration file). 

# Data


1. ISIMIP daily meteo forcing data for historical period and all climate change scenarios. Download [this](https://data.isimip.org/api/v1/datasets/filelist/?page=1&climate_scenario=ssp119&climate_scenario=ssp126&climate_scenario=ssp245&climate_scenario=ssp370&climate_scenario=ssp460&climate_scenario=ssp534-over&climate_scenario=ssp585&climate_scenario=historical&query=&ISIMIP3b=time_step&simulation_round=ISIMIP3b&time_step=daily) file list and download the data files with `wget -c -i isimip3b.txt`.
2. MSWX
3. W5E5
4. GMTED2010 1-km mean surface elevation data from the [EarthEnv website](
https://data.earthenv.org/topography/elevation_1KMmn_GMTEDmn.tif). Put the file in the `gmted2010_folder` specified in the configuration file.

 

Download GMTED2010 elevation data from https://www.earthenv.org/topography (dataset:Elevation, Aggregation: Median

# Instructions

conda create -n <env> -c conda-forge pcraster h5py pandas numpy netcdf4 matplotlib rasterio scikit-image


