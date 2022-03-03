# Overview

Scripts to produce csv files with station and discharge data needed for the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). The scripts require a database containing discharge and catchment boundary data in HDF5 format as input. This database cannot be made publicly available but can be requested from the author of this script.

ISIMIP daily meteo forcing data for all scenarios: Download file list using [this search criteria](https://data.isimip.org/search/climate_scenario/ssp119/climate_scenario/ssp126/climate_scenario/ssp245/climate_scenario/ssp370/climate_scenario/ssp460/climate_scenario/ssp534-over/climate_scenario/ssp585/climate_scenario/historical/query//simulation_round/ISIMIP3b/time_step/daily/) and then download with `wget -c -i isimip3b.txt`
 

Download GMTED2010 elevation data from https://www.earthenv.org/topography

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/lisflood-create-obs-discharge-csv
cd lisflood-create-obs-discharge-csv
```
Enter the correct paths and folders in `config.cfg`.

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python step_1_snap_locations.py
python step_2_create_csv.py
```

