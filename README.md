# Overview

Two scripts to produce csv files with station and discharge data needed for the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). A discharge database in HDF5 format is used as input. The location of the database is specified using `database_dir` in `config.ini`. The database cannot be made publicly available but can be requested from the author of this script.

The first script (`step_1_snap_locations.py`) snaps the stations to the 'correct' grid-cell using an automatic procedure based on catchment area and catchment centroid location. If the automatic procedure fails, a window will open to manually select the most appropriate grid-cell. If the window is closed without making a selection, the station will be skipped. The script generates, for each station, a file with the correct station row and column in the folder `corrected_locations_dir` (set using `config.ini`).

The second script (`step_2_create_csv.py`) loads the corrected station location and the discharge data, computes the record length, and generates the csv files for the calibration tool (`Qss.tss` and `stations.csv`). The minimum record length can be set using `config.cfg`. 

A local drainage direction (LDD) map can be produced using https://github.com/hylken/create_ldd_from_merit.

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/lisflood-create-obs-discharge-csv
cd lisflood-create-obs-discharge-csv
```
Modify `config.cfg` and enter the correct paths and folders.

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python step_1_snap_locations.py
python step_2_create_csv.py
```

