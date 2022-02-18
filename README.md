# Overview

Scripts to produce csv files with station and discharge data needed for the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). The scripts require a database containing discharge and catchment boundary data in HDF5 format as input. This database cannot be made publicly available but can be requested from the author of this script.

The first script (`step_1_snap_locations.py`) snaps the stations to the 'correct' grid-cell using an automatic procedure based on catchment area and catchment centroid location. If the automatic procedure fails, a window will open to manually select the most appropriate grid-cell. The corrected station locations are stored in the folder `corrected_locations_dir` (set using `config.ini`). If the window is closed without making a selection, the station will be omitted.

The second script (`step_2_create_csv.py`) loads the corrected station location and discharge data and computes the record length for each station, and generates the csv files for the calibration tool (`Qtss.csv` and `stations.csv`). The minimum record length can be set using `config.cfg`. 

Local drainage direction (LDD) and upstream area maps can be produced using https://github.com/hylken/create_ldd_from_merit.

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

