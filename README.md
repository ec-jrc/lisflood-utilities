# Overview

These scripts produce csv files to be used as input to the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). 
The first script (`step_1_snap_locations.py`) snaps the stations to the 'correct' grid-cell using an automatic proximity search procedure based on catchment area and catchment centroid location. If the procedure fails, a window will open to manually select the most appropriate grid-cell. If the window is closed without making a selection, the station will be skipped. This script generates, for each station, a file with the correct grid-cell row and column in the folder `corrected_locations_dir` (set using `config.ini`).

The second script (`step_2_create_csv.py`) selects stations with a sufficiently long record and produces the csv files (`Qss.tss` and `stations.csv`). The minimum record length can be set using `config.cfg`.

A local drainage direction map for the study area can be produced using https://github.com/hylken/create_ldd_from_merit.

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

