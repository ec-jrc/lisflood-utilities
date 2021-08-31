# Overview

These scripts produce csv files to be used as input to the LISFLOOD calibration tool (https://github.com/ec-jrc/lisflood-calibration). The first script (`step_1_snap_locations.py`) aligns the stations with the river network of the LDD () using an automatic proximity-search procedure and a manual point-and-click procedure as fallback. scripts requireThe script is entirely self-containing and only requires one to specify the path of the clone map and several output folders.

The script:
1. downloads and extracts the MERIT Hydro upstream area data;
1. resamples the data using a maximum filter (to retain the river network);
1. uses the resampled upstream area (inverted) as elevation; and
1. computes the LDD.

# System requirements

The script can be run on a normal desktop PC (Windows and Linux). The download of the MERIT data requires ~45 GB of free disk space. The download will be performed only once.

# Instructions

Register at the [MERIT Hydro](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/) website using the Google form (important).

If you're using Windows, install the cross-platform [wget](http://gnuwin32.sourceforge.net/packages/wget.htm) and [tar](http://gnuwin32.sourceforge.net/packages/gtar.htm) utilities to download and extract the MERIT Hydro data. Add the location of the utilities to the environment variables, so they can be executed from any location.

Clone the repository:
```
git clone https://github.com/hylken/create_ldd_from_merit
cd create_ldd_from_merit
```
Modify `config.cfg` with the correct paths and folders. `res` should match the resolution of the clone map. The clone map should contain `lat` and `lon` variables and a data variable (any name).

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python create_ldd_from_merit.py
```

