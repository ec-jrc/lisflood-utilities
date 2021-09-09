# Overview

This script produces monthly dynamic land use/cover fraction maps to be used as input to the LISFLOOD model with the `Transient land use change option` option enabled. The maps are resampled and subsetted to match resolution and area of the clone map (the location of which is specified in the configuration file). 

# Data sources

1. The `fracforest` and `fracsealed' are based on [HILDA+](https://doi.org/10.1038/s41467-021-22702-2) V1.0. Download the file `hildap_vGLOB-1.0-f_netcdf.zip` from [the PANGAEA data repository](https://doi.org/10.1594/PANGAEA.921846) and extract it to the `hildaplus_folder` specified in the configuration file.
1. The `fracirrigation` and `fracrice` are based on [HYDE](https://doi.org/10.5194/essd-9-927-2017) V3.2. Download the `baseline' and `general_files' folders from the [DANS data portal](https://doi.org/10.17026/dans-25g-gez3) and put them in the `hyde_folder` specified in the configuration file.
1. The `fracwater` is based on [GSWE](https://doi.org/10.1038/nature20584) V4.


# Overview


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

