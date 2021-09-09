# Overview

This script produces monthly dynamic land use/cover fraction maps to be used as input to the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model with the `TransientLandUseChange` option enabled. The maps are resampled and subsetted to match resolution and area of the clone map (located at the `clonemap_path` specified in the configuration file). 

# Data sources

The script uses three data sources (HILDA+, HYDE, and GSWE):

1. The `fracforest` and `fracsealed` are based on HILDA+ V1.0 ([Winkler et al., 2021](https://doi.org/10.1038/s41467-021-22702-2)). Download the file `hildap_vGLOB-1.0-f_netcdf.zip` from [the PANGAEA data repository](https://doi.org/10.1594/PANGAEA.921846) and extract it to the `hildaplus_folder` specified in the configuration file.
1. The `fracirrigation` and `fracrice` are based on HYDE V3.2 ([Klein Goldewijk et al., 2017](https://doi.org/10.5194/essd-9-927-2017)). Download the `baseline` and `general_files` folders from the [DANS data portal](https://doi.org/10.17026/dans-25g-gez3) and put them in the `hyde_folder` specified in the configuration file.
1. The `fracwater` is based on GSWE V4 ([Pekel et al., 2016](https://doi.org/10.1038/nature20584)).

# Overview

The script carries out the following steps to obtain the fractions for each year based on HILDA+ and HYDE:
1. Loads and resamples the global 1-km HILDA+ data to the clone map resolution and calculates the water, forest, and sealed fractions.
1. Loads and resamples the global 0.083° HYDE data to the clone map resolution.
1. Fixes HYDE fractions exceeding 1.
1. Makes sure the HYDE cropland fraction does not exceed the HILDA+ 'other' fraction and calculates the rice and irrigation (no rice) fractions.
1. Makes sure the sum of the five non-other fractions is <1 and calculates the other fraction.
1. The water fraction will be replaced with GSWE data and the other five fractions will be rescaled accordingly. However, if the the water fraction is 1, the other fractions cannot be adjusted, as they will all be 0. As a workaround, we reduce the initial water fraction by a tiny amount while increasing the other fractions by a tiny amount using interpolated (non-zero) values.
1. Makes sure the sum of the fractions is 1.

The scripts carries out the following steps to adjust the yearly fractions using GSWE on a monthly basis:
1. Loads and resamples the global 1-km GSWE data to the clone map resolution.
1. Updates the water fraction using the GSWE data.
1. Adjusts the other five fractions according to the updated water fraction.
1. Makes sure the fractions sum to 1.
1. Subsets the global maps to the clone map area.
1. Saves the data to netCDF-4 files (one for each fraction).

# System requirements

The script can be run on a normal desktop PC (Windows and Linux).

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

