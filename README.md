# Overview

This script produces monthly future land use/cover fraction maps to be used as input to the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model with the `TransientLandUseChange` option enabled. The maps account for yearly changes in the forest, sealed, irrigation (no rice), and irrigated rice fractions, and monthly changes in the water fraction (with the other fractions adjusted accordingly). The maps are resampled and subsetted to match the resolution and area of the template map (located at the `templatemap_path` specified in the configuration file). 

# Data

The script is based on five data sources (GAIA, GSWE, HILDA+, MEaSUREs VCF, and HYDE):
1. The sealed fraction was derived form the GAIA impervious area dataset (30-m resolution; 1985–2018; [Gong et al., 2020](https://doi.org/10.1016/j.rse.2019.111510)). We did not use data prior to 1990, due to a lower reliability (e.g., Monterrey in northeastern Mexico is almost non-existant). The script to download and resample GAIA data to 1-km resolution is included (`GAIA_download_and_process.py`).
1. The water fraction was derived from GSWE V4 (30-m resolution; 1985–2019; [Pekel et al., 2016](https://doi.org/10.1038/nature20584)). We used a version of the GSWE resampled to 1-km resolution by Susann Guenther available on the JRC network. The code used to generate this version from the raw data is also included [waiting for a response from Susann].
1. The water fraction at high latitudes was based on HILDA+ V1.0 (1-km resolution; 1900–2019; [Winkler et al., 2021](https://doi.org/10.1038/s41467-021-22702-2)) available via the [PANGAEA data repository](https://doi.org/10.1594/PANGAEA.921846).
1. The forest fraction was derived from the MEaSUREs VCF dataset (0.05° resolution; 1982–2016; [Song et al., 2018](https://doi.org/10.1038/s41586-018-0411-9)) available from the [NASA Earthdata website](https://doi.org/10.5067/MEaSUREs/VCF/VCF5KYR.001).
1. The irrigation (no rice) and irrigated rice fractions were based on HYDE V3.2 (0.083° resolution; 10,000 BC to 2015; [Klein Goldewijk et al., 2017](https://doi.org/10.5194/essd-9-927-2017)). The `baseline` and `general_files` folders should be downloaded from the [DANS data portal](https://doi.org/10.17026/dans-25g-gez3).

The locations of the data are specified in the configuration file. The script also requires a template map which defines the output resolution and area. The template map should be in netCDF-4 format and contain `lat` and `lon` variables and a data variable (any name). The location of the template map is specified using `templatemap_path` in the configuration file.

# Methods

The script carries out the following steps for each month:
1. Resamples the global 1-km GAIA data to the template map resolution and calculates the sealed fraction.
1. Resamples the global 1-km GSWE data to the template map resolution and calculates the water fraction.
1. Resamples the global 1-km HILDA+ data to the template map resolution and fills gaps in the GSWE water fraction map.
1. Makes sure the sealed and water fractions sum to <=1.
1. Computes the other fraction by subtracting the sealed and water fractions from 1.
1. Resamples the global 0.05° MEaSUREs VCF data and calculates the forest fraction.
1. Resamples the global 0.083° HYDE data and calculates the irrigation (no rice) and irrigated rice fractions.
1. Reduces the forest, irrigation (no rice), and irrigated rice fractions if their sum >1.
1. Makes sure the sum of all fractions (excluding other) is <=1.
1. Recalculates the other fraction as the residual.
1. Subsets the global maps to the template map area.
1. Saves the data to compressed Numpy files (one for each fraction each month).

Finally, the fractions are converted to netCDF-4 format (one file for each fraction). The data are produced for the period spanning `year_start` to `year_end` and saved to `output_folder` (all specified in the configuration file). For years without data, the script uses data from the closest year. For 1984, for example, HYDE data from 1980 (the closest year with HYDE data) are used, while for March 1981, GSWE data from March 1985 (the first year with GSWE data) are used.

# System requirements

The script can be run on a normal desktop PC (Windows and Linux) with 16 GB or more of physical memory.

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/lisflood-land-cover-future
cd lisflood-land-cover-future
```
Produce a configuration file with the correct paths and folders based on the provided template. 

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python main.py <config file>
```
If the environment creation step fails, we recommend creating the environment and installing the packages as follows:
```
conda create -n <env> -c conda-forge geopandas h5py pandas numpy netcdf4 matplotlib rasterio scikit-image
```