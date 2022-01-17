# Overview

This module produces monthly future land use/cover fraction maps to be used as input to the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model with the `TransientLandUseChange` option enabled. The maps account for yearly changes in the forest, sealed, irrigation (no rice), and irrigated rice fractions. The water fraction represents a monthly climatology and thus does not change yearly. The maps are resampled and subsetted to match the resolution and area of the template map (located at the `templatemap_path` specified in the configuration file). 

# Data

The module requires four datasets (Chen et al., 2020, GCAM-Demeter, GSWE, and HILDA+):
1. [Chen et al. (2020)](https://doi.org/10.1594/PANGAEA.905890) urban land expansion projections under SSP1–5 (1-km resolution; 2020–2100). Each map was reprojected to 0.01° using `gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -r near -te -180.0 -90.0 180.0 90.0 -of GTiff <ssp>/global_<ssp>_<year>.tif <ssp>/global_<ssp>_<year>_reprojected.tif`, where `<ssp>` is SSP1, SSP2, ..., SSP5 and year is 2020, 2030, ..., 2100.
1. [GCAM-Demeter](https://dx.doi.org/10.25584/data.2020-04.1190/1615771) land use change projections for 15 SSP and RCP combinations (0.05° resolution; 2015–2100). Downloading all files at once did not work as the download kept getting interrupted. Downloading each SSP separately did work but required a lot of clicking.
1. [JRC GSWE](https://doi.org/10.1038/nature20584) open water extent data (30-m resolution; 1985–2019). We used a version of the GSWE resampled to 1-km resolution by Susann Guenther available on the JRC network and on request from [Hylke Beck](mailto:hylke.beck@gmail.com).
1. [HILDA+](https://doi.org/10.1594/PANGAEA.921846) historic land use change dataset (1-km resolution; 1900–2019).

The locations of the data are specified in the configuration file. The script also requires a template map which defines the output resolution and area. The template map should be in netCDF-4 format and contain `lat` and `lon` variables and a data variable (any name). The location of the template map is specified using `templatemap_path` in the configuration file.

# Methods

The module consists of four scripts: 
1. `step1_Chen_2020_urban.py`: Resamples the reprojected Chen et al. (2020) urban land expansion projections under SSP1–5 to the template resolution.
1. `step2_GCAM_Demeter.py`: Resamples the GCAM-Demeter land use change projections for 15 SSP and RCP combinations to the template resolution.
1. `step3_JRC_GSWE.py`: Resamples the JRC GSWE open water extent data to the template resolution and computes a climatology. HILDA+ data are used at high latitudes.
1. `step4_harmonization.py`: 
	Carries out the following steps for each month:
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