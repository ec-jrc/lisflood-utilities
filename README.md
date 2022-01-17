# Overview

This module produces monthly land use/cover fraction maps for 15 SSP-RCP combinations to be used as input to the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model with the `TransientLandUseChange` option enabled. The maps account for yearly changes in the forest, sealed, irrigation (no rice), and irrigated rice fractions. The water fraction represents a monthly climatology and thus does not change yearly. The maps are resampled and subsetted to match the resolution and area of the template map (located at the `templatemap_path` specified in the configuration file). The data are produced for the period spanning `year_start` to `year_end` and saved to `output_folder` (all specified in the configuration file). The raw data sources start in 2015 and 2020 and are provided at 5- and 10-year intervals. The module linearly interpolates between subsequent maps to ensure gradual changes from year to year. For years prior to 2015, the module performs nearest-neighbour extrapolation and simply uses the earliest available data from 2015 and 2020.

# Data

The module requires four datasets (Chen et al., 2020, GCAM-Demeter, GSWE, and HILDA+) as input:
1. [Chen et al. (2020)](https://doi.org/10.1594/PANGAEA.905890) urban land expansion projections under SSP1–5 (1-km resolution; 2020–2100). Each map was reprojected to 0.01° using `gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -r near -te -180.0 -90.0 180.0 90.0 -of GTiff <ssp>/global_<ssp>_<year>.tif <ssp>/global_<ssp>_<year>_reprojected.tif`, where `<ssp>` is SSP1, SSP2, ..., SSP5 and year is 2020, 2030, ..., 2100. The files should be put in `chen_2020_urban_folder`.
1. [GCAM-Demeter](https://dx.doi.org/10.25584/data.2020-04.1190/1615771) land use change projections for 15 SSP-RCP combinations (0.05° resolution; 2015–2100). Downloading all files at once did not work as the download kept getting interrupted. Downloading each SSP separately did work but required a lot of clicking. The files should be put in `gcam_demeter_folder`.
1. [JRC GSWE](https://doi.org/10.1038/nature20584) open water extent data (30-m resolution; 1985–2019). We used a version of the GSWE resampled to 1-km resolution by Susann Guenther available on the JRC network and on request from [Hylke Beck](mailto:hylke.beck@gmail.com). The files should be put in `gswe_folder`.
1. [HILDA+](https://doi.org/10.1594/PANGAEA.921846) historic land use change dataset (1-km resolution; 1900–2019). The files should be put in `hildaplus_folder`.

The locations of the datasets are specified in the configuration file. The module also requires a template map which defines the output resolution and area. The template map should be in netCDF-4 format and contain `lat` and `lon` variables and a data variable (any name). The location of the template map is specified using `templatemap_path` in the configuration file.

# Methods

The module consists of four scripts: 
1. `step1_Chen_2020_urban.py`: Resamples the reprojected Chen et al. (2020) sealed fraction projections under SSP1–5 to the template resolution.
1. `step2_GCAM_Demeter.py`: Resamples the GCAM-Demeter forest, irrigation, and rice fraction projections for 15 SSP and RCP combinations to the template resolution.
1. `step3_JRC_GSWE.py`: Resamples the JRC GSWE water fraction to the template resolution and computes a monthly climatology. HILDA+ data are used at high latitudes.
1. `step4_harmonization.py`: Loops through the 15 SSP-RCP combinations, years, and months and harmonizes the different data sources such that the sum of all fractions equals one. For each loop, the script:
	1. loads the Chen et al. (2020) fracsealed data (interpolating between years);
	1. loads the GSWE fracwater data (monthly climatology);
	1. reduces fracwater if fracsealed+fracwater >1 (fracsealed overrides fracwater);
	1. computes fracother as residual of fracsealed and fracwater;
	1. loads GCAM-Demeter fracforest, fracirrigation, and fracrice data (interpolating between years);
	1. reduces fracforest, fracirrigation, and fracrice if sum exceeds fracother;
	1. recalculates fracother as the residual; and
	1. subsets all data to the template map area.	
	Finally, the fractions are converted to netCDF-4 format (one file for each RCP-SSP combination and fraction).

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