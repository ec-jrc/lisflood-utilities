# Overview

This module produces monthly historic sectoral water withdrawal maps to be used as input to the [LISFLOOD](https://github.com/ec-jrc/lisflood-code) hydrological model. The maps cover four water use sectors: domestic, thermoelectric, manufacturing, and livestock. The maps are resampled and subsetted to match the resolution and area of the template map (located at `templatemap_path` specified in the configuration file). 

# Data

The module consists of five scripts which require the following datasets and files:
1. [FAO AQUASTAT](http://fao.org/aquastat/statistics/query/index.html?lang=en) sectoral water withdrawal estimates. Select "All Countries" and the following seven fields: "Gross Domestic Product (GDP)", "Industry, value added to GDP", "Agricultural water withdrawal", "Industrial water withdrawal", "Municipal water withdrawal", "Total water withdrawal", and "Irrigation water withdrawal". Export as "CSV (Flat)". The generated AQUASTAT CSV is, unfortunately, not readily usable. First, the CSV contains a long list of footnotes at the end, which should be manually deleted. Second, the header does not include the last column, which should be fixed as well (e.g., by appending ,"Unknown" to the header line). The fixed CSV should be named `aquastat_clean.csv` and put in `aquastat_folder` specified in the configuration file.
1. [Global Human Settlement Layer (GHSL)](https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php) POP R2019A residential population estimates for target years 1975, 1990, 2000, and 2015. For each year, download the 0.0025° (9 arcsec) maps in WGS84 projection and resample them to 0.01° using `gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -r average -of GTiff GHS_POP_E<year>_GLOBE_R2019A_4326_9ss_V1_0.tif GHS_POP_E<year>_GLOBE_R2019A_4326_9ss_V1_0_reprojected.tif`, where `<year>` is 1975, 1990, 2000, or 2015. We resample the maps using `average` instead of `sum` as the latter is broken.
1. [Global Change Analysis Model (GCAM)](https://github.com/JGCRI/gcam-core/releases) regional water withdrawal and electricity consumption estimates. Download the Windows release package, replace `configuration.xml` by `configuration_ssp.xml` in `run-gcam.bat`, execute `run-gcam.bat`, wait for the model to finish. Execute `run-model-interface.bat`, click "File" > "Open" > "DB Open", select `output/database_basexdb`, and select all scenarios and all regions. Select `water demand`, select `water withdrawals by sector`, and click "Run query". Select `energy transformation`, select `electricity`, select `elec consumption by demand sector`, and click "Run query". For each tab, select all data with ctrl-a and click "Edit" > "Copy". Open a blank [Google Sheets](http://sheets.google.com/) spreadsheet, press ctrl-v, manually add the headers, click "File" > "Download" > "Comma-separated values", save as `water_withdrawals.csv` and `elec_consumption.csv`, respectively, and put both files in `gcam_folder`.
1. [Gridded Livestock of the World (GLW 3)](https://doi.org/10.1038/sdata.2018.227) species distribution dataset. Download the eight zip files (each representing a different species) and extract them to `glw_folder`.
1. [Huang et al. (2018)](https://doi.org/10.5281/zenodo.1209296) global gridded water withdrawal estimates. These estimates are not incorporated in our dataset, but are only used for the sake of comparison. Download all 7z files and extract them to `huang_folder`.
1. [Multi-Source Weather (MSWX)](http://www.gloh2o.org/mswx) daily and monthly mean air temperature. Download using rclone as explained in the FAQ on the web page. Put the daily and monthly netCDFs in `<mswx_folder>/Past/Temp/Daily` and `<mswx_folder>/Past/Temp/Monthly`, respectively, where `mswx_folder` is specified in the configuration file.
1. [US Census Bureau](https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html) state borders shape file. Download `cb_2018_us_state_500k.zip`, unzip it, open the shape file in QGIS, open the Field Calculator, enter `STATEFP` in "Expression", enter `STATEFP_` in "Output field name", and click OK. Rasterize to 0.01° using `gdal_rasterize -l cb_2018_us_state_500k -a STATEFP_ -tr 0.01 0.01 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff cb_2018_us_state_500k.shp cb_2018_us_state_500k_rasterized.tif`. Put `cb_2018_us_state_500k_rasterized.tif` in `us_states_folder`.
1. [USGS NWIS](https://waterdata.usgs.gov/nv/nwis/wu) water withdrawal estimates for 1985--present. For each state in the "Geographic Area" drop-down menu, select "State Data", "ALL Years", "State Total", and "ALL Categories" and click "Submit". Then select "Tab-separated data" and click "Submit". Do this for each state and put files in `usgs_water_use_folder`.
1. USGS water withdrawal data files for [1985](https://water.usgs.gov/watuse/data/1985/index.html) and [1990](https://water.usgs.gov/watuse/data/1990/index.html) (to supplement the NWIS data). Download "Data file of state-level data" for each year (`us85st.txt` and `us90st.txt`) and put the files in `usgs_water_use_folder`.
1. [Vassolo and Döll (2005)](https://doi.org/10.1029/2004WR003360) industrial and thermoelectric water withdrawal maps requested from Petra Döll. The industrial withdrawal shape file is rasterized using `gdal_rasterize -l "industry_paper_vassolo&doell" -a MANUF_WWD -tr 0.5 0.5 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff "industry_paper_vassolo&doell.shp" manuf_wwd.tif`. Thermoelectric shape file rasterized using same command but with `WWD_PS` and `wwd_ps.tif`. Put the files in `vassolo_doll_folder`.
1. [Thematic Mapping](https://thematicmapping.org/downloads/world_borders.php) country borders shape file. Download `TM_WORLD_BORDERS-0.3.zip`, unzip, and rasterize to 0.01° using `gdal_rasterize -l TM_WORLD_BORDERS-0.3 -a UN -tr 0.01 0.01 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff TM_WORLD_BORDERS-0.3.shp TM_WORLD_BORDERS_UN_rasterized_.tif`. Put the result in `world_borders_folder`.
1. [World Bank](https://data.worldbank.org/) manufacturing value added and gross domestic product data. Search for "Manufacturing, value added (constant 2010 US$)" and "GDP (constant 2010 US$)", download as CSV, and put in `world_bank_folder`.
<!---
1. USGS Circular xxx 1980 (https://pubs.usgs.gov/circ/1983/1001/report.pdf) Table 7 provided in ancillary_data folder
1. [Global Power Plant Database](https://datasets.wri.org/dataset/globalpowerplantdatabase)
1. [Lohrmann et al. (2019) power plant database](https://doi.org/10.1038/s41560-019-0501-4) Supplementary Data 1
1. [UN Industrial Commodity Statistics Database](http://data.un.org/), enter for example <beer>, then under "Select filters" select "Industry" and click "Apply filter", click "Download", click "Comma", extract and rename.
1. EIA [electricy capacity](https://www.eia.gov/international/data/world/electricity/), Select World, select Electricity, select Capacity and Generation, click View Data, click Sort by Energy Source/Activity, click Download Options, select Export CSV (table).
-->

The locations of the datasets and files are specified in the configuration file (see the `config_GloFAS_2021.cfg` example). The script also requires a template map which defines the output resolution and area. The template map should be in netCDF-4 format and contain `lat` and `lon` variables and a data variable (any name). The location of the template map is specified using `templatemap_path` in the configuration file. The data are produced for the period spanning `year_start` to `year_end` and saved in netCDF-4 format to `output_folder` (all specified in the configuration file). 

# Methods

The module consists of five scripts:
1. `step1_population_density.py`: 
	1. Resamples the 1-km GHSL population grids to the template map resolution.
	1. Computes population grids (in density per km<sup>2</sup>) for each year using linear interpolation and nearest-neighbor extrapolation.
1. `step2_domestic_demand.py`:
	1. Produces a map of the R parameter (measures the relative withdrawal difference between the warmest and coldest months).
	1. Computes annual time series of population for each country.
	1. Loads the country-scale AQUASTAT industrial withdrawal estimates and produces annual time series using linear interpolation and population-based forward and backward extrapolation.
	1. Loads the state-level USGS withdrawal data and computes annual values using linear interpolation and nearest-neighbour extrapolation.
	1. For each year, computes country- and state-scale per capita water demand and produces a global map.
	1. Gap-fills the annual per capita water demand map using nearest-neighbor interpolation. Necessary for countries without any AQUASTAT estimates, such as Sudan.
	1. Spatially disaggregates the country- and state-scale annual per capita water demand estimates using the population estimates to obtain global annual withdrawal grids.
	1. Temporally downscales the annual withdrawal grids using monthly MSWX air temperature grids and the map of the R parameter.
1. `step3a_industrial_demand.py`: 
	1. Computes country-scale values of the Vassolo and Döll (2005) industrial and thermoelectric withdrawal grids.
	1. Processes the country-scale AQUASTAT industrial withdrawal data.
	1. Processes the country-scale World Bank manufacturing value added (MVA) data.
	1. Processes the GCAM regional industry and thermoelectric withdrawals and spatially downscales them to the country scale based on population.
	1. Processes the GCAM regional electricity consumptions and assigns values to each country.
	1. Computes annual country-scale manufacturing and thermoelectric withdrawal time series as follows:
		1. If AQUASTAT and GCAM data are available, temporally gap-fills the AQUASTAT industrial withdrawals using linear interpolation and MVA- or population-based extrapolation, temporally gap-fills the GCAM manufacturing and thermoelectric withdrawals using linear interpolation and nearest-neighbor extrapolation, rescales the GCAM manufacturing and thermoelectric withdrawals to match the Vassolo and Döll (2005) estimates, and temporally disaggregates the AQUASTAT industrial withdrawals into manufacturing and thermoelectric using GCAM.
		1. If only GCAM data are available, simply use the rescaled GCAM manufacturing and thermoelectric withdrawals.
		1. If neither AQUASTAT nor GCAM data are available, set the manufacturing and thermoelectric withdrawals to zero.
	1. Computes monthly heating degree days (HDD) and cooling degree days (CDD) grids using daily MSWX air temperature.
1. `step3b_industrial_demand.py`:
	1. Loads the annual country-scale manufacturing and thermoelectric withdrawal time series, the USGS manufacturing and thermoelectric withdrawal estimates, and the GCAM electricity consumptions from the preceding scripts.
	1. Gap-fills the annual USGS manufacturing and thermoelectric withdrawal time series using linear interpolation and nearest-neighbor extrapolation.
	1. Spatially downscales the country- and state-scale manufacturing and thermoelectric withdrawals using population grids to obtain annual withdrawal grids.
	1. Temporally downscales the annual thermoelectric withdrawal grids based on the GCAM electricity consumptions and the HDD and CDD grids.
1. `step4_livestock_demand.py`:
	1. Computes country-scale livestock withdrawals by taking the difference between the AQUASTAT agriculture and irrigation withdrawals.
	1. Processes the raw GLW grids for the different species and estimates the total livestock mass.
	1. Computes the total livestock mass for each country.
	1. Spatially downscales the GCAM regional livestock withdrawals to the country scale using the GLW-based total livestock mass estimates.
	1. Loops through countries and gap-fills the annual GCAM time series using linear interpolation and nearest-neighbor extrapolation, and rescales the time series to match the AQUASTAT estimate (if available).
	1. Gap-fills the annual USGS manufacturing and thermoelectric withdrawal time series using linear interpolation and nearest-neighbor extrapolation.
	1. For each year, spatially downscales the country- and state-scale livestock withdrawals using the GLW total livestock mass grid.

# System requirements

The script can be run on a normal desktop PC (Windows and Linux) with 16 GB or more of physical memory.

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/lisflood-water-demand-historic
cd lisflood-water-demand-historic
```
Produce a configuration file with the correct paths and folders (based on the included example). 

Create and activate a Conda environment and run the script as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python main.py <config file>
```
If the environment creation step fails, the following might work:
```
conda create -n <env> -c conda-forge geopandas h5py pandas numpy netcdf4 matplotlib rasterio scikit-image openpyxl geopy
```