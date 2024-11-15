# Overview

This utility allows to create sectoral water demand maps at the desired spatial resolution and for the desired geographical areas. The maps indicate, for each pixel, the time-varying water demand map to supply for domestic, livestock, industrial, and energetic (cooling) water consumption. The temporal discretization of the information is monthly for domestic and energy demand, yearly for industrial and livestock demand. The users can select the historical temporal span of the output datasets (start year and end year). The maps of sectoral water demand are required by the LISFLOOD OS [water use module](https://ec-jrc.github.io/lisflood-model/2_18_stdLISFLOOD_water-use). Clearly, the scripts and the outputs of this utility can be used also for other applications, as well as for stand-alone analysis of historical water demand for anthropogenic use. Due to the large uncertainty of the global data sets, the distinction between demand and withdrawal is here omitted.

This README file provides detailed technical information about the input datasets and the usage of this utility. The methodology is explained in the manuscript: Choulga, M., Moschini, F., Mazzetti, C., Grimaldi, S., Disperati, J., Beck, H., Salamon, P., and Prudhomme, C.: Technical note: Surface fields for global environmental modelling, EGUsphere, 2023 ([preprint](https://doi.org/10.5194/egusphere-2023-1306)).

The global sectoral water demand maps at 3 arcmin (or 0.05 degrees) resolution, 1979-2019, produced using the scripts of this utility can be downloaded from [Joint Research Centre Data Catalogue - LISFLOOD static and parameter maps for GloFAS - European Commission (europa.eu)](https://data.jrc.ec.europa.eu/dataset/68050d73-9c06-499c-a441-dc5053cb0c86)

# Data

The utility consists of five scripts which require several following datasets and files. The full list of datasets is provided below. For each dataset, this documentation provides the link and the instructions relative to available version at the time of creation of the scripts. The correct functioning of the scripts is guaranteed for those datasets. Clearly, updated versions are regularly published: one relevant example is [Global Human Settlement - Datasets - European Commission (europa.eu)](https://ghsl.jrc.ec.europa.eu/datasets.php) (see Global Human Settlement Layer (GHSL)). Albeit the functioning of the utility is expected to be maintained, the users are invited to report problems or suggestions via GitHub issues.
1. [FAO AQUASTAT](http://fao.org/aquastat/statistics/query/index.html?lang=en) sectoral water withdrawal estimates. Select "All Countries". Then Select the variable groups "Geography and population" and "Water use" (if needed select the subgroubps "Economy, development and food security" and "Water withdrawal by sector"). From these variable groups, select the following seven fields: "Gross Domestic Product (GDP)", "Industry, value added to GDP", "Agricultural water withdrawal", "Industrial water withdrawal", "Municipal water withdrawal", "Total water withdrawal", and "Irrigation water withdrawal". Select All Years, or manually select the years from 1976 to 2021. Select "Design your own" to add "Area","M49","Variable","Value","Year" fields if not present, then click on the "Show Data" button. Finally click on the "Download" button on the top right corner to download the csv file containing all data (or copy and paste the Table in CSV format shown in the web page to a csv file). Save the csv file as `aquastat_clean.csv` in `aquastat_folder` specified in the configuration file.
1. [Global Human Settlement Layer (GHSL)](https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php) POP R2023A residential population estimates for target years from 1975 to 2020. You will find data as an archive [here](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/). For each year, download the 3 arcsec (folders ending by "_3ss"), e.g. "GHS_POP_E2015_GLOBE_R2023A_4326_3ss") maps in WGS84 projection and resample them to 0.01° using `gdalwarp -t_srs EPSG:4326 -tr 0.01 0.01 -r average -of GTiff GHS_POP_E<year>_GLOBE_R2023A_4326_3ss_V1_0.tif GHS_POP_E<year>_GLOBE_R2023A_4326_3ss_V1_0_reprojected.tif`, where `<year>` is 1975 to 2020. We resample the maps using `average` instead of `sum` as the latter is broken. Put resampled files in `ghsl_folder`. You will need to specify also the 'ghsl_folder_factor = 144' to account for resampling from 0.000833 (3 arcsec) to 0.01 degrees.
1. [Global Change Analysis Model (GCAM)](https://github.com/JGCRI/gcam-core/releases) regional water withdrawal and electricity consumption estimates. Download the Windows release package version 7.1, execute `run-gcam.bat`, wait for the model to finish. Execute `run-model-interface.bat`, click "File" > "Open" > "DB Open", select `output/database_basexdb`, and select all scenarios and all regions. Select `water demand`, select `water withdrawals by sector`, and click "Run query". Select `energy transformation`, select `electricity`, select `elec consumption by demand sector`, and click "Run query". For each tab, select all data with ctrl-a and click "Edit" > "Copy". Open a blank [Google Sheets](http://sheets.google.com/) spreadsheet, press ctrl-v, manually add the headers, click "File" > "Download" > "Comma-separated values", save as `water_demand.csv` and `elec_consumption.csv`, respectively, and put both files in `gcam_folder`.
1. [Gridded Livestock of the World (GLW 3)](https://doi.org/10.1038/sdata.2018.227) species distribution dataset. Download the eight zip files (each representing a different species) and extract them to `glw_folder`.
1. [Huang et al. (2018)](https://doi.org/10.5281/zenodo.1209296) global gridded water withdrawal estimates. These estimates are not incorporated in our dataset, but are only used for the sake of comparison. Download all 7z files and extract them to `huang_folder`.
1. [Multi-Source Weather (MSWX)](http://www.gloh2o.org/mswx) daily and monthly mean air temperature. Download using rclone as explained in the FAQ on the web page. Put the daily and monthly netCDFs in `<mswx_folder>/Past/Temp/Daily` and `<mswx_folder>/Past/Temp/Monthly`, respectively, where `mswx_folder` is specified in the configuration file.
1. [US Census Bureau](https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html) state borders shape file. Download `cb_2023_us_state_500k.zip`, unzip it, open the shape file in QGIS, open the Field Calculator, enter `STATEFP` in "Expression", enter `STATEFP_` in "Output field name", and click OK. Rasterize to 0.01° using `gdal_rasterize -l cb_2023_us_state_500k -a STATEFP_ -tr 0.01 0.01 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff cb_2023_us_state_500k.shp cb_2023_us_state_500k_rasterized.tif`. Put `cb_2023_us_state_500k_rasterized.tif` in `us_states_folder`.
1. [USGS NWIS](https://waterdata.usgs.gov/nv/nwis/wu) water withdrawal estimates for 1985--present. For each state in the "Geographic Area" drop-down menu, select "State Data", "ALL Years", "State Total", and "ALL Categories" and click "Submit". Then select "Tab-separated data" and click "Submit". Do this for each state and put files in `usgs_water_use_folder`.
1. USGS water withdrawal data files for [1985](https://water.usgs.gov/watuse/data/1985/index.html) and [1990](https://water.usgs.gov/watuse/data/1990/index.html) (to supplement the NWIS data). Download "Data file of state-level data" for each year (`us85st.txt` and `us90st.txt`) and put the files in `usgs_water_use_folder`.
1. [Vassolo and Döll (2005)](https://doi.org/10.1029/2004WR003360) industrial and thermoelectric water withdrawal maps (included in this repository with permission from Petra Döll). The industrial withdrawal shape file is rasterized using `gdal_rasterize -l "industry_paper_vassolo&doell" -a MANUF_WWD -tr 0.5 0.5 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff "industry_paper_vassolo&doell.shp" manuf_wwd.tif`. Thermoelectric shape file rasterized using same command but with `WWD_PS` and `wwd_ps.tif`. Put the files in `vassolo_doll_folder`.
1. [GISCO](https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries) country borders shape file. Select year 2024, File format SHP, Geometry type Polygons (RG), Scale 01M, Coordinate system EPGS 4326, click Download. Open QGIS, then open the shape file and the csv file "un_country_codes.csv" in ancillary_data folder. Open Preprocessing->Toolbox and search for "Join Attributes by Field Value". Select the shape file as Input Layer with Table field "ISO3_CODE" and the csv file as Input layer 2 with Table field "alpha-3". Select "country-code" into "Layer 2 fields to copy" and run the preprocessing. Save the reulting layer as "CNTR_RG_01M_2024_4326".
Rasterize to 0.01° using `gdal_rasterize -l CNTR_RG_01M_2024_4326 -a country-co -tr 0.01 0.01 -a_nodata 0.0 -te -180.0 -90.0 180.0 90.0 -ot Float32 -of GTiff CNTR_RG_01M_2024_4326.shp CNTR_RG_01M_2024_4326_rasterized.tif`. Put the result in `world_borders_folder`.

1. [World Bank](https://data.worldbank.org/) manufacturing value added and gross domestic product data. Search for "Manufacturing, value added (constant 2015 US$)" and "GDP (constant 2015 US$)", download as CSV, and put in `world_bank_folder` (Remove Metadata csv files, if any).
<!---
1. USGS Circular xxx 1980 (https://pubs.usgs.gov/circ/1983/1001/report.pdf) Table 7 provided in ancillary_data folder
1. [Global Power Plant Database](https://datasets.wri.org/dataset/globalpowerplantdatabase)
1. [Lohrmann et al. (2019) power plant database](https://doi.org/10.1038/s41560-019-0501-4) Supplementary Data 1
1. [UN Industrial Commodity Statistics Database](http://data.un.org/), enter for example <beer>, then under "Select filters" select "Industry" and click "Apply filter", click "Download", click "Comma", extract and rename.
1. EIA [electricy capacity](https://www.eia.gov/international/data/world/electricity/), Select World, select Electricity, select Capacity and Generation, click View Data, click Sort by Energy Source/Activity, click Download Options, select Export CSV (table).
-->

The locations of the datasets and files are specified in the configuration file (see the `config.cfg` example). The script also requires a template map which defines the output resolution and area. The template map should be in netCDF-4 format and contain `lat` and `lon` variables and a data variable (any name). The location of the template map is specified using `templatemap_path` in the configuration file. The data are produced for the period spanning `year_start` to `year_end` and saved in netCDF-4 format to `output_folder` (all specified in the configuration file). 

# Methods

The module consists of five scripts:
1. `step1_population_density.py`: 
	1. Resamples the native resolution of the GHSL population grids to the template map resolution.
	1. Computes population grids (in density per km<sup>2</sup>) for each year using linear interpolation and nearest-neighbor extrapolation.
1. `step2_domestic_demand.py`:
	1. Produces a map of the R parameter (measures the relative withdrawal difference between the warmest and coldest months).
	1. Computes annual time series of population for each country.
	1. Loads the country-scale AQUASTAT municipal water withdrawal estimate and produces annual time series using linear interpolation and population-based forward and backward extrapolation.
	1. Loads the state-level USGS withdrawal data and computes annual values using linear interpolation and nearest-neighbour extrapolation.
	1. For each year, computes country- and state-scale per capita water demand and produces a global map.
	1. Annual per capita water demand for countries without any AQUASTAT information is estimated using nearest-neighbour interpolation.
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

The script can be run on a normal desktop PC, support is provided only for Linux. Physical memory requirements depend on the desired geographical extent and resolution of the maps dataset. For instance, the generation of  the sectoral water demand dataset for the global domain at 3arcmin (or 0.05 degrees) resolution required 32 GB.

# Instructions
Install main lisflood-utilities package in a conda environment:
```bash
conda create --name <env> python=3.7 -c conda-forge
conda activate <env>
conda install -c conda-forge pcraster eccodes gdal
pip install lisflood-utilities
```
Then clone the repository:
```
git clone --single-branch --branch master https://github.com/ec-jrc/lisflood-utilities
cd lisflood-utilities/src/lisfloodutilities/water-demand-historic
```
Produce a configuration file with the correct paths and folders (based on the [included example](config_global_0p1deg.cfg)). 

Run the scripts as follows:
```
python step1_population_density.py <config file>
python step2_domestic_demand.py <config file>
...
```