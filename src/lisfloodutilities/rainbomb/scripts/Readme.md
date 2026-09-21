### Generation of forcings for GloFAS v5

This folder provides all the data needed for generating the GloFAS v5 forcings.

The forcings are based on ERA5. The steps needed are:
1. Retrieve hourly data from MARS for the variables of interest, using the native spatial resolution of ERA5.
2. Get daily values from the hourly data.
3. Process rainfall daily fields for correcting the ERA5 rainbombs. (auxiliary data available in *rainbomb_correction* subdirectory)
4. Interpolate the variables to the grid of GloFAS with pyg2p ([**v3.2.7**](https://github.com/ec-jrc/pyg2p/tree/v3.2.7)). (auxiliary data available in *glofas_execution_templates* subdirectory)
5. Use LISVAP ([**v1.3.1**](https://github.com/ec-jrc/lisflood-lisvap/releases/tag/v1.3.1)) for generating the evapotranspiration variables that are needed for LISFLOOD.
6. Mask all grid cells that are not used in GloFAS (e.g., sea).
7. Save final data as 1 nc file per year.

## Useful links:
1. [Confluence page for input data from MARS](https://smhi-efas.atlassian.net/wiki/spaces/ECC/pages/32789274/Forcing+generation)
2. [Confluence page for ERA5 rainbomb correction](https://smhi-efas.atlassian.net/wiki/spaces/ECC/pages/32785194/ERA5+rainbombs+correction)


## Workflow
The user should run the script *GloFAS_v5_ERA5_forcings_generation.ksh*. This will call the template script *GloFAS_v5_ERA5_forcings_gen_template.ksh*, that will produce for each year a different task for forcings generation.
The user should have access to MARS, so the raw forcing data on the native resolution can be downloaded.
Also there is a need to have various packages installed, as LISVAP, pyg2p, pcraster, etc.

## Notes
- The pyg2p templates have the original paths for the elevation and templates files. Copies of these 2 files are also available in this repository (elv_Global_03min.nc, template_Global_03min.nc)
- The scipts keep also the original rainfall forcings without the rainbomb correction, in case changes are needed, or an assessment on the hydrological impact due to the updated forcings is requested. The files are kept both in the native ERA5 spatial resolution, and the interpolated GloFAS grid.