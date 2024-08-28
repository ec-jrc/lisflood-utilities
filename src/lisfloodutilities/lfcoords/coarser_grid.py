import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray
import pyflwdir
from pathlib import Path
from tqdm import tqdm
from typing import Optional
import logging
import warnings
warnings.filterwarnings("ignore")

from lisfloodutilities.lfcoords import Config
from lisfloodutilities.lfcoords.utils import catchment_polygon

# set logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def coordinates_coarse(
    cfg: Config,
    stations: pd.DataFrame,
    save: bool = True
) -> Optional[pd.DataFrame]:
    """
    Transforms station coordinates from a high-resolution grid to a corresponding location in a coarser grid, aiming to match the shape of the catchment area derived from the high-resolution map. It updates the station coordinates and exports the catchment areas as shapefiles in the coarser grid.

    The function reads the upstream area map and local drainage direction (LDD) map in the coarse grid. It then finds the pixel in the coarse grid that best matches the catchment shape derived from the high-resolution map. The match is evaluated based on the intersection-over-union of catchment shapes and the ratio of upstream areas.

    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in the configuration file.
    stations: pandas.DataFrame
        DataFrame containing station coordinates and upstream areas in the finer grid. It's the result of coarse_grid.coarse_grid()
    save: bool
        If True, the updated stations table is saved to a CSV file.
        If False, the updated stations DataFrame is returned without saving.

    Returns:
    --------
    staions: pandas.DataFrame
        If save is False, returns a pandas DataFrame with updated station coordinates and upstream areas in the coarser grid. Otherwise, the function returns None, and the results are saved directly to a CSV file.
    """

    ### READ INPUTS

    # read upstream area map of coarse grid
    upstream_coarse = rioxarray.open_rasterio(cfg.UPSTREAM_COARSE).squeeze(dim='band')
    logger.info(f'Map of upstream area corretly read: {cfg.UPSTREAM_COARSE}')

    # read local drainage direction map
    ldd_coarse = rioxarray.open_rasterio(cfg.LDD_COARSE).squeeze(dim='band')
    logger.info(f'Map of local drainage directions correctly read: {cfg.LDD_COARSE}')


    ### PROCESSING
    
    # create river network
    fdir_coarse = pyflwdir.from_array(ldd_coarse.data,
                                      ftype='ldd',
                                      transform=ldd_coarse.rio.transform(),
                                      check_ftype=False,
                                      latlon=True)

    # boundaries of the input maps
    lon_min, lat_min, lon_max, lat_max = np.round(ldd_coarse.rio.bounds(), 6)

    # resolution of the input maps
    cellsize = np.round(np.mean(np.diff(ldd_coarse.x)), 6) # degrees
    cellsize_arcmin = int(np.round(cellsize * 60, 0)) # arcmin
    suffix_coarse = f'{cellsize_arcmin}min'
    logger.info(f'Coarse resolution is {cellsize_arcmin} arcminutes')

    # extract resolution of the finer grid from 'stations'
    suffix_fine = sorted(set([col.split('_')[1] for col in stations.columns if '_' in col]))[0]
    cols_fine = [f'{col}_{suffix_fine}' for col in ['area', 'lat', 'lon']]

    # add new columns to 'stations'
    cols_coarse = [f'{col}_{suffix_coarse}' for col in ['area', 'lat', 'lon']]
    stations[cols_coarse] = np.nan

    # output folders
    SHAPE_FOLDER_FINE = cfg.SHAPE_FOLDER / suffix_fine
    SHAPE_FOLDER_COARSE = cfg.SHAPE_FOLDER / suffix_coarse
    SHAPE_FOLDER_COARSE.mkdir(parents=True, exist_ok=True)

    # search range of 5x5 array -> this is where the best point can be found in the coarse grid
    rangexy = np.linspace(-2, 2, 5) * cellsize # arcmin
    for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):

        # real upstream area
        area_ref = attrs['area']

        # coordinates and upstream area in the fine grid
        lat_fine, lon_fine, area_fine = attrs[[f'{col}_{suffix_fine}' for col in ['lat', 'lon', 'area']]]

        if (area_ref < cfg.MIN_AREA) or (area_fine < cfg.MIN_AREA):
            logger.warning(f'The catchment area of station {ID} is smaller than the minimum of {cfg.MIN_AREA} km2')
            continue

        # import shapefile of catchment polygon
        shapefile = SHAPE_FOLDER_FINE / f'{ID}.shp'
        try:
            basin_fine = gpd.read_file(shapefile)
            logger.info(f'Catchment polygon correctly read: {shapefile}')
        except OSError as e:
            logger.error(f'Error reading {shapefile}: {e}')
            continue
        except Exception as e:  # This will catch other exceptions that might occur.
            logger.error(f'An unexpected error occurred while reading {shapefile}: {e}')
            continue

        # find ratio
        logger.debug('Start search')
        inter_vs_union, area_ratio, area_lisf = [], [], []
        for Δlat in rangexy:
            for Δlon in rangexy:
                lon = lon_fine + Δlon
                lat = lat_fine + Δlat
                basin = catchment_polygon(fdir_coarse.basins(xy=(lon, lat)).astype(np.int32),
                                          transform=ldd_coarse.rio.transform(), 
                                          crs=ldd_coarse.rio.crs)

                # calculate union and intersection of shapes
                intersection = gpd.overlay(basin_fine, basin, how='intersection')
                union = gpd.overlay(basin_fine, basin, how='union')
                inter_vs_union.append(intersection.area.sum() / union.area.sum())

                # get upstream area (km2) of coarse grid (LISFLOOD)
                area = upstream_coarse.sel(x=lon, y=lat, method='nearest').item() * 1e-6
                area_lisf.append(area)

                # ratio between reference and coarse area
                if area_ref == 0 or area == 0:
                    ratio = 0
                else:
                    ratio = area_ref / area if area_ref < area else area / area_ref
                area_ratio.append(ratio)
        logger.debug('End search')

        # maximum of shape similarity and upstream area accordance
        i_shape = np.argmax(inter_vs_union)
        area_shape = area_lisf[i_shape]
        i_centre = int(len(rangexy)**2 / 2) # middle point
        area_centre = area_lisf[i_centre]           
        # use middle point if errors are small
        abs_error = abs(area_shape - area_centre)
        pct_error = 100 * abs(1 - area_centre / area_shape)
        if (abs_error <= cfg.ABS_ERROR) and (pct_error <= cfg.PCT_ERROR):
            i_shape = i_centre
            area_shape = area_centre

        #i_ratio = np.argmax(area_ratio)          

        # coordinates in the fine resolution
        i = i_shape // len(rangexy)
        j = i_shape % len(rangexy)
        lat = lat_fine + rangexy[i]
        lon = lon_fine + rangexy[j]

        # coordinates and upstream area on coarse resolution
        area = upstream_coarse.sel(x=lon, y=lat, method='nearest')
        area_coarse = area.item() * 1e-6
        lon_coarse = area.x.item()
        lat_coarse = area.y.item()

        # derive catchment polygon from the selected coordinates
        basin_coarse = catchment_polygon(fdir_coarse.basins(xy=(lon_coarse, lat_coarse)).astype(np.int32),
                                         transform=ldd_coarse.rio.transform(), 
                                         crs=ldd_coarse.rio.crs,
                                         name='ID')
        basin_coarse['ID'] = ID
        basin_coarse.set_index('ID', inplace=True)
        basin_coarse[cols_fine] = area_fine, lat_fine, lon_fine
        basin_coarse[cols_coarse] = area_coarse, lat_coarse, lon_coarse

        # export shapefile
        output_shp = SHAPE_FOLDER_COARSE / f'{ID}.shp'
        basin_coarse.to_file(output_shp)
        logger.info(f'Catchment {ID} exported as shapefile: {output_shp}')

        # update new columns in 'stations'
        stations.loc[ID, cols_coarse] = [int(area_coarse), round(lat_coarse, 6), round(lon_coarse, 6)]
    
    # return/save
    stations.sort_index(axis=1, inplace=True)
    if save:
        output_csv = cfg.STATIONS.parent / f'{cfg.STATIONS.stem}_{suffix_coarse}.csv'
        stations.to_csv(output_csv)
        logger.info(f'The updated stations table in the coarser grid has been exported to: {output_csv}')
    else:
        return stations
