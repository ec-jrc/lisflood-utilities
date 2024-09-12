import os
os.environ['USE_PYGEOS'] = '0'
import yaml
from pathlib import Path
from typing import Union, Dict
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import rioxarray
import logging

# set logger
logger = logging.getLogger('lfcoords')


class Config:
    
    def __init__(self, config_file):
        """
        Reads the configuration from a YAML file and sets default values if not provided.

        Parameters:
        -----------
        config_file: string or pathlib.Path
            The path to the YAML configuration file.
        """
        
        # read configuration file
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            config = yaml.load(ymlfile, Loader=yaml.FullLoader)
            
        # input
        self.POINTS = Path(config['input']['points'])
        self.LDD_FINE = Path(config['input']['ldd_fine'])
        self.UPSTREAM_FINE = Path(config['input']['upstream_fine'])
        self.LDD_COARSE = Path(config['input']['ldd_coarse'])
        self.UPSTREAM_COARSE = Path(config['input']['upstream_coarse'])
        
        # resolutions
        self.FINE_RESOLUTION = None
        self.COARSE_RESOLUTION = None
        
        # output
        self.OUTPUT_FOLDER = Path(config.get('output_folder', './shapefiles'))
        self.OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
        
        # conditions
        self.MIN_AREA = config['conditions'].get('min_area', 10)
        self.ABS_ERROR = config['conditions'].get('abs_error', 50)
        self.PCT_ERROR = config['conditions'].get('pct_error', 1)
        
    def update_config(
        self,
        fine_grid: xr.DataArray,
        coarse_grid: xr.DataArray
    ):
        """It extracts the resolution of the finer and coarser grid, updates the respective attributes in the configuration object, and it creates the necessary structure of directories

        Parameters:
        -----------
        fine_grid: xarray.DataArray
            Any map in the fine grid
        coarse_grid: xarray.DataArray
            Any map in the coarse grid    
        """

        # resolution of the finer grid
        cellsize = np.mean(np.diff(fine_grid.x)) # degrees
        cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
        logger.info(f'The resolution of the finer grid is {cellsize_arcsec} arcseconds')
        self.FINE_RESOLUTION = f'{cellsize_arcsec}sec'

        # resolution of the input maps
        cellsize = np.round(np.mean(np.diff(coarse_grid.x)), 6) # degrees
        cellsize_arcmin = int(np.round(cellsize * 60, 0)) # arcmin
        logger.info(f'The resolution of the coarser grid is {cellsize_arcmin} arcminutes')
        self.COARSE_RESOLUTION = f'{cellsize_arcmin}min'
        
        
def read_input_files(
    cfg: Config
) -> Dict:
    """Reads the input files and saves them in a dictionary. It also updates attributes in the Config object based on the resolution of the input maps
    
    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in the configuration file.
        
    Returns:
    --------
    inputs: dictionary
        * 'points': geopandas.GeoDataFrame of input points
        * 'ldd_fine': xarray.DataArray of local drainaige directions in the fine grid
        * 'upstream_fine': xarray.DataArray of upstream area (km2) in the fine grid
        * 'ldd_coarse': xarray.DataArray of local drainaige directions in the coarse grid
        * 'upstream_coarse': xarray.DataArray of upstream area (m2) in the coarse grid
    """
    
    # read upstream map with fine resolution
    upstream_fine = rioxarray.open_rasterio(cfg.UPSTREAM_FINE).squeeze(dim='band')
    logger.info(f'Map of upstream area in the finer grid corretly read: {cfg.UPSTREAM_FINE}')

    # read local drainage direction map
    ldd_fine = rioxarray.open_rasterio(cfg.LDD_FINE).squeeze(dim='band')
    logger.info(f'Map of local drainage directions in the finer grid corretly read: {cfg.LDD_FINE}')
    
    # read upstream area map of coarse grid
    upstream_coarse = rioxarray.open_rasterio(cfg.UPSTREAM_COARSE).squeeze(dim='band')
    logger.info(f'Map of upstream area in the coarser grid corretly read: {cfg.UPSTREAM_COARSE}')

    # read local drainage direction map
    ldd_coarse = rioxarray.open_rasterio(cfg.LDD_COARSE).squeeze(dim='band')
    logger.info(f'Map of local drainage directions in the coarser grid correctly read: {cfg.LDD_COARSE}')
    
    # read points text file
    points = pd.read_csv(cfg.POINTS, index_col='ID')
    points.columns = points.columns.str.lower()
    logger.info(f'Table of points correctly read: {cfg.POINTS}')
    points = check_points(cfg, points, ldd_fine)
    
    # convert to geopandas and export as shapefile
    points = gpd.GeoDataFrame(points,
                              geometry=[Point(xy) for xy in zip(points['lon'], points['lat'])],
                              crs=ldd_coarse.rio.crs)
    point_shp = cfg.OUTPUT_FOLDER / f'{cfg.POINTS.stem}.shp'
    points.to_file(point_shp)
    logger.info(f'The original points table has been exported to: {point_shp}')
    
    inputs = {
        'points': points,
        'ldd_fine': ldd_fine,
        'upstream_fine': upstream_fine,
        'ldd_coarse': ldd_coarse,
        'upstream_coarse': upstream_coarse,
    }
    
    # update Config
    cfg.update_config(ldd_fine, ldd_coarse)
    
    return inputs


def check_points(
    cfg: Config,
    points: pd.DataFrame,
    ldd: xr.DataArray    
    ) -> pd.DataFrame:
    """Removes input points with missing value, catchment area smaller than the predefined threshold, or outside the extent of the input map
    
    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in the configuration file.
    points: pandas.DataFrame
        Table of input points with fields 'lat', 'lon' and 'area' (km2)
    ldd: xarray.DataArray
        Map of local drainage directions
        
    Returns:
    --------
    points: pandas.DataFrame
        The input table in which points with conflicts have been removed
    """
    
    # remove points with missing values
    mask_nan = points.isnull().any(axis=1)
    if mask_nan.sum() > 0:
        points = points[~mask_nan]
        logger.warning(f'{mask_nan.sum()} points were removed because of missing values')
        
    # remove points with small catchment area
    mask_area = points['area'] < cfg.MIN_AREA
    if mask_area.sum() > 0:
        points = points[~mask_area]
        logger.info(f'{mask_area.sum()} points were removed due to their small catchment area')
        
    # remove points outside the input LDD map
    lon_min, lat_min, lon_max, lat_max = np.round(ldd.rio.bounds(), 6)
    mask_lon = (points.lon < lon_min) | (points.lon > lon_max)
    mask_lat = (points.lat < lat_min) | (points.lat > lat_max)
    mask_extent = mask_lon | mask_lat
    if mask_extent.sum() > 0:
        points = points[~mask_extent]
        logger.info(f'{mask_extent.sum()} points were removed because they are outside the input LDD map')
        
    return points
