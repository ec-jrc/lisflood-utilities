import os
import logging
import yaml
from pathlib import Path
from typing import Dict
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray as rxr

# set logger
logger = logging.getLogger(__name__)

# Use a module-level constant for the environment variable setting
os.environ['USE_PYGEOS'] = '0'

class Config:
    """
    Manages the application's configuration by reading a YAML file
    and setting default values.
    """
    
    def __init__(self, config_file: Path):
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
            
        # input file paths
        self.points = Path(config['input']['points'])
        self.ldd_fine = Path(config['input']['ldd_fine'])
        self.upstream_fine = Path(config['input']['upstream_fine'])
        self.ldd_coarse = Path(config['input']['ldd_coarse'])
        self.upstream_coarse = Path(config['input']['upstream_coarse'])
        
        # resolutions
        self.fine_resolution = None
        self.coarse_resolution = None
        
        # output folder
        self.output_folder = Path(config.get('output_folder', './shapefiles'))
        self.output_folder.mkdir(parents=True, exist_ok=True)
        
        # conditions
        self.min_area = config['conditions'].get('min_area', 10)
        self.abs_error = config['conditions'].get('abs_error', 50)
        self.pct_error = config['conditions'].get('pct_error', 1)
        
    def update_config(
        self,
        fine_grid: xr.DataArray,
        coarse_grid: xr.DataArray
    ):
        """
        Extracts the resolution from the finer and coarser grids and updates the configuration object.

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
        self.fine_resolution = f'{cellsize_arcsec}sec'

        # resolution of the input maps
        cellsize = np.round(np.mean(np.diff(coarse_grid.x)), 6) # degrees
        cellsize_arcmin = int(np.round(cellsize * 60, 0)) # arcmin
        logger.info(f'The resolution of the coarser grid is {cellsize_arcmin} arcminutes')
        self.coarse_resolution = f'{cellsize_arcmin}min'
        
 
def read_input_files(cfg: Config) -> Dict:
    """
    Reads input files, updates the Config object, and returns a dictionary
    of the loaded data.
    
    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters.
        
    Returns:
    --------
    Dict
        A dictionary containing the loaded data:
        * 'points': geopandas.GeoDataFrame of input points
        * 'ldd_fine': xarray.DataArray of local drainage directions in the fine grid
        * 'upstream_fine': xarray.DataArray of upstream area (km2) in the fine grid
        * 'ldd_coarse': xarray.DataArray of local drainage directions in the coarse grid
        * 'upstream_coarse': xarray.DataArray of upstream area (m2) in the coarse grid
    """

    # a helper function to reduce code repetition
    def open_raster(path):
        """Helper to open and squeeze a raster file."""
        return rxr.open_rasterio(path).squeeze(dim='band')
        
    # read upstream map with fine resolution
    upstream_fine = rxr.open_rasterio(cfg.upstream_fine).squeeze(dim='band')
    logger.info(f'Map of upstream area in the finer grid corretly read: {cfg.upstream_fine}')

    # read local drainage direction map
    ldd_fine = rxr.open_rasterio(cfg.ldd_fine).squeeze(dim='band')
    logger.info(f'Map of local drainage directions in the finer grid corretly read: {cfg.ldd_fine}')
    
    # read upstream area map of coarse grid
    upstream_coarse = rxr.open_rasterio(cfg.upstream_coarse).squeeze(dim='band')
    logger.info(f'Map of upstream area in the coarser grid corretly read: {cfg.upstream_coarse}')

    # read local drainage direction map
    ldd_coarse = rxr.open_rasterio(cfg.ldd_coarse).squeeze(dim='band')
    logger.info(f'Map of local drainage directions in the coarser grid correctly read: {cfg.ldd_coarse}')
    
    # read points text file
    points = pd.read_csv(cfg.points, index_col='ID')
    points.columns = points.columns.str.lower()
    logger.info(f'Table of points correctly read: {cfg.points}')
    points = check_points(cfg, points, ldd_fine)
    
    # convert to geopandas and export as shapefile
    points = gpd.GeoDataFrame(
        points,
        geometry=gpd.points_from_xy(points['lon'], points['lat']),
        crs=ldd_coarse.rio.crs
    )
    point_shp = cfg.output_folder / f'{cfg.points.stem}.shp'
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
    """
    Removes input points that have missing values, a small catchment area,
    or are outside the map extent.
    
    Parameters:
    -----------
    cfg: Config
        Configuration object.
    points: pandas.DataFrame
        Table of input points with fields 'lat', 'lon' and 'area' (km2)
    ldd: xarray.DataArray
        Map of local drainage directions
        
    Returns:
    --------
    pandas.DataFrame
        The input table with points with conflicts removed.
    """
    
    # remove points with missing values
    mask_nan = points.isnull().any(axis=1)
    if mask_nan.sum() > 0:
        points = points[~mask_nan]
        logger.warning(f'{mask_nan.sum()} points were removed because of missing values.')
        
    # remove points with small catchment area
    mask_area = points['area'] < cfg.min_area
    if mask_area.sum() > 0:
        points = points[~mask_area]
        logger.info(f'{mask_area.sum()} points were removed due to their small catchment area.')
        
    # remove points outside the input LDD map
    lon_min, lat_min, lon_max, lat_max = np.round(ldd.rio.bounds(), 6)
    mask_lon = (points.lon < lon_min) | (points.lon > lon_max)
    mask_lat = (points.lat < lat_min) | (points.lat > lat_max)
    mask_extent = mask_lon | mask_lat
    if mask_extent.sum() > 0:
        points = points[~mask_extent]
        logger.info(f'{mask_extent.sum()} points were removed because they are outside the input LDD map.')
        
    return points