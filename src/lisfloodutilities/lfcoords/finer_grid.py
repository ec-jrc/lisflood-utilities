import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import rioxarray
import pyflwdir
from tqdm import tqdm
from pathlib import Path
from typing import Optional, Union
import logging
import warnings
warnings.filterwarnings("ignore")

from lisfloodpreprocessing import Config
from lisfloodpreprocessing.utils import find_pixel, catchment_polygon 

# set logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def coordinates_fine(
    cfg: Config,
    points: pd.DataFrame,
    ldd_fine: xr.DataArray,
    upstream_fine: xr.DataArray,
    save: bool = False
) -> Optional[gpd.GeoDataFrame]:
    """
    Processes point coordinates to find the most accurate pixel in a high-resolution map, 
    based on a reference value of catchment area. It updates the station coordinates and 
    exports the catchment areas as shapefiles.

    The function reads the upstream area map and local drainage direction (LDD) map in 
    fine resolution, as well as the station coordinates with their reference upstream areas.
    It then iteratively searches for the best-matching pixel within an increasing search 
    range and applies penalties and factors to determine the closest match. The function 
    creates a boolean map of each station's catchment area and vectorizes it into a 
    GeoDataFrame for export as a shapefile.

    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in  the configuration file.
    points: pandas.DataFrame
        DataFrame containing reference point coordinates and upstream areas
    ldd_fine: xarray.DataArray
        Map of local drainaige directions in the fine grid
    upstream_fine: xarray.DataArray
        Map of upstream area (km2) in the fine grid
    save: boolean
        If True, the updated table of points is exported as a shapefile.
        
    Returns:
    --------
    points: geopandas.GeoDataFrame
        A table with updated station coordinates and upstream areas in the finer grid.
    polygons: geopandas.GeoDataFrame
        A table with the catchment polygons in the finer grid.
    """
    
    # resolution of the input map
    cellsize = np.mean(np.diff(upstream_fine.x)) # degrees
    cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
    # cfg.FINE_RESOLUTION = f'{cellsize_arcsec}sec'
    fine_resolution = f'{cellsize_arcsec}sec'
    
    # add columns to the table of points
    new_cols = sorted([f'{col}_{fine_resolution}' for col in ['lat', 'lon', 'area']])
    points[new_cols] = np.nan

    # create river network
    fdir_fine = pyflwdir.from_array(ldd_fine.data,
                                    ftype='d8', 
                                    transform=ldd_fine.rio.transform(),
                                    check_ftype=False,
                                    latlon=True)
    
    polygons = []
    for ID, attrs in tqdm(points.iterrows(), total=points.shape[0], desc='points'):  

        # reference coordinates and upstream area
        lat_ref, lon_ref, area_ref = attrs[['lat', 'lon', 'area']]

        # search new coordinates in an increasing range
        ranges = [55, 101, 151]
        penalties = [500, 500, 1000]
        factors = [2, .5, .25]
        acceptable_errors = [50, 80, np.nan]
        for rangexy, penalty, factor, max_error in zip(ranges, penalties, factors, acceptable_errors):
            logger.debug(f'Set range to {rangexy}')
            lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=penalty, factor=factor)
            if error <= max_error:
                break

        # update new columns in 'points'
        points.loc[ID, new_cols] = [int(upstream_fine.sel(y=lat, x=lon).item()), round(lat, 6), round(lon, 6)]

        # boolean map of the catchment associated to the corrected coordinates
        basin_arr = fdir_fine.basins(xy=(lon, lat)).astype(np.int32)

        # vectorize the boolean map into geopandas
        basin_gdf = catchment_polygon(basin_arr.astype(np.int32),
                                      transform=ldd_fine.rio.transform(),
                                      crs=ldd_fine.rio.crs,
                                      name='ID')
        basin_gdf['ID'] = ID
        basin_gdf.set_index('ID', inplace=True)
        basin_gdf[attrs.index] = attrs.values
        
        # save polygon
        polygons.append(basin_gdf)
        
    # concatenate polygons shapefile
    polygons = pd.concat(polygons)
    
    # convert points to geopandas
    geometry = [Point(xy) for xy in zip(points[f'lon_{fine_resolution}'], points[f'lat_{fine_resolution}'])]
    points = gpd.GeoDataFrame(points, geometry=geometry, crs=4326)
    points.sort_index(axis=1, inplace=True)
    
    if save is True:
        # polygons
        polygon_shp = cfg.OUTPUT_FOLDER / f'catchments_{fine_resolution}.shp'
        polygons.to_file(polygon_shp)
        logger.info(f'Catchments in the finer grid have been exported to: {polygon_shp}')
        # points
        point_shp = cfg.OUTPUT_FOLDER / f'{cfg.POINTS.stem}_{fine_resolution}.shp'
        points.to_file(point_shp)
        logger.info(f'The updated points table in the finer grid has been exported to: {point_shp}')
        
    return points, polygons