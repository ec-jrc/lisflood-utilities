import os
import logging
from pathlib import Path
from typing import Optional, Union, Tuple
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import pyflwdir
from tqdm import tqdm

from lisfloodpreprocessing import Config
from lisfloodpreprocessing.utils import find_pixel, catchment_polygon

os.environ['USE_PYGEOS'] = '0'
warnings.filterwarnings("ignore")

# set logger
logger = logging.getLogger(__name__)


def coordinates_fine(
    cfg: Config,
    points: pd.DataFrame,
    ldd_fine: xr.DataArray,
    upstream_fine: xr.DataArray,
    save: bool = False
) -> Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]:
    """
    Processes point coordinates to find the most accurate pixel in a high-resolution
    map, based on a reference value of catchment area. It updates the station
    coordinates and exports the catchment areas as shapefiles.

    Parameters
    ----------
    cfg : Config
        Configuration object containing file paths and parameters.
    points : pd.DataFrame
        DataFrame containing reference point coordinates and upstream areas.
    ldd_fine : xr.DataArray
        Map of local drainage directions in the fine grid.
    upstream_fine : xr.DataArray
        Map of upstream area (km2) in the fine grid.
    save : bool, optional
        If True, the updated table of points and catchments are exported.

    Returns
    -------
    Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]
        A tuple containing:
        - A table with updated station coordinates and upstream areas.
        - A table with the catchment polygons in the finer grid.
        Returns None if no polygons are generated.
    """
    
    # add columns to the table of points
    points_fine = points.copy()
    n_points = points.shape[0]
    cols = ['lat', 'lon', 'area']
    new_cols = sorted([f'{col}_{cfg.fine_resolution}' for col in cols])
    points_fine[new_cols] = np.nan

    # create river network
    fdir_fine = pyflwdir.from_array(
        ldd_fine.data,
        ftype='d8',
        transform=ldd_fine.rio.transform(),
        check_ftype=False,
        latlon=True
    )
    
    polygons_fine = []
    for point_id, attrs in tqdm(points.iterrows(), total=n_points, desc='points'):    
        try:
            # reference coordinates and upstream area
            lat_ref, lon_ref, area_ref = attrs[cols]

            # search new coordinates in an increasing range
            ranges = [55, 101, 151]
            penalties = [500, 500, 1000]
            factors = [2, 0.5, 0.25]
            acceptable_errors = [50, 80, np.nan]
            for range_xy, penalty, factor, max_error in zip(ranges, penalties, factors, acceptable_errors):
                logger.debug(f'Set range to {range_xy}')
                lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, range_xy=range_xy, penalty=penalty, factor=factor)
                if error <= max_error:
                    break

            # update new columns in 'points_fine'
            points_fine.loc[point_id, new_cols] = [int(upstream_fine.sel(y=lat, x=lon).item()), round(lat, 6), round(lon, 6)]

            # boolean map of the catchment
            basin_arr = fdir_fine.basins(xy=(lon, lat)).astype(np.int32)

            # vectorize the boolean map into geopandas
            basin_gdf = catchment_polygon(
                basin_arr.astype(np.int32),
                transform=ldd_fine.rio.transform(),
                crs=ldd_fine.rio.crs,
                name='ID'
            )
            basin_gdf['ID'] = point_id
            basin_gdf[cols] = attrs[cols]
            basin_gdf.set_index('ID', inplace=True)

            # save polygon
            polygons_fine.append(basin_gdf)

            # logger.info(f'Point {point_id} located in the finer grid')
        except Exception as e:
            logger.error(f'Point {point_id} could not be located in the finer grid: {e}')

    # handle case where no polygons were generated
    if not polygons_fine:
        logger.warning('No points could be located in the finer grid. Returning empty dataframes.')
        return gpd.GeoDataFrame(), gpd.GeoDataFrame()
        
    # concatenate polygons shapefile
    polygons_fine = pd.concat(polygons_fine)
    
    # convert points to geopandas
    points_fine = gpd.GeoDataFrame(
        points_fine, 
        geometry=gpd.points_from_xy(points_fine[f'lon_{cfg.fine_resolution}'], points_fine[f'lat_{cfg.fine_resolution}']), 
        crs=ldd_fine.rio.crs
    )
    points_fine.sort_index(axis=1, inplace=True)
    
    # compute error
    points_fine['abs_error'] = abs(points_fine[f'area_{cfg.fine_resolution}'] - points_fine['area'])
    points_fine['pct_error'] = points_fine.abs_error / points_fine['area'] * 100
    
    if save:
        # polygons
        polygon_shp = cfg.output_folder / f'catchments_{cfg.fine_resolution}.shp'
        polygons_fine.to_file(polygon_shp)
        logger.info(f'Catchments in the finer grid have been exported to: {polygon_shp}')
        
        # points
        point_shp = cfg.output_folder / f'{cfg.points.stem}_{cfg.fine_resolution}.shp'
        points_fine.to_file(point_shp)
        logger.info(f'The updated points table in the finer grid has been exported to: {point_shp}')
        
    return points_fine, polygons_fine