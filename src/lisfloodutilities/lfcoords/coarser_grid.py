import logging
from typing import Union, Tuple
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from . import Config
from .io import read_flowdir, save_results
from .utils import check_points, define_neighbourhood, compute_area_error, compute_shape_error


# set logger
logger = logging.getLogger(__name__)


def coordinates_coarse(
    cfg: Config,
    points_fine: Union[pd.DataFrame, gpd.GeoDataFrame],
    polygons_fine: gpd.GeoDataFrame,
    n_cell: int = 3,
    save: bool = False
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Transforms point coordinates from a high-resolution grid to a corresponding
    location in a coarser grid, aiming to match the shape of the catchment
    area. It updates station coordinates and exports catchments as shapefiles.

    Parameters
    ----------
    cfg : Config
        Configuration object with file paths and parameters.
    points_fine : pd.DataFrame or gpd.GeoDataFrame
        Table with updated station coordinates and upstream areas from a finer grid.
    polygons_fine : gpd.GeoDataFrame
        Table with the catchment polygons from a finer grid.
    flwdir_coarse : pyflwdir.FlwdirRaster
        Map of local drainage directions in the coarse grid.
    n_cell: integer
        Distance (in number of cells) to explore around the input point.
    save : bool, optional
        If True, the updated tables are exported as shapefiles.

    Returns
    -------
    Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        A tuple containing:
        - A table with updated station coordinates in the coarser grid.
        - A table with the catchment polygons in the coarser grid.
    """

    # read flow direction map
    try: 
        flwdir, resolution = read_flowdir(cfg.flwdir_coarse, ftype='ldd') # type: ignore
        cfg.coarse_resolution = resolution
        logger.info(f'Flow direction map at {resolution} resolution correctly read: {cfg.flwdir_coarse}')
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        raise RuntimeError("Aborting computation due to missing or corrupt flow direction map.") from e

    # check input points
    # Preserve CRS before check_points (which returns pd.DataFrame, losing CRS info)
    points_crs = points_fine.crs if isinstance(points_fine, gpd.GeoDataFrame) else None
    points_fine = check_points(cfg, points_fine, flwdir.extent)
    n_points = points_fine.shape[0]

    # upstream area in km²
    uparea = flwdir.upstream_area(unit='km2')

    # create output point layer
    points_coarse = points_fine.copy()
    cols = ['area', 'lat', 'lon']
    cols_fine = [f'{col}_{cfg.fine_resolution}' for col in cols]
    cols_coarse = [f'{col}_{cfg.coarse_resolution}' for col in cols]
    points_coarse[cols_coarse] = np.nan

    polygons_coarse = []
    for ID, point in tqdm(points_fine.iterrows(), total=n_points, desc='points'):

        try:
            # real upstream area and coordinates
            area_ref, lat_ref, lon_ref = point[cols]

            # coordinates and upstream area in the fine grid
            area_fine, lat_fine, lon_fine = point[cols_fine]

            # define pixels to explore
            idxs, xs, ys = define_neighbourhood(flwdir, lon_fine, lat_fine, n_cell)
        
            # compute error in catchment area
            area, area_error = compute_area_error(uparea, idxs, area_fine)

            # compute error in shape similarity
            # Use .iloc with the position and explicitly copy to avoid SettingWithCopyWarning
            idx_position = points_fine.index.get_loc(ID)
            polygon_ref = polygons_fine.iloc[[idx_position]].copy() # type: ignore
            basins, shape_error = compute_shape_error(flwdir, idxs, polygon_ref)

            # select best pixel based on the product of area and shape error
            i = np.argmin(area_error + shape_error)
            lon_coarse, lat_coarse, area_coarse = round(xs[i].item(), 6), round(ys[i].item(), 6), round(area[i].item(), 3)

            # update new columns in 'points_coarse'
            points_coarse.loc[ID, cols_coarse] = [area_coarse, lat_coarse, lon_coarse]

            # get the associated basin polygon
            basin_coarse = basins.iloc[[i]] # type: ignore
            basin_coarse['ID'] = ID
            basin_coarse.set_index('ID', inplace=True)
            basin_coarse[cols] = area_ref, lat_ref, lon_ref
            basin_coarse[cols_fine] = area_fine, lat_fine, lon_fine
            basin_coarse[cols_coarse] = area_coarse, lat_coarse, lon_coarse

            # save polygon
            polygons_coarse.append(basin_coarse)

        except Exception as e:
            logger.error(f'Point {ID} could not be located in the coarser grid: {e}')

    # handle case where no polygons were generated
    if not polygons_coarse:
        logger.warning('No points could be located in the finer grid. Returning empty dataframes.')
        return gpd.GeoDataFrame(), gpd.GeoDataFrame()
        
    # concatenate polygons shapefile
    polygons_coarse = pd.concat(polygons_coarse)
    
    # regenerate point geometry to the new coordinates
    points_coarse = gpd.GeoDataFrame(
        points_coarse,
        geometry=gpd.points_from_xy(
            points_coarse[f'lon_{cfg.coarse_resolution}'],
            points_coarse[f'lat_{cfg.coarse_resolution}']
        ),
        crs=points_crs
    )
    points_coarse.sort_index(axis=1, inplace=True)
    
    # compute error
    points_coarse['abs_error'] = abs(points_coarse[f'area_{cfg.coarse_resolution}'] - points_coarse['area'])
    points_coarse['pct_error'] = points_coarse.abs_error / points_coarse['area'] * 100
    
    if save:
        save_results(cfg, points_coarse, polygons_coarse, cfg.coarse_resolution)

    return points_coarse, polygons_coarse