import logging
from typing import Tuple
import warnings
from tqdm import tqdm

import numpy as np
import pandas as pd
import geopandas as gpd

from . import Config
from .io import read_flowdir, read_shapefiles, save_results
from .utils import check_points, catchment_polygon, define_neighbourhood, compute_area_error, compute_distance


logger = logging.getLogger(__name__)


def coordinates_fine(
    cfg: Config,
    save: bool = False
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Processes point coordinates to find the most accurate pixel in a high-resolution
    map, based on a reference value of catchment area. It updates the station
    coordinates and exports the catchment areas as shapefiles.

    Parameters
    ----------
    cfg : Config
        Configuration object containing file paths and parameters.
    save : bool, optional
        If True, the updated table of points and catchments are exported.

    Returns
    -------
    Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        A tuple containing:
        - A table with updated station coordinates and upstream areas.
        - A table with the catchment polygons in the finer grid.
        Returns empty GeoDataFrames if no polygons are generated.
    """
    
    # read flow direction map
    try:
        flwdir, resolution = read_flowdir(cfg.flwdir_fine, ftype='d8') # type: ignore
        cfg.fine_resolution = resolution
        logger.info(f'Flow direction map at {resolution} resolution correctly read: {cfg.flwdir_fine}')
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        raise RuntimeError("Aborting computation due to missing or corrupt flow direction map.") from e

    # read points shapefile
    try:
        points = read_shapefiles(cfg.points) # type: ignore
        logger.info(f'Shapefile of points correctly read: {cfg.points}')
        points = check_points(cfg, points, flwdir.extent)
        n_points = len(points)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        raise RuntimeError("Aborting computation due to missing or corrupt points shapefile.") from e

    # Compute upstream area map in km²
    uparea = flwdir.upstream_area(unit='km2')

    # create output point layer
    points = points.copy()
    cols = ['area', 'lat', 'lon']
    cols_fine = [f'{col}_{resolution}' for col in cols]
    points[cols_fine] = np.nan

    polygons = []
    for ID, point in tqdm(points.iterrows(), total=n_points, desc='points'):
        try:
            # real upstream area and coordinates
            area_ref, lat_ref, lon_ref = point[cols]

            # search new coordinates in an increasing range
            ranges = [55, 101, 151] # number of cells
            penalties = [500, 500, 1000] # number of cells
            dist_factors = [2, 0.5, 0.25] # multiplier of the distanceS
            max_errors = [50, 80, np.nan] # allowed absolute error (%)
            
            # initialize variables in case no match is found
            area_fine = lat_fine = lon_fine = np.nan
            idxs = None
            i = None
            
            for n_cell, penalty, factor, error_thr in zip(ranges, penalties, dist_factors, max_errors):
                # define pixels to explore
                idxs, xs, ys = define_neighbourhood(flwdir, lon_ref, lat_ref, n_cell)

                # compute error in catchment area
                area, area_error = compute_area_error(uparea, idxs, area_ref)
                area_error *= 100 # %

                # compute distance
                distance = compute_distance(n_cell)

                # penalise if error is too big
                distance = np.where(area_error <= error_thr, distance, distance + penalty)

                # update error based on distance
                area_error += factor * distance

                if np.min(area_error) <= error_thr:
                    # select best pixel and extract coordinates and area
                    i = np.argmin(area_error)
                    lon_fine, lat_fine, area_fine = round(xs[i].item(), 6), round(ys[i].item(), 6), round(area[i].item(), 3)
                    break

            # update new columns in 'points'
            points.loc[ID, cols_fine] = [area_fine, lat_fine, lon_fine]

            # delineate basin only if a valid pixel was found
            if idxs is not None and i is not None:
                basin_fine = catchment_polygon(
                    flwdir.basins(idxs=idxs[i]),
                    transform=flwdir.transform,
                    crs=points.crs,
                    name='ID'
                )
                basin_fine['ID'] = ID
                basin_fine.set_index('ID', inplace=True)
                basin_fine[cols] = point[cols].values
                basin_fine[cols_fine] = area_fine, lat_fine, lon_fine

                # save polygon
                polygons.append(basin_fine)

        except Exception as e:
            logger.error(f'Point {ID} could not be located in the finer grid: {e}')

    # handle case where no polygons were generated
    if not polygons:
        logger.warning('No points could be located in the finer grid. Returning empty dataframes.')
        return gpd.GeoDataFrame(), gpd.GeoDataFrame()
        
    # concatenate polygons shapefile
    polygons = pd.concat(polygons)
    
    # regenerate point geometry to the new coordinates
    points = gpd.GeoDataFrame(
        points, 
        geometry=gpd.points_from_xy(
            points[f'lon_{resolution}'], 
            points[f'lat_{resolution}']
            ), 
        crs=points.crs
    )
    points.sort_index(axis=1, inplace=True)
    
    # compute error
    points['abs_error'] = abs(points[f'area_{resolution}'] - points['area'])
    points['pct_error'] = points.abs_error / points['area'] * 100
    
    if save:
        save_results(cfg, points, polygons, resolution)

    return points, polygons