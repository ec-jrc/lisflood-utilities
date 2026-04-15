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
from lisfloodpreprocessing.utils import catchment_polygon

os.environ['USE_PYGEOS'] = '0'
warnings.filterwarnings("ignore")

# set logger
logger = logging.getLogger(__name__)


def coordinates_coarse(
    cfg: Config,
    points_fine: Union[pd.DataFrame, gpd.GeoDataFrame],
    polygons_fine: gpd.GeoDataFrame,
    ldd_coarse: xr.DataArray,
    upstream_coarse: xr.DataArray,
    save: bool = False
) -> Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]:
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
    ldd_coarse : xr.DataArray
        Map of local drainage directions in the coarse grid.
    upstream_coarse : xr.DataArray
        Map of upstream area (m2) in the coarse grid.
    save : bool, optional
        If True, the updated tables are exported as shapefiles.

    Returns
    -------
    Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]
        A tuple containing:
        - A table with updated station coordinates in the coarser grid.
        - A table with the catchment polygons in the coarser grid.
    """
    
    points_coarse = points_fine.copy()
    n_points = points_coarse.shape[0]
    
    # create river network
    fdir_coarse = pyflwdir.from_array(
        ldd_coarse.data,
        ftype='ldd',
        transform=ldd_coarse.rio.transform(),
        check_ftype=False,
        latlon=True
    )

    # get resolution of the coarse grid
    cellsize = np.round(np.mean(np.diff(ldd_coarse.x)), 6) # degrees

    # extract resolution of the finer grid
    cols = ['area', 'lat', 'lon']
    cols_fine = [f'{col}_{cfg.fine_resolution}' for col in cols]

    # add new columns to 'points_coarse'
    cols_coarse = [f'{col}_{cfg.coarse_resolution}' for col in cols]
    points_coarse[cols_coarse] = np.nan

    # search range of 5x5 array
    range_xy = np.linspace(-2, 2, 5) * cellsize # arcmin
    polygons_coarse = []
    for point_id, attrs in tqdm(points_coarse.iterrows(), total=n_points, desc='points'):
        try:
            # real upstream area
            area_ref, lat_ref, lon_ref = attrs[cols]

            # coordinates and upstream area in the fine grid
            area_fine, lat_fine, lon_fine = attrs[cols_fine]

            # find ratio
            logger.debug('Start search')
            inter_vs_union, area_ratio, area_lisf = [], [], []
            for delta_lat in range_xy:
                for delta_lon in range_xy:
                    lon = lon_fine + delta_lon
                    lat = lat_fine + delta_lat
                    basin = catchment_polygon(
                        fdir_coarse.basins(xy=(lon, lat)).astype(np.int32),
                        transform=ldd_coarse.rio.transform(),
                        crs=ldd_coarse.rio.crs
                    )

                    # calculate union and intersection of shapes
                    intersection = gpd.overlay(polygons_fine.loc[[point_id]], basin, how='intersection')
                    union = gpd.overlay(polygons_fine.loc[[point_id]], basin, how='union')
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
            i_centre = int(len(range_xy)**2 / 2) # middle point
            area_centre = area_lisf[i_centre]
            
            # use middle point if errors are small
            abs_error = abs(area_shape - area_centre)
            pct_error = 100 * abs(1 - area_centre / area_shape)
            if (abs_error <= cfg.abs_error) and (pct_error <= cfg.pct_error):
                i_shape = i_centre
                area_shape = area_centre
                
            # coordinates in the fine resolution
            i = i_shape // len(range_xy)
            j = i_shape % len(range_xy)
            lat = lat_fine + range_xy[i]
            lon = lon_fine + range_xy[j]

            # coordinates and upstream area on coarse resolution
            area = upstream_coarse.sel(x=lon, y=lat, method='nearest')
            area_coarse = area.item() * 1e-6
            lon_coarse = area.x.item()
            lat_coarse = area.y.item()

            # derive catchment polygon from the selected coordinates
            basin_coarse = catchment_polygon(
                fdir_coarse.basins(xy=(lon_coarse, lat_coarse)).astype(np.int32),
                transform=ldd_coarse.rio.transform(),
                crs=ldd_coarse.rio.crs,
                name='ID'
            )
            basin_coarse['ID'] = point_id
            basin_coarse.set_index('ID', inplace=True)
            basin_coarse[cols] = area_ref, lat_ref, lon_ref
            basin_coarse[cols_fine] = area_fine, lat_fine, lon_fine
            basin_coarse[cols_coarse] = area_coarse, lat_coarse, lon_coarse

            # save polygon
            polygons_coarse.append(basin_coarse)

            # update new columns in 'points_coarse'
            points_coarse.loc[point_id, cols_coarse] = [int(area_coarse), round(lat_coarse, 6), round(lon_coarse, 6)]

            # logger.info(f'Point {point_id} located in the coarser grid')
        except Exception as e:
            logger.error(f'Point {point_id} could not be located in the coarser grid: {e}')

    # handle case where no polygons were generated
    if not polygons_coarse:
        logger.warning('No points could be located in the finer grid. Returning empty dataframes.')
        return gpd.GeoDataFrame(), gpd.GeoDataFrame()
        
    # concatenate polygons shapefile
    polygons_coarse = pd.concat(polygons_coarse)
    
    # convert points to geopandas
    points_coarse = gpd.GeoDataFrame(
        points_coarse, 
        geometry=gpd.points_from_xy(
            points_coarse[f'lon_{cfg.coarse_resolution}'], 
            points_coarse[f'lat_{cfg.coarse_resolution}']
        ), 
        crs=ldd_coarse.rio.crs
    )
    points_coarse.sort_index(axis=1, inplace=True)
    
    # compute error
    points_coarse['abs_error'] = abs(points_coarse[f'area_{cfg.coarse_resolution}'] - points_coarse['area'])
    points_coarse['pct_error'] = points_coarse.abs_error / points_coarse['area'] * 100
    
    if save:
        # polygons
        polygon_shp = cfg.output_folder / f'catchments_{cfg.coarse_resolution}.shp'
        polygons_coarse.to_file(polygon_shp)
        logger.info(f'Catchments in the coarser grid have been exported to: {polygon_shp}')
        
        # points
        point_shp = cfg.output_folder / f'{cfg.points.stem}_{cfg.coarse_resolution}.shp'
        points_coarse.to_file(point_shp)
        logger.info(f'The updated points table in the coarser grid has been exported to: {point_shp}')

    return points_coarse, polygons_coarse
