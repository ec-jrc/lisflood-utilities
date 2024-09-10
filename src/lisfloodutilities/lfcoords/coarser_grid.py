import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import rioxarray
import pyflwdir
from pathlib import Path
from tqdm import tqdm
from typing import Optional, Union
import logging
import warnings
warnings.filterwarnings("ignore")

from lisfloodpreprocessing import Config
from lisfloodpreprocessing.utils import catchment_polygon, downstream_pixel

# set logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def coordinates_coarse(
    cfg: Config,
    points_fine: Union[pd.DataFrame, gpd.GeoDataFrame],
    polygons_fine: gpd.GeoDataFrame,
    ldd_coarse: xr.DataArray,
    upstream_coarse: xr.DataArray,
    reservoirs: bool = False,
    save: bool = False
) -> Optional[gpd.GeoDataFrame]:
    """
    Transforms point coordinates from a high-resolution grid to a corresponding location in a coarser grid, aiming to match the shape of the catchment area derived from the high-resolution map. It updates the station coordinates and exports the catchment areas as shapefiles in the coarser grid.

    The function reads the upstream area map and local drainage direction (LDD) map in the coarse grid. It then finds the pixel in the coarse grid that best matches the catchment shape derived from the high-resolution map. The match is evaluated based on the intersection-over-union of catchment shapes and the ratio of upstream areas.

    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in the configuration file.
    points_fine: pandas.DataFrame or geopandas.GeoDataFrame
        A table with updated station coordinates and upstream areas in the finer grid. It's the result of coarse_grid.coarse_grid()`
    polygons_fine: geopandas.GeoDataFrame
        A table with the catchment polygons in the finer grid. It's the result of `coarse_grid.coarse_grid()`
    ldd_coarse: xarray.DataArray
        Map of local drainaige directions in the coarse grid
    upstream_coarse: xarray.DataArray
        Map of upstream area (m2) in the coarse grid
    reservoirs: bool
        Whether the points are reservoirs or not. If True, the resulting coordinates refer to one pixel downstream of the actual solution, a deviation required by the LISFLOOD reservoir simulation
    save: boolean
        If True, the updated table of points is exported as a shapefile.

    Returns:
    --------
    points: geopandas.GeoDataFrame
        A table with updated station coordinates and upstream areas in the coarser grid.
    polygons: geopandas.GeoDataFrame
        A table with the catchment polygons in the coarser grid.
    """
    
    points = points_fine.copy()
    
    # create river network
    fdir_coarse = pyflwdir.from_array(ldd_coarse.data,
                                      ftype='ldd',
                                      transform=ldd_coarse.rio.transform(),
                                      check_ftype=False,
                                      latlon=True)

    # boundaries and resolution of the input maps
    lon_min, lat_min, lon_max, lat_max = np.round(ldd_coarse.rio.bounds(), 6)
    cellsize = np.round(np.mean(np.diff(ldd_coarse.x)), 6) # degrees

    # extract resolution of the finer grid from 'points'
    cols_fine = [f'{col}_{cfg.FINE_RESOLUTION}' for col in ['area', 'lat', 'lon']]

    # add new columns to 'points'
    cols_coarse = [f'{col}_{cfg.COARSE_RESOLUTION}' for col in ['area', 'lat', 'lon']]
    points[cols_coarse] = np.nan

    # search range of 5x5 array -> this is where the best point can be found in the coarse grid
    rangexy = np.linspace(-2, 2, 5) * cellsize # arcmin
    polygons_coarse = []
    for ID, attrs in tqdm(points.iterrows(), total=points.shape[0], desc='points'):

        # real upstream area
        area_ref = attrs['area']

        # coordinates and upstream area in the fine grid
        lat_fine, lon_fine, area_fine = attrs[[f'{col}_{cfg.FINE_RESOLUTION}' for col in ['lat', 'lon', 'area']]]

        if (area_ref < cfg.MIN_AREA) or (area_fine < cfg.MIN_AREA):
            logger.warning(f'The catchment area of station {ID} is smaller than the minimum of {cfg.MIN_AREA} km2')
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
                intersection = gpd.overlay(polygons_fine.loc[[ID]], basin, how='intersection')
                union = gpd.overlay(polygons_fine.loc[[ID]], basin, how='union')
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
        
        # save polygon
        polygons_coarse.append(basin_coarse)
        
        # move the result one pixel downstream, in case of reservoir
        if reservoirs:
            lat_coarse, lon_coarse = downstream_pixel(lat_coarse, lon_coarse, ldd_coarse)
            
        # update new columns in 'points'
        points.loc[ID, cols_coarse] = [int(area_coarse), round(lat_coarse, 6), round(lon_coarse, 6)]
        
    # concatenate polygons shapefile
    polygons_coarse = pd.concat(polygons_coarse)
    
    # convert points to geopandas
    geometry = [Point(xy) for xy in zip(points[f'lon_{cfg.COARSE_RESOLUTION}'], points[f'lat_{cfg.COARSE_RESOLUTION}'])]
    points = gpd.GeoDataFrame(points, geometry=geometry, crs=4326)
    points.sort_index(axis=1, inplace=True)
    
    if save is True:
        # polygons
        polygon_shp = cfg.OUTPUT_FOLDER / f'catchments_{cfg.COARSE_RESOLUTION}.shp'
        polygons_coarse.to_file(polygon_shp)
        logger.info(f'Catchments in the coarser grid have been exported to: {polygon_shp}')
        # points
        point_shp = cfg.OUTPUT_FOLDER / f'{cfg.POINTS.stem}_{cfg.COARSE_RESOLUTION}.shp'
        points.to_file(point_shp)
        logger.info(f'The updated points table in the coarser grid has been exported to: {point_shp}')

    return points, polygons_coarse