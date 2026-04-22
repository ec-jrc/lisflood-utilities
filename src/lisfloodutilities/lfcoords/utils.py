import logging
from typing import Tuple, Optional, Union
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from affine import Affine
from rasterio import features
from pyproj.crs.crs import CRS
import pyflwdir

from . import Config

# set logger
logger = logging.getLogger(__name__)


def check_points(
    cfg: Config,
    points: pd.DataFrame,
    extent: np.ndarray
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
    extent: tuple
        Extent of the flow direction map
        
    Returns:
    --------
    pandas.DataFrame
        The input table with points with conflicts removed.
    """
    
    # remove points with missing values
    mask_nan = points.isnull().any(axis=1)
    if mask_nan.sum() > 0:
        points = points.loc[~mask_nan].copy()
        logger.warning(f'{mask_nan.sum()} points were removed because of missing values.')
        
    # remove points with small catchment area
    mask_area = points['area'] < cfg.min_area
    if mask_area.sum() > 0:
        points = points.loc[~mask_area].copy()
        logger.info(f'{mask_area.sum()} points were removed due to their small catchment area.')
        
    # remove points outside the input flow direction map
    lon_min, lat_min, lon_max, lat_max = np.round(extent, 6)
    mask_lon = (points.lon < lon_min) | (points.lon > lon_max)
    mask_lat = (points.lat < lat_min) | (points.lat > lat_max)
    mask_extent = mask_lon | mask_lat
    if mask_extent.sum() > 0:
        points = points.loc[~mask_extent].copy()
        logger.info(f'{mask_extent.sum()} points were removed because they are outside the input flow direction map.')
        
    return points


def catchment_polygon(
    data: np.ndarray,
    transform: Affine,
    crs: Optional[CRS] = None,
    name: str = "catchment"
) -> gpd.GeoDataFrame:
    """
    Converts a boolean 2D array of catchment extent to a GeoDataFrame.
    
    Parameters:
    -----------
    data : numpy.ndarray
        The input raster data to be vectorized.
    transform : affine.Affine
        The affine transform matrix for the raster data.
    crs : dict or str
        The Coordinate Reference System (CRS) of the input data.
    name : str, optional
        The name to be used for the value column in the resulting GeoDataFrame.
    
    Returns:
    --------
    geopandas.GeoDataFrame
        A GeoDataFrame containing the vectorized geometries.
    """
    
    # Convert data to a dtype compatible with rasterio.features.shapes
    # (uint32 is not supported, so convert to uint16 if possible)
    if data.dtype == np.uint32:
        # Check if values fit in uint16
        if data.max() <= np.iinfo(np.uint16).max:
            data = data.astype(np.uint16)
        else:
            data = data.astype(np.uint8)
    
    # Generate shapes and associated values from the raster data
    feats_gen = features.shapes(
        data,
        mask=data != 0,
        transform=transform,
        connectivity=8,
    )
    
    # Create a list of features with geometry and properties
    feats = [
        {"geometry": geom, "properties": {name: val}} for geom, val in feats_gen
    ]
    
    # Check if no features were found
    if not feats:
        raise ValueError("No features found in the raster data.")
    
    # Create a GeoDataFrame from the features
    gdf = gpd.GeoDataFrame.from_features(feats, crs=crs)
    
    # Convert the value column to the same data type as the input raster data
    gdf[name] = gdf[name].astype(data.dtype)
    
    return gdf

        
def find_conflicts(
    points: gpd.GeoDataFrame,
    resolution: str,
    pct_error: float = 30, 
    save: Optional[Union[Path, str]] = None
) -> gpd.GeoDataFrame:
    """
    Finds conflicts in the new point layer, either due to points that overlap,
    or large catchment area errors.
    
    Parameters:
    -----------
    points: geopandas.GeoDataFrame
        Point layer resulting from `coordinates_fine` or `coordinates_coarse`
    resolution: str
        Spatial resolution of the fields in "points" to be checked. For instance, 
        '3min' will check the coordinates in the fields 'lat_3min' and 'lon_3min', 
        and the catchment area in the field 'area_3min'
    pct_error: float
        Maximum percentage error allowed in the derived catchment area compared 
        with the reference. It must be a value between 0 and 100
    save: pathlib.Path or string (optional)
        If provided, file name of the shapefile of conflicting points
        
    Returns:
    --------
    geopandas.GeoDataFrame
        Subset of "points" with conflicts. Only if "save" is None.
    """  
        
    # find overlaping points
    columns = [f'{col}_{resolution}' for col in ['lat', 'lon']]
    mask = points.duplicated(subset=columns, keep=False)
    duplicates = points[mask].copy()
    duplicates['conflict'] = 'points overlap'
    if duplicates.shape[0] > 0:
        n_duplicates = len(duplicates[columns[0]].unique())
        logger.warning(f'There are {n_duplicates} conflicts in which points are located at the same pixel in the finer grid')
        
    # errors in the delineated area
    assert 0 <= pct_error <= 100, '"pct_error" must be a value between 0 and 100'
    if 'pct_error' not in points.columns:
        points['pct_error'] = abs(points['area'] - points[f'area_{resolution}']) / points['area'] * 100
    mask = points.pct_error >= pct_error
    deviations = points[mask].copy()
    deviations['conflict'] = 'large area error'
    n_deviations = deviations.shape[0]
    if n_deviations > 0:
        logger.warning(f'There are {n_deviations} conflicts in which the new points area has a large error')
    
    # combine
    if (duplicates.shape[0] > 0) or (n_deviations > 0):
        conflicts = pd.concat((duplicates, deviations), axis=0)
        if save is not None:
            conflicts.to_file(save)
            logger.info(f'The conflicting points were saved in {save}')
        return conflicts

    return gpd.GeoDataFrame() # return an empty GeoDataFrame if no conflicts found


def create_input_file(geometry: gpd.GeoSeries, area: pd.Series, file: str):
    """
    Creates the input GeoJSON for the `basin-delineation` tool.

    Parameters:
    -----------
    geometry: geopandas.GeoSeries
        Coordinates of the points
    area: pd.Series
        Catchment area (km²)
    """
    
    # reformat input data
    gdf = gpd.GeoDataFrame(
        data={
            'lat': geometry.y,
            'lon': geometry.x,
            'area': area
        },
        geometry=geometry,
        crs=geometry.crs
        )
    gdf.index.name = 'ID'
    gdf.dropna(axis=0, how='any', inplace=True)

    # export
    gdf.to_file(file, driver='GeoJSON')


def define_neighbourhood(
        flwdir: pyflwdir.FlwdirRaster,
        x: float,
        y: float,
        n_cell: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Defines the indices of the pixels at a given distance range from the input coordinates.

    Parameters:
    -----------
    flwdir : pyflwdir.FlwdirRaster
        Flow direction raster object used for basin delineation.
    x: float
        Longitude or X coordinate of the input point.
    y: float
        Latitude or Y coordinate of the input point.
    n_cell: integer
        Distance (in number of cells) to explore around the input point.
    
    Returns:
    --------
    idxs : np.ndarray
        Flat array of indices.
    xs: np.ndarray
        Flat array of X coordinates.
    ys: np.ndarray
        Flat array of Y coordinates.
    """

    # find pixel in the grid associated to the input coords
    idx_o = flwdir.index(xs=x, ys=y)
    x_o, y_o = flwdir.xy(idx_o)

    # define search range
    cellsize = abs(flwdir.transform[0])
    range_xy = np.arange(-n_cell, n_cell + 1) * cellsize

    # define pixels to explore
    xs, ys = np.meshgrid(x_o + range_xy, y_o + range_xy)
    xs, ys = xs.flatten(), ys.flatten()
    idxs = flwdir.index(xs=xs, ys=ys)

    return idxs, xs, ys


def compute_area_error(
        uparea: np.ndarray,
        idxs: np.ndarray,
        area_ref: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the absolute difference between the delineated upstream area 
    and a reference area for a set of given pixel indices.

    Parameters:
    -----------
    uparea: nd.ndarray
        Upstream area raster in km².
    idxs : np.ndarray
        Array of flat indices representing the candidate outlet locations.
    area_ref : float
        Reference catchment area to compare against (in km²).

    Returns:
    --------
    area : np.ndarray
        Delineated upstream area for each index, rounded to 1 decimal place (km²).
    area_error : np.ndarray
        Absolute relative error between the delineated area and the reference area (km²).
    """

    # find pixel with closest upstream area
    area = uparea.ravel()[idxs].round(1)
    area_error = abs(area - area_ref) / area_ref

    return area, area_error


def compute_distance(
        n_cell: int,
        scaler: float = 0.92
) -> np.ndarray:
    """Computes distances (in number of cells) in a neighbourhood window.

    Parameters:
    -----------
    n_cell: integer
        Number of cells that define the neighbourhood window.
    scaler: float
        Scaling factor.

    Returns:
    --------
    distance: np.ndarray
        Flat array of distances in number of cells.
    """

    range = np.arange(-n_cell, n_cell + 1)
    distance = np.sqrt(range[np.newaxis, :]**2 + range[:, np.newaxis]**2) * scaler

    return distance.flatten()


def compute_shape_error(
        flwdir: pyflwdir.FlwdirRaster,
        idxs: np.ndarray,
        polygon_ref: gpd.GeoDataFrame,
) -> Tuple[gpd.GeoDataFrame, np.ndarray]:
    """
    Evaluates the spatial fit of delineated basins compared to a reference 
    polygon using the Jaccard Index (Intersection over Union) logic.

    The error is calculated as 1 minus the ratio of the intersection area to 
    the union area. A result of 0 indicates a perfect spatial match.

    Parameters:
    -----------
    flwdir : pyflwdir.FlwdirRaster
        Flow direction raster object used for basin delineation.
    idxs : np.ndarray
        Flat indices of the outlet points to be delineated and evaluated.
    polygon_ref : gpd.GeoDataFrame
        The "truth" or reference catchment boundary.

    Returns:
    --------
    basins : gpd.GeoDataFrame
        Concatenated table containing the delineated basin geometries.
    shape_error : np.ndarray
        The shape mismatch value (1 - IoU) for each outlet. 
        Range: [0, 1], where 0 is perfect overlap.
    """
    
    # the Mollweide projection will be used to compute area
    proj = 'ESRI:54009'
    
    basins, shape_error = [], []
    for idx in idxs:
        # delineate basin
        basin = catchment_polygon(
            flwdir.basins(idxs=idx),
            transform=flwdir.transform,
            crs=polygon_ref.crs,
            name='ID'
        )
        basins.append(basin)

        # intersection between reference and delineation
        intersection = gpd.overlay(polygon_ref, basin, how='intersection')
        intersection['area'] = intersection.to_crs(proj).area * 1e-6 # km²

        # union between reference and delineation
        union = gpd.overlay(polygon_ref, basin, how='union')
        union['area'] = union.to_crs(proj).area * 1e-6 # km²

        # compute intersection/union ratio
        shape_error.append(intersection['area'].sum() / union['area'].sum())
    basins = pd.concat(basins)
    shape_error = 1 - np.array(shape_error)

    return basins, shape_error