import os
import logging
from typing import Tuple, Optional, Union
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from affine import Affine
from rasterio import features
from pyproj.crs import CRS

os.environ['USE_PYGEOS'] = '0'

# set logger
logger = logging.getLogger(__name__)


def find_pixel(
    upstream: xr.DataArray,
    lat: int,
    lon: int,
    area: float,
    range_xy: int = 55,
    penalty: int = 500,
    factor: int = 2,
    distance_scaler: float = .92,
    error_threshold: int = 50
) -> Tuple[float, float, float]:
    """
    Finds the coordinates of the pixel in the upstream map with a smaller 
    error compared with a reference area.
    
    Parameters:
    -----------
    upstream: xr.DataArray
        The upstream data containing latitude and longitude coordinates.
    lat: float
        The original latitude value.
    lon: float
        The original longitude value.
    area: float
        The reference area to calculate percent error.
    range_xy: int, optional
        The range in both x and y directions to search for the new location.
    penalty: int, optional
        The penalty value to add to the distance when the percent error is too high.
    error_threshold: float, optional
        The threshold for the percent error to apply the penalty.
    factor: int, optional
        The factor to multiply with the distance for the error calculation.
    distance_scaler: float, optional
        The scaling factor for the distance calculation in pixels.
    
    Returns:
    --------
    Tuple[float, float, float]
        A tuple containing:
        - lat_new : float
            The latitude of the new location.
        - lon_new : float
            The longitude of the new location.
        - min_error : float
            The minimum error value at the new location.
    """

    # find coordinates of the nearest pixel in the map
    nearest_pixel = upstream.sel(y=lat, x=lon, method='nearest')
    lat_orig, lon_orig = (nearest_pixel[coord].item() for coord in ['y', 'x'])

    # extract subset of the upstream map
    cellsize = np.mean(np.diff(upstream.x.data))
    delta = range_xy * cellsize + 1e-6
    upstream_sel = upstream.sel(y=slice(lat_orig + delta, lat_orig - delta),
                                x=slice(lon_orig - delta, lon_orig + delta))
    
    # distance from the original pixel (in pixels)
    i = np.arange(-range_xy, range_xy + 1)
    ii, jj = np.meshgrid(i, i)
    distance = xr.DataArray(
        data=np.sqrt(ii**2 + jj**2) * distance_scaler, 
        coords=upstream_sel.coords, 
        dims=upstream_sel.dims
    )

    # percent error in catchment area
    error = 100 * abs(area - upstream_sel) / area

    # penalise if error is too big
    distance = distance.where(error <= error_threshold, distance + penalty)

    # update error based on distance
    error += factor * distance

    # the new location is that with the smallest error
    min_error = error.where(error == error.min(), drop=True)
    if min_error.size == 1:
        lat_new, lon_new = [min_error[coord].item() for coord in ['y', 'x']]
    else:
        idx = np.argwhere(~np.isnan(min_error.data.flatten()))[0][0]
        x_idx, y_idx = np.unravel_index(idx, min_error.shape)
        lat_new, lon_new = min_error.y.data[y_idx], min_error.x.data[x_idx]
    
    return lat_new, lon_new, min_error.item()


def catchment_polygon(
    data: np.ndarray,
    transform: Affine,
    crs: CRS,
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