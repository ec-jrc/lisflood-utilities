import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import geopandas as gpd
import xarray as xr
from affine import Affine
from rasterio import features
from pyproj.crs import CRS
from typing import Tuple



def find_pixel(
    upstream: xr.DataArray,
    lat: int,
    lon: int,
    area: float,
    rangexy: int = 55,
    penalty: int = 500,
    factor: int = 2,
    distance_scaler: float = .92,
    error_threshold: int = 50
) -> Tuple:
    """
    Find the coordinates of the pixel in the upstream map with a smaller error compared with a reference area.
    
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
    rangexy: int, optional
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
    lat_new : float
        The latitude of the new location.
    lon_new : float
        The longitude of the new location.
    min_error : float
        The minimum error value at the new location.
    """

    # find coordinates of the nearest pixel in the map
    nearest_pixel = upstream.sel(y=lat, x=lon, method='nearest')
    lat_orig, lon_orig = [nearest_pixel[coord].item() for coord in ['y', 'x']]

    # extract subset of the upstream map
    cellsize = np.mean(np.diff(upstream.x.data))
    delta = rangexy * cellsize + 1e-6
    upstream_sel = upstream.sel(y=slice(lat_orig + delta, lat_orig - delta),
                                x=slice(lon_orig - delta, lon_orig + delta))
    
    # distance from the original pixel (in pixels)
    i = np.arange(-rangexy, rangexy + 1)
    ii, jj = np.meshgrid(i, i)
    distance = xr.DataArray(data=np.sqrt(ii**2 + jj**2) * distance_scaler, coords=upstream_sel.coords, dims=upstream_sel.dims)

    # percent error in catchment area
    error = 100 * abs(area - upstream_sel) / area

    # penalise if error is too big
    distance = distance.where(error <= error_threshold, distance + penalty)

    # update error based on distance (every pixel deviation is a 'factor' increase in the error)
    error += factor * distance

    # the new location is that with the smallest error
    min_error = error.where(error == error.min(), drop=True)
    lat_new, lon_new = [min_error[coord].item() for coord in ['y', 'x']]
    
    return lat_new, lon_new, min_error.item()



def catchment_polygon(
    data: np.ndarray,
    transform: Affine,
    crs: CRS,
    name: str = "catchment"
) -> gpd.GeoDataFrame:
    """
    Convert a boolean 2D array of catchment extent to a GeoDataFrame.
    
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
    gdf : geopandas.GeoDataFrame
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