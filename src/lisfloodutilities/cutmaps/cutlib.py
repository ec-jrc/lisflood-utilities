"""
Copyright 2019-2026 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import os
import sys
import datetime
import warnings
from pathlib import Path
from typing import List, Tuple, Union, Optional

import xarray as xr
import numpy as np

from lisfloodutilities.readers.pcr import PCRasterMap

from .helpers import (col2netcdf, array_to_nc_from_clone, bbox_from_netcdf,
                      get_river_network_from_map, LabelledProgressBar,
                      COORDINATE_NAMES, LATITUDE_NAMES, LATITUDE_NAME_PAIR)
from .. import version, logger

from earthkit.hydro import catchments

# Performance: Reduced compression level from 4 to 1 for faster writes
# Compression can be increased later if file size is a concern
encoding_netcdf_vars = {'zlib': True, 'complevel': 4}

# Default chunk size for time dimension processing (number of timesteps per chunk)
DEFAULT_TIME_CHUNK_SIZE = 100

# Value used to identify the masked cells
MASK_VALUE = 1

# Internal output filenames
FULL_MASK_FILENAME = 'mask_full.nc'
SMALL_MASK_FILENAME = 'mask_small.nc'
OUTLETS_FILENAME = 'outlets.nc'


def cutmap(f, fileout, x_min, x_max, y_min, y_max, use_coords = True):
    """
    Cut a map from a NetCDF file using coordinates or indices.

    Parameters
    ----------
    f : str
        Input NetCDF file path
    fileout : str
        Output NetCDF file path
    x_min, x_max, y_min, y_max : float
        Bounding box coordinates or indices
    use_coords : bool, optional
        If True, use coordinates; if False, use indices (default: True)

    Returns
    -------
    None
        Writes the cropped NetCDF file to disk
    """
    logger.info(f"## Starting to cut map from {f} to {fileout}")
    if not os.path.exists(f):
        raise IOError(f'File not found: {f}')
    nc, num_dims = open_dataset(f)
    var = str([v for v in nc.variables if len(nc.variables[v].dims) == num_dims][0])
    logger.info('Variable: %s', var)
    logger.info('Cropping %s from %s to %s', var, (x_min, y_min), (x_max, y_max))
    if 'time' in nc.dims:
        logger.info('Time dimension found: %s', nc.sizes['time'])
    logger.info('Using coordinates: %s', use_coords)
    sliced_var = None
    if use_coords:
        if any(isinstance(val, int) for val in [x_min, x_max, y_min, y_max]):
            raise ValueError('box values must be float when using cut_from_coords')
        sliced_var = cut_from_coords(nc, var, x_min, x_max, y_min, y_max)
    elif any(isinstance(val, float) for val in [x_min, x_max, y_min, y_max]):
        raise ValueError('box values must be integer when using cut_indices')
    else:
        # user provides with indices directly (not coordinates)
        sliced_var = cut_from_indices(nc, var, x_min, x_max, y_min, y_max)

    if sliced_var is not None:
        if 'missing_value' in sliced_var.encoding:
            sliced_var.encoding['_FillValue'] = sliced_var.encoding['missing_value']
        logger.info('Creating: %s', fileout)
        # encoding_netcdf_vars['scale_factor'] = sliced_var.attrs.get('scale_factor')
        # encoding_netcdf_vars['add_offset'] = sliced_var.attrs.get('add_offset')
        delayed_obj = sliced_var.to_netcdf(fileout, compute=False,
                                           encoding={var: encoding_netcdf_vars})
        if hasattr(delayed_obj, 'compute'):
            with LabelledProgressBar(label=os.path.basename(fileout), dt=0.1):
                _ = delayed_obj.compute()

    grid_mapping = sliced_var.attrs.get('grid_mapping')
    if grid_mapping in nc.variables:
        varname = grid_mapping
        varproj = nc.variables[varname]
        logger.info('Writing projection variable: %s - %s', varname, varproj.attrs)
        # Write CRS as a scalar variable to avoid duplicate dimension name warnings
        # CRS/projection variables are metadata-only (no meaningful array data)
        dims = varproj.dims
        logger.info('Writing projection variable with dims: %s', dims)
        # Handle potential duplicate dimension names in CRS variables 
        if len(dims) > 1:
            unique_dims = list(dict.fromkeys(dims))  # Preserves order while removing duplicates
            if len(unique_dims) < len(dims):
                dims = tuple(unique_dims)
                logger.info('Deduplicated CRS dimensions to: %s', dims)
        # Handle case where CRS variable has repeated dimension names (GDAL/CDO artifacts)
        if len(dims) > 1 and dims[0] == dims[1]:
            dims = (dims[0],)
            logger.info('Fixed duplicate CRS dims to: %s', dims)
        # Write the projection variable with proper dimensions
        proj_data = varproj.data
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Duplicate dimension names present",
                category=UserWarning,
            )
            del_res = xr.DataArray(name=varname, data=proj_data, dims=dims,
                                   attrs=varproj.attrs).to_netcdf(fileout, mode='a', compute=False)
        if hasattr(del_res, 'compute'):
            with LabelledProgressBar(label=os.path.basename(fileout), dt=0.1):
                _ = del_res.compute()

    # adding global attributes only if the file contains variables
    try:
        # Open the newly created file in append mode
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Duplicate dimension names present",
                category=UserWarning,
            )
            with xr.open_dataset(fileout, mode='a', decode_cf=False) as nc_out:
                if nc_out.data_vars:
                    nc_out.attrs = nc.attrs
                    nc_out.attrs['history'] = 'lisfloodutilities cutmaps {} {} \n {}'.format(
                        version, datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), nc_out.attrs.get('history', '')
                    )
                    nc_out.attrs['conventions'] = 'CF-1.6'
                    nc_out.attrs['institution'] = 'JRC E1'
                    nc_out.attrs['source_software'] = 'lisfloodutilities cutmaps {}'.format(version)
                    nc_out.attrs.pop('Source_Software', None)
                    nc_out.attrs.pop('Institution', None)
                    nc_out.attrs.pop('Conventions', None)
                    logger.info('Writing additional attrs to: %s - %s', fileout, nc_out.attrs)
                    del_res = nc_out.to_netcdf(fileout, mode='a', compute=False, encoding={var: encoding_netcdf_vars})
                    if hasattr(del_res, 'compute'):
                        with LabelledProgressBar(label=os.path.basename(fileout), dt=0.1):
                            _ = del_res.compute()
    except Exception as e:
        logger.warning('Cannot add global attributes to %s - %s', fileout, e)
    finally:
        nc.close()
    logger.info(f"## Finished processing variable '{var}' and saved to {fileout}")


def open_dataset(file_path: Union[Path, str], time_chunk_size: int = DEFAULT_TIME_CHUNK_SIZE) -> Tuple[xr.Dataset, int]:
    """
    Open a NetCDF dataset with optimized chunking for time dimension.
    
    Parameters
    ----------
    file_path : Union[Path, str]
        Path to the NetCDF file
    time_chunk_size : int
        Number of timesteps per chunk (default: 100)
        
    Returns
    -------
    Tuple[xr.Dataset, int]
        Tuple of (opened dataset, number of dimensions)
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Duplicate dimension names present",
            category=UserWarning,
        )
        try:
            # Performance: Use explicit chunk size for better memory management
            # 'auto' can sometimes create chunks that are too large
            nc = xr.open_dataset(file_path, chunks={'time': time_chunk_size}, decode_cf=False)
            if 'time' in nc.coords:
                num_dims = 3
            else:
                num_dims = 2
        except Exception:  # file has no time component
            num_dims = 2
            nc = xr.open_dataset(file_path, decode_cf=False)
    return nc, num_dims


def cut_time_chunks(
    nc: xr.Dataset,
    var: str,
    y_indices,
    x_indices,
    chunk_size: int = DEFAULT_TIME_CHUNK_SIZE
) -> xr.DataArray:
    """
    Process time dimension in chunks for better memory efficiency.
    
    This function is useful for large time-series files where loading all
    timesteps at once may cause memory issues. It processes the data in
    chunks and concatenates the results.
    
    Parameters
    ----------
    nc : xr.Dataset
        The xarray dataset containing the data
    var : str
        Variable name to slice
    y_indices
        Slice or array of indices for the y (latitude) dimension
    x_indices
        Slice or array of indices for the x (longitude) dimension
    chunk_size : int
        Number of timesteps to process at once (default: 100)
        
    Returns
    -------
    xr.DataArray
        The sliced data array with all timesteps
    """
    if 'time' not in nc.variables:
        # No time dimension, return regular slice
        return nc[var][y_indices, x_indices]
    
    # Use 'time' directly as the dimension name, not nc.dims.get() which returns sizes
    time_dim_name = 'time'
    if time_dim_name not in nc.sizes:
        return nc[var][y_indices, x_indices]
    
    total_timesteps = nc.sizes[time_dim_name]
    
    # For small files, process all at once
    if total_timesteps <= chunk_size:
        return nc[var][:, y_indices, x_indices]
    
    # Process in chunks
    chunks = []
    for start in range(0, total_timesteps, chunk_size):
        end = min(start + chunk_size, total_timesteps)
        chunk = nc[var][start:end, y_indices, x_indices]
        chunks.append(chunk)
    
    # Concatenate chunks along time dimension
    result = xr.concat(chunks, dim='time')
    # Ensure we return a DataArray (xr.concat can return DataTree in some cases)
    if isinstance(result, xr.DataTree):
        result = result.to_dataset().get(var, result)
    return result  # type: ignore[return-value]


def cut_from_indices(nc: xr.Dataset, var: str, x_min: float, x_max: float, y_min: float, y_max: float) -> xr.DataArray:
    """
    Slice a variable from a dataset using indices.
    
    Note: netcdf has lats on first dimension e.g. y_min:y_max are Y/lat dimension indices
    that in nc file are stored on first dimension: ta(time, lat, lon)
    you can always adjust indices in input in order to match your nc files structure
    
    Performance: Uses chunk-based processing for time dimension to reduce memory usage.
    """
    if 'time' in nc.variables:
        # Use chunk-based processing for better memory efficiency
        y_slice = slice(y_min, y_max + 1)
        x_slice = slice(int(x_min), int(x_max) + 1)
        sliced_var = cut_time_chunks(nc, var, y_slice, x_slice)
    else:
        sliced_var = nc[var][y_min:y_max + 1, x_min:x_max + 1]
    return sliced_var


def cut_from_coords(nc: xr.Dataset, var: str, x_min: float, x_max: float, y_min: float, y_max: float) -> xr.DataArray:
    """
    Slice a variable from a dataset using coordinate bounds.
    
    Performance: Uses chunk-based processing for time dimension to reduce memory usage.
    """
    lats = None
    lons = None
    # Find the latitude and longitude variables
    for lat_var_name in LATITUDE_NAME_PAIR:
        if lat_var_name in nc.variables:
            lats = nc.variables[lat_var_name][:]
            lons = nc.variables[LATITUDE_NAME_PAIR[lat_var_name]][:]
            break
    if lats is None or lons is None:
        raise ValueError(f'Latitude or Y should be one of: {LATITUDE_NAMES}')
    # find indices
    # np float values comparisons are very sensitive to machine representation
    # so we add some "space" around bounding box when coordinates matches exactly
    buffer_y = abs(lats[0] - lats[1]) / 1000
    buffer_x = abs(lons[0] - lons[1]) / 1000
    # Handle both ascending and descending coordinate orders
    y_min_bound = min(y_min, y_max)
    y_max_bound = max(y_min, y_max)
    x_min_bound = min(x_min, x_max)
    x_max_bound = max(x_min, x_max)
    ys = np.where((lats > y_min_bound - buffer_y) & (lats < y_max_bound + buffer_y))[0]
    xs = np.where((lons > x_min_bound - buffer_x) & (lons < x_max_bound + buffer_x))[0]
    
    if 'time' in nc.variables:
        # Use chunk-based processing for better memory efficiency
        sliced_var = cut_time_chunks(nc, var, ys, xs)
    else:
        try:
            sliced_var = nc[var][ys, xs]
        except IndexError:
            # it happens when cutting lat or lon netcdf files (not 2D array)
            if var not in COORDINATE_NAMES:
                raise ValueError('Cannot find main variable to cut')

            sliced_var = nc[var][ys] if var in LATITUDE_NAMES else nc[var][xs]
    return sliced_var


def get_filelist(input_folder: str = '', static_data_folder: str = '', input_file: str = '') -> List[Path]:
    list_to_cut = []
    if input_folder and len(input_folder) > 0:
        list_to_cut = [f for f in Path(input_folder).glob('**/*.nc')]
    elif static_data_folder and len(static_data_folder) > 0:
        list_to_cut = [f for f in Path(static_data_folder).glob('**/*') if '/.git/' not in f.as_posix()]
    if input_file and len(input_file) > 0:
        list_to_cut = [Path(input_file)]
    logger.info('==================> Going to cut %d files', len(list_to_cut))
    return list_to_cut


def get_cuts(cuts=None, cuts_indices=None, mask=None):
    if mask:
        if not os.path.isfile(mask):
            raise FileNotFoundError('Wrong input mask: %s not a file' % mask)
        maskname, ext = os.path.splitext(mask)
        if ext == '.map':
            mask = PCRasterMap(mask)
            lats = mask.lats
            lons = mask.lons
            mask.close()
            x_min, x_max = float(np.min(lons)), float(np.max(lons))
            y_min, y_max = float(np.min(lats)), float(np.max(lats))
        elif ext == '.nc':
            x_min, x_max, y_min, y_max = bbox_from_netcdf(mask)
        else:
            logger.error('Mask map format not recognized. Must be either .map or .nc. Found %s', ext)
            sys.exit(1)
        logger.info('MASK: \nmin x: %s \nmax x: %s \nmin y: %s \nmax y: %s', x_min, x_max, y_min, y_max)
    elif cuts:
        # user provided coordinates bounds
        x_min, x_max, y_min, y_max = cuts
        logger.info('CUTS: \nmin x: %s \nmax x: %s \nmin y: %s \nmax y: %s', x_min, x_max, y_min, y_max)
    elif cuts_indices:
        # user provided indices bounds
        x_min, x_max, y_min, y_max = cuts_indices
        logger.info('CUTS_INDICES: \nmin x: %s \nmax x: %s \nmin y: %s \nmax y: %s', x_min, x_max, y_min, y_max)
    else:
        logger.error('You must provide either cuts (in the format "lonmin lonmax latmin latmax") or cuts_indices (in the format "imin imax jmin jmax") or a mask map')
        sys.exit(1)
    return x_min, x_max, y_min, y_max


def mask_from_ldd(ldd_map: Union[Path, str], stations: Union[Path, str]) -> Tuple[Path, Path, Path]:
    """
    Generate a mask map from a LDD where the outlets are identified in the stations file
    """
    ldd_map_path = Path(ldd_map) if isinstance(ldd_map, str) else ldd_map
    stations_path = Path(stations) if isinstance(stations, str) else stations
    path = stations_path.parent
    masknc_path = Path(path, FULL_MASK_FILENAME)
    outlets_path = Path(path, OUTLETS_FILENAME)
    smallmask_path = Path(path, SMALL_MASK_FILENAME)
    # clean existing files from previuos executions
    for out_file in (masknc_path, smallmask_path, outlets_path):
        if out_file.exists():
            os.unlink(out_file)

    # Default format for output netcdf file is NETCDF3_CLASSIC
    metadata = {'variable': {'description': 'stations id', 'longname': 'platform_id', 'units': '',
                             'shortname': 'outlets', 'mv': 0},
                'source': 'JRC E.1 Space, Security, Migration',
                'reference': 'JRC E.1 Space, Security, Migration',
                'geographical': {'datum': ''}}

    # Identify the outlets on the ldd map marked by the station coordinates
    outlets_raster = col2netcdf(stations_path, outlets_path, ldd_map, metadata, quiet=False)

    # Get outlets as row, col indexes of elements greater than 0
    # Note: The fill value is -1 (np.iinfo(np.int32).min), so we use > 0 to filter
    # earthkit expects (row, col) format, not (col, row)
    rows, cols = np.where(outlets_raster > 0)
    outlets = np.column_stack((rows, cols))  # (row, col) format for earthkit

    # Obtain the catchments that contain these outlets creating a boolean mask where 1 identifies
    # the catchment cells
    network = get_river_network_from_map(ldd_map_path)
    catchments_mask = catchments.find(network, outlets)
    # Convert xarray DataArray to numpy array if needed
    if hasattr(catchments_mask, 'values'):
        catchments_mask = catchments_mask.values

    # Setup the mask values to return 1 as masked and np.nan as unmasked
    # (following the convention of mask maps used in cutmaps)
    catchments_mask[catchments_mask<0] = np.nan
    catchments_mask[catchments_mask>=0] = MASK_VALUE

    # Mask map for netCDF format
    nc_metadata = {'variable': {'description': 'Mask Area', 'longname': 'area', 'units': '',
                                'shortname': 'area', 'mv': 0},
                   'source': 'JRC E.1 Space, Security, Migration',
                   'reference': 'JRC E.1 Space, Security, Migration',
                   'geographical': {'datum': ''}
                   }
    array_to_nc_from_clone(masknc_path, ldd_map_path, catchments_mask, metadata=nc_metadata)

    # In order to keep the same functionality as before we need to generate also the smaller mask map 
    x_min, x_max, y_min, y_max = bbox_from_netcdf(masknc_path)
    cutmap(masknc_path, smallmask_path, x_min, x_max, y_min, y_max, use_coords = True)

    return smallmask_path, outlets_path, masknc_path
