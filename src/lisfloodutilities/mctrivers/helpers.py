import xarray as xr
from typing import Tuple, Optional
import numpy as np

# Define common coordinate name patterns for automatic detection
X_COORD_NAMES = ('lon', 'x', 'rlon')
Y_COORD_NAMES = ('lat', 'y', 'rlat')


def copy_coordinates_and_attributes(source: xr.Dataset, target: xr.DataArray, y_proj: str, x_proj: str) -> xr.DataArray:
    """
    Copy coordinates and attributes from source Dataset to target DataArray.
    
    Parameters:
    -----------
    source : xr.Dataset
        The Dataset from which to copy attributes.
    target : xr.DataArray
        The DataArray to which attributes will be copied.
    
    Returns:
    --------
    xr.DataArray
        The target DataArray with copied attributes.
    """
    # target.attrs = source.attrs.copy()
    target = target.assign_coords({y_proj: source.coords[y_proj], x_proj: source.coords[x_proj]})
    
    # Get the main variable from the source Dataset (first data variable)
    if len(source.data_vars) > 0:
        main_var_name = list(source.data_vars)[0]
        main_var = source[main_var_name]
        
        # Copy attributes from the main variable (including grid_mapping, esri_pe_string, etc.)
        for attr_name in ['esri_pe_string', 'spatial_ref', 'crs_wkt', 'grid_mapping']:
            if attr_name in main_var.attrs:
                target.attrs[attr_name] = main_var.attrs[attr_name]
    
    # Also check for grid_mapping coordinate and copy if present
    if 'grid_mapping' in source.coords:
        target = target.assign_coords({'grid_mapping': source.coords['grid_mapping']})
        if 'grid_mapping' not in target.attrs:
            target.attrs['grid_mapping'] = source.coords['grid_mapping'].values
    return target
