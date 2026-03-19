import os
import argparse
import sys
import time
from pathlib import Path
from typing import List, Tuple, Union, Optional
import numpy as np
from netCDF4 import Dataset
import xarray as xr
from pyproj import CRS
from earthkit.hydro import river_network, data_structures


LATITUDE_VARIABLES = ['y', 'lat', 'latitude', 'nlat', 'lats', 'latitudes']
LONGITUDE_VARIABLES = ['x', 'lon', 'longitude', 'nlon', 'lons', 'longitudes']
# Used to verify if the name of a variable is a lat variable
LATITUDE_NAMES = set(LATITUDE_VARIABLES)
# Used to verify if the name of a variable is a lon variable
LONGITUDE_NAMES = set(LONGITUDE_VARIABLES)
# Used to verify if the name of a variable is a lat or lon variable
COORDINATE_NAMES = LATITUDE_NAMES | LONGITUDE_NAMES
# Used to verify if the name of a variable is a time variable
TIME_NAMES = {'time', 't'}
# Used to find the corresponding lon var name to a given lat var name
LATITUDE_NAME_PAIR = dict(zip(LATITUDE_VARIABLES, LONGITUDE_VARIABLES))
# Used to find the corresponding lat var name to a given lon var name
LONGITUDE_NAME_PAIR = dict(zip(LONGITUDE_VARIABLES, LATITUDE_VARIABLES))


def verify_existing_netcdf(input_file: Optional[Union[Path, str]] = None, file_id: str = '') -> str:
    """
    Check if input_file is an existing netCDF.
    ----------
    Parameters:
    
    input_file: Path to the file
    file_id: String identifying the file for the error message
    ----------
    Returns:

    Error message if is invalid
    Empty string if is valid
    """
    msg = ''
    if input_file is None:
        msg = f'Error: The {file_id} must be provided.'
    else:
        input_file_path = input_file if isinstance(input_file, Path) else Path(input_file)
        if (not input_file_path.exists() or
            not input_file_path.is_file() or
            input_file_path.suffix != '.nc'):
            msg = f'Error: The {file_id} must exist and be a netCDF file.'
    return msg


def get_from_metadata(metadata: dict[str, Union[str, int, float, bool, dict[str, Union[str, int, float, bool]]]],
                      main_key: str, sub_key: str,
                      default_value: Union[str, int, float, bool]) -> Union[str, int, float, bool]:
    """
    Retrieve a nested value from a metadata dictionary.

    Safely accesses a two-level nested dictionary, returning a default
    value if either the main key or sub-key is missing.

    Parameters
    ----------
    metadata : dict
        The dictionary containing metadata with nested structure.
    main_key : str
        The top-level key to look up in the metadata.
    sub_key : str
        The second-level key to look up within the main_key's value.
    default_value : Any
        The value to return if either key is not found.

    Returns
    -------
    Any
        The value associated with sub_key if found, otherwise default_value.
    """
    # Validate inputs to fail fast with clear error messages
    if metadata is None:
        raise TypeError("metadata cannot be None")
    if not isinstance(metadata, dict):
        raise TypeError(f"metadata must be a dict, got {type(metadata).__name__}")
    if not isinstance(main_key, str):
        raise TypeError(f"main_key must be a str, got {type(main_key).__name__}")
    if not isinstance(sub_key, str):
        raise TypeError(f"sub_key must be a str, got {type(sub_key).__name__}")

    main_value = metadata.get(main_key)

    # Return default if main_key doesn't exist
    if main_value is None:
        return default_value

    if not isinstance(main_value, dict):
        return main_value

    return main_value.get(sub_key, default_value)


def get_river_network_from_map(ldd_map: Union[Path, str],
                               precomputed_network: Optional[Union[Path, str]] = None,
                               export: bool = True)-> data_structures.RiverNetwork:
    """
    Get the river network from a LDD map. If a precomputed network file exists it is used,
    otherwise it is created and saved in the ldd_map folder for future use.
    Parameters
    ----------
    ldd_map: Path to the LDD map file.
    precomputed_network: Optional path to a precomputed network file.
    export: If True, the network will be exported to a file.
    Returns
    -------
    data_structures.RiverNetwork
        The river network object.
    """
    ldd_map_path = Path(ldd_map) if isinstance(ldd_map, str) else ldd_map
    if precomputed_network is not None:
        precomputed_network_path = Path(precomputed_network) if isinstance(precomputed_network, str) else precomputed_network
    else:
        precomputed_network_path = Path(ldd_map_path.parent, f'{ldd_map_path.stem}.joblib')
    # If the precomputed network file exists, load and return it
    if precomputed_network_path.exists():
        network = river_network.create(precomputed_network_path.as_posix(), "precomputed", "file")
        return network
    # Creates the new network and saves it for future use
    network = river_network.create(ldd_map_path.as_posix(), "pcr_d8", "file")
    if export:
        network.export(precomputed_network_path.as_posix())
    return network


def read_column_file(
    file_path: Path,
    delimiter: str = ' ',
    skip_header: bool = False,
) -> np.ndarray[tuple[int, int, int], np.dtype[np.float64]]:
    """
    Load a three‑column numeric file and return an (N, 3) array:
        column 0: x coordinate
        column 1: y coordinate
        column 2: value
    Parameters
    ----------
    file_path: Path to the column file.
    delimiter: Delimiter to use (default: any whitespace).
    skip_header: If True, the first line is ignored (useful for a header row).
    Returns
    -------
    np.ndarray
        Shape (N, 3), dtype = float64 (values are cast later as needed).

    x, y, value
    values in column file that have x,y coordinates at the upper and left margins
    of a cell come into that cell, values at the bottom and right margins come
    into neighbouring cells. So, cell values with x, y coordinates at vertexes
    of cells come into the cell at the lower right side of the vertex.
    """
    data = None
    try:
        with open(file_path) as f:
            data = np.loadtxt(
                (r.replace('\t', delimiter) for r in f),
                delimiter=delimiter,
                skiprows=1 if skip_header else 0,
            )
    except Exception as exc:
        raise RuntimeError(f"Failed to read column file {file_path}: {exc}")

    if data.ndim == 1:
        # Only one line, make it a 2‑D row vector
        data = data[np.newaxis, :]

    if data.shape[1] < 3:
        raise ValueError(
            f"Column file {file_path} must contain at least three columns (x, y, value)."
        )
    # Keep only the first three columns
    return data[:, :3].astype(np.float64)


def copy_clone_geometry(src_ds: Dataset, dst_ds: Dataset):
    """
    Replicate dimensions, coordinate variables and global attributes
    from the clone (src) to the new dataset (dst).
    """
    # ----------- dimensions -------------------------------------------------
    for dim_name, dim in src_ds.dimensions.items():
        dst_ds.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))
    total_dimensions = len(src_ds.dimensions.items())
    # ----------- coordinate variables ---------------------------------------
    # We copy everything that is not the main data variable (i.e. any
    # variable whose dimensions are a subset of the dimensions we just created).
    for var_name, var in src_ds.variables.items():
        # Skip data variables that have more than one dimension (e.g. the raster itself).
        # Most clones only have coordinate variables (e.g. lat, lon, x, y, time).
        if len(var.dimensions) < total_dimensions:
            # Create the variable with the same datatype and attributes.
            new_var = dst_ds.createVariable(
                var_name,
                var.datatype,
                var.dimensions,
                zlib=True,  # compress by default
            )
            # Copy attributes (including units, standard_name, etc.)
            new_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            # Copy the data itself.
            new_var[:] = var[:]
    # ----------- global attributes ------------------------------------------
    dst_ds.setncatts({k: src_ds.getncattr(k) for k in src_ds.ncattrs()})


def get_crs(ds: Dataset) -> CRS:
    """
    Try to obtain the coordinate system/projection from a dataset.
    ------------
    Returns a CRS object containing the information 
    """
    crs_var = None
    # If it exists, identify the variable containing the grid mapping
    for v in ds.variables.values():
        if (getattr(v, "grid_mapping_name", None) is not None or
            getattr(v, "grid_mapping", None) is not None):
            crs_var = v
            break

    proj_wkt = None
    proj4   = None
    # Identify in the variable properties the string containing
    # the projection definitions
    if crs_var is not None:
        proj_wkt = (getattr(crs_var, "spatial_ref", None) or
                    getattr(crs_var, "crs_wkt", None) or
                    getattr(crs_var, "esri_pe_string", None)
                    )
        proj4   = getattr(crs_var, "proj4_params", None)
    # If Identify the string containing the projection definitions
    # in the dataset itself.
    if not proj_wkt and not proj4:
        proj_wkt = (getattr(ds, "spatial_ref", None) or
                    getattr(ds, "esri_pe_string", None)
                    )
        proj4   = getattr(ds, "proj4_params", None)
    # Build the CRS object from the projection definitions
    if proj_wkt:
        crs = CRS.from_wkt(proj_wkt)
    elif proj4:
        crs = CRS.from_proj4(proj4)
    else:
        # As a last resort assume generic WGS84 geographic
        crs = CRS.from_epsg(4326)
    return crs


def copy_datum(source: xr.Dataset, target: xr.Dataset) -> None:
    """
    Copy geographical datum attributes from source Dataset to target Dataset.
    
    Copies the following attributes from the source to the target:
    - esri_pe_string
    - spatial_ref
    - crs_wkt
    - grid_mapping
    
    Parameters:
    -----------
    source : xr.Dataset
        The source Dataset from which to copy datum attributes.
    target : xr.Dataset
        The target Dataset to which datum attributes will be copied.
    
    Returns:
    --------
    xr.Dataset
        The target Dataset with copied datum attributes.
    """
    # List of datum attributes to copy
    datum_attrs = ['esri_pe_string', 'spatial_ref', 'crs_wkt']
    
    # Get the main variable from the source Dataset (first data variable)
    if len(source.data_vars) > 0:
        main_var_name = list(source.data_vars)[0]
        main_var = source[main_var_name]
        
        # Copy attributes from the main variable
        for attr_name in datum_attrs:
            if attr_name in main_var.attrs:
                target.attrs[attr_name] = main_var.attrs[attr_name]
    
    # Also check for grid_mapping coordinate and copy if present
    if 'grid_mapping' in source.coords:
        target = target.assign_coords({'grid_mapping': source.coords['grid_mapping']})
        if 'grid_mapping' not in target.attrs and 'grid_mapping' in source.data_vars:
            target.attrs['grid_mapping'] = source['grid_mapping'].values
        elif 'grid_mapping' not in target.attrs:
            # Try to get from the main variable's attributes
            if len(source.data_vars) > 0:
                main_var_name = list(source.data_vars)[0]
                if 'grid_mapping' in source[main_var_name].attrs:
                    target.attrs['grid_mapping'] = source[main_var_name].attrs['grid_mapping']


def read_spatial_dimensions(ds: Union[Dataset, xr.Dataset]) -> Tuple[str, str]:
    """
    Identify the two spatial dimensions (usually y, x)
    Detect the dimension variables where the order is important: (y, x) or (x, y)
    -------------
    Returns (dim_x, dim_y)
    """
    
    if isinstance(ds, Dataset):
        spatial_dims = [d for d in ds.dimensions if not ds.dimensions[d].isunlimited()]
    else:
        spatial_dims = [d for d in ds.dims]
    if len(spatial_dims) < 2:
        raise ValueError(
            f"The data set does not contain at least two spatial dimensions {spatial_dims}."
        )
    # Transpose the dimension names to get x, y order
    dim_y, dim_x = [str(c) for c in spatial_dims[:2]]
    if not dim_y.lower() in LATITUDE_NAMES:
        tmp = dim_y
        dim_y = dim_x
        dim_x = tmp
    return dim_x, dim_y


def bbox_from_netcdf(path: Path, time_index: int = 0) -> Tuple[float, float, float, float]:
    """
    Get the minimum bounding box possible that includes all data from a netCDF grid.
    If there is more than 1 raster, selects the raster at time_index position.
    --------------
    Returns:
    min_x, max_x, min_y, max_y
    """
    ds = xr.open_dataset(path)
    var_names = [n for n,da in ds.data_vars.items() if len(da.dims)>1]
    var_name = var_names[0]
    da = ds[var_name]

    var_time_names = [n for n,_ in ds.data_vars.items() if n in TIME_NAMES]
    time_dim = None if len(var_time_names) == 0 else var_time_names[0]
    if time_dim and time_dim in da.dims:
        da = da.isel({time_dim: time_index})

    x_dim, y_dim = read_spatial_dimensions(ds)

    # guarantee (y, x) order
    da = da.transpose(..., y_dim, x_dim)

    # Boolean mask of valid data
    y_idx, x_idx = np.where(da.notnull().values)

    if y_idx.size == 0:
        raise ValueError("All values are NaN / masked.")

    # min and max indexes
    i_min, i_max = x_idx.min(), x_idx.max()
    j_min, j_max = y_idx.min(), y_idx.max()

    # coordinate values corresponding to the min and max indexes
    coord_idx_min_x = float(ds[x_dim].values[i_min])
    coord_idx_max_x = float(ds[x_dim].values[i_max])
    coord_idx_min_y = float(ds[y_dim].values[j_min])
    coord_idx_max_y = float(ds[y_dim].values[j_max])
    
    # min and max coordinates
    min_x = min(coord_idx_min_x, coord_idx_max_x)
    max_x = max(coord_idx_min_x, coord_idx_max_x)
    min_y = min(coord_idx_min_y, coord_idx_max_y)
    max_y = max(coord_idx_min_y, coord_idx_max_y)
    
    ds.close()

    return min_x, max_x, min_y, max_y


def get_fill_value_packed(ds: Dataset,
                          default_fill_value: Optional[Union[np.int32, int, float] | None] = None
                        ) -> Tuple[Union[np.int32, int, float], Union[np.int32, int, float]]:
    """
    Identifies the main variable in the given netCDF4 dataset and returns its fill value, 
    accounting for scale_factor and add_offset if present.
    default_fill_value can be provided to override the fill value from the dataset; if not provided, it defaults to -9999.
    
    Parameters:
    dataset (netCDF4.Dataset): The netCDF4 dataset object.

    Returns:
    A tuple containing the original fill value and the packed fill value after applying scale and offset.
    """
    main_variable = None
    max_dimensions = 0

    # Iterate through the variables to identify the main variable
    for var_name, var in ds.variables.items():
        # Assuming the main variable has the most dimensions
        if len(var.dimensions) > max_dimensions:
            main_variable = var
            max_dimensions = len(var.dimensions)

    if main_variable is None:
        raise ValueError("No variables found in the dataset.")

    # Retrieve scaling metadata (if any)
    scale_factor = getattr(main_variable, 'scale_factor', None)
    add_offset = getattr(main_variable, 'add_offset', None)
    fill_value = getattr(main_variable, '_FillValue', -9999) if default_fill_value is None else default_fill_value
    fill_value_packed = fill_value
    # Apply scale/offset to the fill value so that NaNs are encoded correctly.
    if fill_value is not None and scale_factor is not None and add_offset is not None:
        fill_value_packed = fill_value * scale_factor + add_offset

    return fill_value, fill_value_packed


def array_to_nc_from_clone(out_path: Path, clone_path: Path, grid: np.ndarray,
                           metadata: dict[str, Union[str, int, float, bool, dict[str, Union[str, int, float, bool]]]] = {}):
    """
    Create a NetCDF raster (out_path) that mirrors clone_path getting the data from grid which
    should have the same shape as clone.

    Parameters
    -------------
    out_path: Output netCDF file path
    clone_path: Clone netCDF file path
    grid: array of data to be written in the out_path file
    metadata: metadata for the out_file
    """
    default_var_name = 'map'
    # Open clone and obtain geometry and coordinate arrays
    with Dataset(clone_path, "r") as src:
        dim_x, dim_y = read_spatial_dimensions(src)

        # Prepare the raster (filled with NetCDF fill value)
        fill_value = np.iinfo(np.int32).min
        fill_value = np.int32(get_from_metadata(metadata, 'variable', 'mv', fill_value))
        fill_value, fill_value_packed = get_fill_value_packed(src, fill_value)
        raster = grid.copy()
        raster[np.isnan(raster)] = fill_value_packed

        # Write the new NetCDF (copy geometry and raster variable)
        with Dataset(out_path, "w", format="NETCDF4") as dst:
            copy_clone_geometry(src, dst)
            dst.history = 'Created {}'.format(time.ctime(time.time()))
            dst.Conventions = 'CF-1.7'
            dst.Source_Software = 'JRC.E1 lisfloodutilities'
            dst.source = metadata.get('source', '')
            dst.reference = metadata.get('reference', '')

            compression_kwargs = {"zlib": True, "complevel": 4}
            nc_var_name = str(get_from_metadata(metadata, 'variable', 'shortname', default_var_name))
            var = dst.createVariable(
                nc_var_name,
                np.int32,
                (dim_y, dim_x),
                fill_value=fill_value_packed,
                **compression_kwargs,
            )
            var.long_name = str(get_from_metadata(metadata, 'variable', 'longname', ''))
            var.standard_name = nc_var_name
            var.units = str(get_from_metadata(metadata, 'variable', 'units', ''))
            var[:, :] = raster


def write_output_nc(out_path: Path, clone_path: Path, points: np.ndarray[tuple[int, int, int], np.dtype[np.float64]],
    metadata: dict[str, Union[str, int, float, bool, dict[str, Union[str, int, float, bool]]]] = {},
    nodata_val: int = 0, var_name: str = "map",
    compress: bool = True, quiet: bool = True) -> np.ndarray[tuple[int, int], np.dtype[np.int32]]:
    """
    Create a NetCDF raster (out_path) that mirrors clone_path and fills it with
    the values supplied in points (x, y, value).

    Parameters
    -------------
    out_path: Path,
    clone_path: Path,
    points: 3D array (N, 3) -> x, y, value,
    metadata: metadata for the output_file
    nodata_val: int = 0,
    var_name: str = "map",
    compress: compress the output file with zlib and compression level 4
    quiet: Do not print output for the user. (default: True)

    Returns
    -------------
    For convenience returns the result also as a np.ndarray 2D array (N, 2) with the raster values.
    """
    # Open clone, obtain geometry and coordinate arrays
    with Dataset(clone_path, "r") as src:
        dim_x, dim_y = read_spatial_dimensions(src)

        ny = len(src.dimensions[dim_y])
        nx = len(src.dimensions[dim_x])

        # Pull the coordinate vectors (1‑D arrays)
        # We assume the clone has coordinate variables named exactly as the
        # dimensions (common convention).
        try:
            coord_x = src.variables[dim_x][:]      # e.g. easting or longitude
            coord_y = src.variables[dim_y][:]      # e.g. northing or latitude
        except KeyError as exc:
            raise KeyError(
                f"Clone file {clone_path} does not contain coordinate variables "
                f"named after its dimensions ({dim_x}, {dim_y})."
            ) from exc

        diff_x = np.diff(coord_x)
        diff_y = np.diff(coord_y)
        # Ensure they are 1‑D and monotonic (required for searchsorted)
        if coord_x.ndim != 1 or coord_y.ndim != 1:
            raise ValueError("Coordinate variables must be 1‑D arrays.")
        if not (np.all(diff_x > 0) or np.all(diff_x < 0)):
            raise ValueError("X coordinate array must be strictly monotonic.")
        if not (np.all(diff_y > 0) or np.all(diff_y < 0)):
            raise ValueError("Y coordinate array must be strictly monotonic.")

        # Determine the order of the y axis (ascending or descending)
        # This is important for correctly mapping the y coordinates to row indices.
        # Initialize defaults in case diff_y is neither strictly increasing nor decreasing
        sorted_y_side = "left"
        sorted_y_offset = 0
        if np.all(diff_y >= 0):
            # Ascending order (e.g. latitude increasing from south to north)
            sorted_y_side = "left"
            sorted_y_offset = 0
        elif np.all(diff_y <= 0):
            # Descending order (e.g. northing decreasing from top to bottom)
            sorted_y_side = "right"
            sorted_y_offset = 1

        resolution_x = np.round(np.abs(coord_x[1] - coord_x[0]), decimals=2)
        resolution_y = np.round(np.abs(coord_y[1] - coord_y[0]), decimals=2)

        # Save the order of the coordinates to be able to find the correct index in the raster
        coord_x_backup_indexes = {}
        for i, x in enumerate(coord_x):
            coord_x_backup_indexes[x] = i
        coord_y_backup_indexes = {}
        for i, y in enumerate(coord_y):
            coord_y_backup_indexes[y] = i

        # Prepare an empty raster (filled with NetCDF fill value)
        fill_value = np.iinfo(np.int32).min if nodata_val is None else nodata_val
        fill_value = np.int32(get_from_metadata(metadata, 'variable', 'mv', fill_value))
        fill_value, fill_value_packed = get_fill_value_packed(src, fill_value)
        raster = np.full((ny, nx), fill_value_packed, dtype=np.int32)

        # From the stations file map each (x, y) -> (col, row) index
        xs = points[:, 0]
        ys = points[:, 1]
        vals = points[:, 2]

        # Find the index of the nearest coordinate in the clone.
        # Because the coordinate arrays already sorted we can use np.searchsorted.
        # We also verify that the coordinate matches exactly (within a tiny tolerance)
        # otherwise we raise a warning; you can change the tolerance as needed.
        # Tolerance settings:
        #   * ``tol``  – absolute tolerance, set to a fraction (10%) of the
        #                larger grid resolution.  This scales the allowed
        #                deviation with the spatial resolution of the clone.
        #   * ``rol_x`` and ``rol_y`` – relative tolerance factors.  The original
        #                implementation used ``0.5 * resolution``; we keep the same
        #                behaviour but expose the values as explicit variables for
        #                clarity and easy tweaking.
        #   * These tolerances are used with ``np.isclose`` to verify that the
        #                input coordinates match the clone grid within an acceptable
        #                margin.
        tol = 0.1 * max(resolution_x, resolution_y)
        rol_x = 1.0  # relative tolerance for X (originally 0.5 * resolution_x)
        rol_y = 1.0  # relative tolerance for Y (originally 0.5 * resolution_y)

        coord_x_sorted = np.sort(coord_x)

        # X values
        col_idx = np.searchsorted(coord_x_sorted, xs, side="left")
        
        # Adjust indices that fall on the right edge
        col_idx = np.clip(col_idx, 0, nx - 1)
        
        # Verify matches
        x_match = np.isclose(coord_x_sorted[col_idx], xs, atol=tol, rtol=rol_x)

        if not np.all(x_match):
            bad = np.where(~x_match)[0]
            raise ValueError(
                f"The following {len(bad)} x‑coordinates do not match any column of the clone:\n"
                + "\n".join(f"  point {i}: x={xs[i]}" for i in bad)
            )

        coord_y_sorted = np.sort(coord_y)

        # Y values
        # NOTE: In many raster conventions Y increases downwards (row 0 = top).
        # NetCDF often stores Y increasing upwards (south to north). The mapping
        # below works for both; we simply locate the index and later use it as the
        # row number.
        # Determine the row index for each y‑coordinate.
        # ``np.searchsorted`` with ``side="right"`` returns the insertion point
        # after any equal values. Subtracting one gives the index of the matching
        # coordinate (or the nearest lower coordinate). This avoids the off‑by‑
        # one shift that occurs when using ``side="left"`` on ascending sorted
        # coordinates.
        row_idx = np.searchsorted(coord_y_sorted, ys, side=sorted_y_side) - sorted_y_offset

        # Clip indices to the valid range [0, ny‑1] for robustness.
        row_idx = np.clip(row_idx, 0, ny - 1)

        # Verify that the located rows are within the tolerance.
        y_match = np.isclose(coord_y_sorted[row_idx], ys, atol=tol, rtol=rol_y)

        if not np.all(y_match):
            bad = np.where(~y_match)[0]
            raise ValueError(
                f"The following {len(bad)} y‑coordinates do not match any row of the clone:\n"
                + "\n".join(f"  point {i}: y={ys[i]}" for i in bad)
            )

        # Insert the values (vectorised)
        # Cast values to int32 (large‑map). Values equal to the nodata
        # are left as the fill value.
        vals_int = vals.astype(np.int32)
        if nodata_val is not None:
            vals_int[vals_int == nodata_val] = fill_value_packed

        # Use advanced indexing to write all points at once
        col_idx_saved = np.array(list(map(lambda idx: coord_x_backup_indexes[coord_x_sorted[idx]], col_idx)))
        row_idx_saved = np.array(list(map(lambda idx: coord_y_backup_indexes[coord_y_sorted[idx]], row_idx)))
        raster[row_idx_saved, col_idx_saved] = vals_int

        # Write the new NetCDF (copy geometry and raster variable)
        with Dataset(out_path, "w", format="NETCDF4") as dst:
            copy_clone_geometry(src, dst)
            dst.history = 'Created {}'.format(time.ctime(time.time()))
            dst.Conventions = 'CF-1.7'
            dst.Source_Software = 'JRC.E1 lisfloodutilities'
            dst.source = metadata.get('source', '')
            dst.reference = metadata.get('reference', '')

            compression_kwargs = {"zlib": True, "complevel": 4} if compress else {}
            nc_var_name = str(get_from_metadata(metadata, 'variable', 'shortname', var_name))
            var = dst.createVariable(
                nc_var_name,
                np.int32,
                (dim_y, dim_x),
                fill_value=fill_value_packed,
                **compression_kwargs,
            )
            var.long_name = str(get_from_metadata(metadata, 'variable', 'longname', ''))
            var.standard_name = nc_var_name
            var.units = str(get_from_metadata(metadata, 'variable', 'units', ''))
            var[:, :] = raster

    if not quiet:
        print(
            f"Created {out_path} ({out_path.stat().st_size/1e6:.2f} MB) "
            f"with {np.count_nonzero(raster != fill_value_packed)} populated cells."
        )
    return raster


def col2netcdf(column_file: Union[Path, str], output_file: Union[Path, str], clone_file: Union[Path, str],
               metadata: dict[str, Union[str, int, float, bool, dict[str, Union[str, int, float, bool]]]] = {},
               var_name: str = 'map', nodata: int = 0, delimiter: str = ' ', skip_header: bool = False,
               quiet: bool = True) -> np.ndarray[tuple[int, int], np.dtype[np.int32]]:
    """
    Re‑implementation of the PCRaster command using netCDF4 and numpy.
    
    col2map F0 F1 -N --clone F2
    
    Objective: Convert a column text file to a NetCDF raster using the geometry of a clone NetCDF file.
    
    Parameters
    ----------
    column_file: Path to the input column file (F0)
    output_file: Path for the output NetCDF raster (F1)
    clone_file: Path to the clone NetCDF file (F2) that defines the grid
    metadata: metadata for the output_file
    var_name: Name of the raster variable inside the NetCDF file (default: 'map')
    nodata: Value in the column file that should be treated as missing data (default: 0)
    delimiter: Delimiter for the column file (default: whitespace)
    skip_header: Skip the first line of the column file (useful for a header) (default: False)
    quiet: Do not print output for the user. (default: True)

    Returns
    -------------
    For convenience returns the result also as a np.ndarray2D array (N, 2) with the raster values.
    """
    column_file_path = column_file if isinstance(column_file, Path) else Path(column_file)
    output_file_path = output_file if isinstance(output_file, Path) else Path(output_file)
    clone_file_path = clone_file if isinstance(clone_file, Path) else Path(clone_file)
    # Load the column values
    column_values = read_column_file(column_file_path, delimiter=delimiter, skip_header=skip_header)
    # Write the output NetCDF
    output_array = write_output_nc(out_path=output_file_path, clone_path=clone_file_path,
                                   points=column_values, metadata=metadata, nodata_val=nodata,
                                   var_name=var_name, quiet=quiet)
    return output_array
