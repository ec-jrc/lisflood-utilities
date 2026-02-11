import os
import argparse
import sys
import time
from pathlib import Path
from typing import List, Tuple
import numpy as np
from netCDF4 import Dataset
from pyproj import CRS


def pcraster_command(cmd, files=None):

    if files is not None:
        # replace placeholders in command with actual filenames
        for alias, realname in files.items():
            cmd = cmd.replace(alias, '"{}"'.format(realname))

    os.system(cmd)
    return cmd

def get_from_metadata(metadata: dict, main_key: str, sub_key: str, default_value):
    return metadata[main_key].get(sub_key, default_value) if main_key in metadata else default_value

def read_column_file(
    file_path: Path,
    delimiter: str = None,
    skip_header: bool = False,
) -> np.ndarray:
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
        # Only one line → make it a 2‑D row vector
        data = data[np.newaxis, :]

    if data.shape[1] < 3:
        raise ValueError(
            f"Column file {file_path} must contain at least three columns (x, y, value)."
        )
    # Keep only the first three columns
    return data[:, :3]


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
        # if set(var.dimensions).issubset(set(src_ds.dimensions)):
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
    crs_var = None
    for v in ds.variables.values():
        if (getattr(v, "grid_mapping_name", None) is not None or
            getattr(v, "grid_mapping", None) is not None):
            crs_var = v
            break

    proj_wkt = None
    proj4   = None
    if crs_var is not None:
        proj_wkt = (getattr(crs_var, "spatial_ref", None) or
                    getattr(crs_var, "crs_wkt", None) or
                    getattr(crs_var, "esri_pe_string", None)
                    )
        proj4   = getattr(crs_var, "proj4_params", None)

    if not proj_wkt and not proj4:
        proj_wkt = (getattr(ds, "spatial_ref", None) or
                    getattr(crs_var, "esri_pe_string", None)
                    )
        proj4   = getattr(ds, "proj4_params", None)

    if proj_wkt:
        crs = CRS.from_wkt(proj_wkt)
    elif proj4:
        crs = CRS.from_proj4(proj4)
    else:
        # As a last resort assume generic WGS84 geographic
        crs = CRS.from_epsg(4326)
    return crs

def read_spatial_dimensions(ds: Dataset) -> Tuple[str, str]:
    """
    Identify the two spatial dimensions (usually y, x)
    Detect the dimension variables where the order is important: (y, x) or (x, y)
    -------------
    Returns (dim_x, dim_y)
    """
    spatial_dims = [d for d in ds.dimensions if not ds.dimensions[d].isunlimited()]
    if len(spatial_dims) < 2:
        raise ValueError(
            f"The clone file {clone_path} does not contain at least two spatial dimensions."
        )
    dim_y, dim_x = spatial_dims[:2]
    if not dim_y.lower() in ['y', 'lat', 'latitude', 'nlat']:
        dim_x, dim_y = spatial_dims[:2]
    return dim_x, dim_y

def array_to_nc_from_clone(out_path: Path, clone_path: Path, grid: np.ndarray, metadata: dict = {}):
    """
    Create a NetCDF raster (out_path) that mirrors clone_path getting the data from grid which
    should have the same shape as clone.

    Parameters
    -------------
    out_path: Path,
    clone_path: Path,
    grid: array of data
    metadata: metadata for the output_file
    """
    default_var_name = 'map'
    # Open clone – obtain geometry and coordinate arrays
    with Dataset(clone_path, "r") as src:
        dim_x, dim_y = read_spatial_dimensions(src)

        # Prepare the raster (filled with NetCDF fill value)
        fill_value = np.iinfo(np.int32).min
        fill_value = np.int32(get_from_metadata(metadata, 'variable', 'mv', fill_value))
        raster = grid
        raster[raster is np.nan] = fill_value

        # Write the new NetCDF (copy geometry + raster variable)
        with Dataset(out_path, "w", format="NETCDF4") as dst:
            copy_clone_geometry(src, dst)
            dst.history = 'Created {}'.format(time.ctime(time.time()))
            dst.Conventions = 'CF-1.7'
            dst.Source_Software = 'JRC.E1 lisfloodutilities'
            dst.source = metadata.get('source', '')
            dst.reference = metadata.get('reference', '')

            compression_kwargs = {"zlib": True, "complevel": 4}
            nc_var_name = get_from_metadata(metadata, 'variable', 'shortname', default_var_name)
            var = dst.createVariable(
                nc_var_name,
                np.int32,
                (dim_y, dim_x),
                fill_value=fill_value,
                **compression_kwargs,
            )
            var.long_name = get_from_metadata(metadata, 'variable', 'longname', '')
            var.standard_name = nc_var_name
            var.units = get_from_metadata(metadata, 'variable', 'units', '')
            var[:, :] = raster

def find_one(keys: List[str], obj):
    """
    Return the first value found for `key`, or None if not present.
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k.lower() in keys:
                return v
            # Dive deeper only if we haven't found it yet
            result = find_one(keys, v)
            if result is not None:
                return result
    elif isinstance(obj, list):
        for item in obj:
            result = find_one(keys, item)
            if result is not None:
                return result
    return None


def write_output_nc(out_path: Path, clone_path: Path, points: np.ndarray,
    metadata: dict = {}, nodata_val: np.int32 = 0, var_name: str = "map",
    compress: bool = True, quiet: bool = True) -> np.ndarray:
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
    For convenience returns the result also as a np.ndarray
    """
    # Open clone – obtain geometry and coordinate arrays
    with Dataset(clone_path, "r") as src:
        dim_x, dim_y = read_spatial_dimensions(src)

        ny = len(src.dimensions[dim_y])
        nx = len(src.dimensions[dim_x])

        # Pull the coordinate vectors (1‑D arrays)
        # We assume the clone has *coordinate variables* named exactly as the
        # dimensions (common convention).
        try:
            coord_x = src.variables[dim_x][:]      # e.g. easting or longitude
            coord_y = src.variables[dim_y][:]      # e.g. northing or latitude
        except KeyError as exc:
            raise KeyError(
                f"Clone file {clone_path} does not contain coordinate variables "
                f"named after its dimensions ({dim_x}, {dim_y})."
            ) from exc

        # Ensure they are 1‑D and monotonic (required for searchsorted)
        if coord_x.ndim != 1 or coord_y.ndim != 1:
            raise ValueError("Coordinate variables must be 1‑D arrays.")
        if not (np.all(np.diff(coord_x) > 0) or np.all(np.diff(coord_x) < 0)):
            raise ValueError("X coordinate array must be strictly monotonic.")
        if not (np.all(np.diff(coord_y) > 0) or np.all(np.diff(coord_y) < 0)):
            raise ValueError("Y coordinate array must be strictly monotonic.")

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
        raster = np.full((ny, nx), fill_value, dtype=np.int32)

        # From the stations file map each (x, y) -> (row, col) index
        xs = points[:, 0]
        ys = points[:, 1]
        vals = points[:, 2]

        # print('minX', xs.min(), 'maxX', xs.max())
        # print('minY', ys.min(), 'maxY', ys.max())
        #
        # crs_obj = get_crs(src)
        # if crs_obj.is_projected:
        #     crs_dict = crs_obj.to_dict()
        #     false_easting = int(crs_dict.get('x_0'))
        #     false_northing = int(crs_dict.get('y_0'))
        #     # xs = xs + false_easting
        #     ys = ys + false_northing
        # print(f'false_easting: {false_easting} false_northing: {false_northing}')
        # print('PROJ minX', xs.min(), 'maxX', xs.max())
        # print('PROJ minY', ys.min(), 'maxY', ys.max())
        # print('PROJ mincoord_x', coord_x.min(), 'maxcoord_x', coord_x.max())
        # print('PROJ mincoord_y', coord_y.min(), 'maxcoord_y', coord_y.max())

        # Find the index of the nearest coordinate in the clone.
        # Because the coordinate arrays already sorted we can use np.searchsorted.
        # We also verify that the coordinate matches exactly (within a tiny tolerance)
        # otherwise we raise a warning; you can change the tolerance as needed.
        tol = 1e-6

        coord_x_sorted = np.sort(coord_x)
        # X values
        col_idx = np.searchsorted(coord_x_sorted, xs, side="left")
        
        # Adjust indices that fall on the right edge
        col_idx[col_idx == nx] = nx - 1
        
        # Verify matches
        x_match = np.isclose(coord_x_sorted[col_idx], xs, atol=tol)
        if not np.all(x_match):
            bad = np.where(~x_match)[0]
            raise ValueError(
                f"The following {len(bad)} x‑coordinates do not match any column of the clone:\n"
                + "\n".join(f"  row {i}: x={xs[i]}" for i in bad)
            )

        coord_y_sorted = np.sort(coord_y)
        # Y values
        # NOTE: In many raster conventions Y increases downwards (row 0 = top).
        # NetCDF often stores Y increasing upwards (south to north). The mapping
        # below works for both; we simply locate the index and later use it as the
        # row number.
        row_idx = np.searchsorted(coord_y_sorted, ys, side="left")

        row_idx[row_idx == ny] = ny - 1

        y_match = np.isclose(coord_y_sorted[row_idx], ys, atol=tol)
        if not np.all(y_match):
            bad = np.where(~y_match)[0]
            raise ValueError(
                f"The following {len(bad)} y‑coordinates do not match any row of the clone:\n"
                + "\n".join(f"  row {i}: y={ys[i]}" for i in bad)
            )

        # Insert the values (vectorised)
        # Cast values to int32 (large‑map). Values equal to the nodata
        # are left as the fill value.
        vals_int = vals.astype(np.int32)
        if nodata_val is not None:
            vals_int[vals_int == nodata_val] = fill_value

        # Use advanced indexing to write all points at once
        col_idx_saved = np.array(list(map(lambda idx: coord_x_backup_indexes[coord_x_sorted[idx]], col_idx)))
        row_idx_saved = np.array(list(map(lambda idx: coord_y_backup_indexes[coord_y_sorted[idx]], row_idx)))
        raster[row_idx_saved, col_idx_saved] = vals_int

        # Write the new NetCDF (copy geometry + raster variable)
        with Dataset(out_path, "w", format="NETCDF4") as dst:
            copy_clone_geometry(src, dst)
            dst.history = 'Created {}'.format(time.ctime(time.time()))
            dst.Conventions = 'CF-1.7'
            dst.Source_Software = 'JRC.E1 lisfloodutilities'
            dst.source = metadata.get('source', '')
            dst.reference = metadata.get('reference', '')

            compression_kwargs = {"zlib": True, "complevel": 4} if compress else {}
            nc_var_name = get_from_metadata(metadata, 'variable', 'shortname', var_name)
            var = dst.createVariable(
                nc_var_name,
                np.int32,
                (dim_y, dim_x),
                fill_value=fill_value,
                **compression_kwargs,
            )
            var.long_name = get_from_metadata(metadata, 'variable', 'longname', '')
            var.standard_name = nc_var_name
            var.units = get_from_metadata(metadata, 'variable', 'units', '')
            var[:, :] = raster

    if not quiet:
        print(
            f"Created {out_path} ({out_path.stat().st_size/1e6:.2f} MB) "
            f"with {np.count_nonzero(raster != fill_value)} populated cells."
        )
    return raster


def col2netcdf(column_file: Path, output_file: Path, clone_file: Path, metadata: dict = {}, var_name: str = 'map',
               nodata: int = 0, delimiter: str = ' ', skip_header: bool = False, quiet: bool = True) -> np.ndarray:
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
    For convenience returns the result also as a np.ndarray
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
