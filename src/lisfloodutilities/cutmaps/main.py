"""
Copyright 2019-2026 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

cutmaps
A tool to cut netcdf files
"""

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Optional, List, Union

from .. import version, logger
from .cutlib import (mask_from_ldd, get_filelist, get_cuts, cutmap,
                     MASK_VALUE, SMALL_MASK_FILENAME, FULL_MASK_FILENAME, OUTLETS_FILENAME)
from .helpers import COORDINATE_NAMES, TIME_NAMES
from netCDF4 import Dataset 
import xarray as xr
import numpy as np

# Variables that should be excluded from mask processing.
# Common NetCDF projection/CRS variable names
EXCLUDED_VAR_NAMES = {
    "crs",
    "wgs_1984",
    "lambert_azimuthal_equal_area",
    "lambert_conformal_conic",
    "spatial_ref",
    "rotated_pole",
    "albers_conical_equal_area",
    "polar_stereographic",
    "transverse_mercator",
    "oblique_mercator",
    "lambert_cylindrical_equal_area",
    "latitude_longitude",
    "geostationary",
    "crs_wkt",
    "projection",
    "grid_mapping",
    "esri_pe_string",
}

# Add coordinate and time coordinate names to the exclusion set
EXCLUDED_VAR_NAMES.update(COORDINATE_NAMES)
EXCLUDED_VAR_NAMES.update(TIME_NAMES)


def parse_and_check_args(parser, cliargs):
    args = parser.parse_args(cliargs)
    if (args.mask!=None) + (args.cuts!=None) + (args.cuts_indices!=None) > 1:
        parser.error('[--mask | --cuts | --cuts_indices] arguments are mutually exclusive')
    if not (args.mask or args.cuts or args.cuts_indices) and not (args.ldd and args.stations):
        parser.error('(--mask | --cuts | --cuts_indices | [--ldd, --stations]) You need to pass mask path or cuts coordinates '
                     'or a list of stations along with LDD path')
    if (args.mask or args.cuts or args.cuts_indices) and (args.ldd or args.stations):
        parser.error('(--mask | --cuts | --cuts_indices | [--ldd, --stations]) '
                     '--mask, --cuts, --cuts_indices and --ldd and --stations arguments are mutually exclusive')
    return args

def get_arg_coords(value):
    """Parse coordinate values from command line argument.
    
    Args:
        value: String containing 4 space-separated values (coordinates or indices)
        
    Returns:
        Map object yielding the parsed values as int or float
        
    Raises:
        argparse.ArgumentTypeError: If the value doesn't contain exactly 4 values
    """
    apply = float if '.' in value else int  # user can provide coords (float) or matrix indices bbox (int)
    values = value.split()
    if len(values) != 4:
        raise argparse.ArgumentTypeError(
            f"Expected 4 values, got {len(values)}: '{value}'"
        )
    values = map(apply, values)
    return values

def load_mask_file(mask_path: str) -> Optional[np.ndarray]:
    """Load mask values from a NetCDF mask file.

    Args:
        mask_path: Path to the mask NetCDF file.

    Returns:
        NumPy array with mask values, or None if loading fails.
    """
    try:
        with xr.open_dataset(mask_path, engine='netcdf4') as mask_ds:
            mask_var = next(
                (mask_ds[var] for var in mask_ds.data_vars if var not in COORDINATE_NAMES),
                None,
            )
            if mask_var is not None:
                return mask_var.values  # Load into memory once.
    except Exception as exc:  # pragma: no cover – defensive programming.
        logger.exception("Failed to read mask file %s: %s", mask_path, exc)
    return None


def _apply_variable_mask(
    variable,
    mask_values: np.ndarray,
) -> None:
    """Apply a mask to a NetCDF variable, replacing masked values with fill_value.

    This function handles both 2D (raster) and 3D (time-series of rasters) data.
    It preserves the variable's scale_factor and add_offset attributes when
    computing the encoded fill value.

    Args:
        variable: NetCDF4 Variable object to mask.
        mask_values: 2D NumPy array with mask values (1 = keep, other = fill).
    """
    # Skip string variables and known coordinate/reference variables.
    if variable.dtype == np.dtype('|S1') or variable._name in EXCLUDED_VAR_NAMES:
        return

    # Retrieve scaling metadata (if any).
    scale_factor = getattr(variable, 'scale_factor', None)
    add_offset = getattr(variable, 'add_offset', None)
    fill_value = getattr(variable, '_FillValue', np.nan)

    # Apply scale/offset to the fill value so that NaNs are encoded correctly.
    if (
        fill_value is not None
        and fill_value == fill_value  # NaN-safe check
        and scale_factor is not None
        and add_offset is not None
    ):
        fill_value = fill_value * scale_factor + add_offset

    # Load data once and apply mask in a single vectorised operation.
    data = variable[:]

    if data.ndim == 2:
        # 2‑D case (e.g., raster).
        masked = np.where(mask_values == MASK_VALUE, data, fill_value)
        variable[:] = masked
    elif data.ndim == 3:
        # 3‑D case (e.g., time‑series of rasters).
        # Reshape mask to broadcast across the time axis, then apply in one call.
        mask_3d = mask_values[np.newaxis, :, :]
        masked = np.where(mask_3d == MASK_VALUE, data, fill_value)
        variable[:] = masked
    else:
        logger.warning(
            "Skipping variable '%s' with unsupported dimensionality %d",
            variable._name,
            data.ndim,
        )

class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):
        group_mask = self.add_argument_group(title='Cut with a provided mask or a bounding box or '
                                                   'create mask cookie-cutter on-fly from stations list and ldd map')
        group_filelist = self.add_mutually_exclusive_group(required=True)
        group_mask.add_argument("-m", "--mask", help='mask file cookie-cutter in pcraster (.map) or netcdf (.nc) format')
        group_mask.add_argument("-c", "--cuts", help='Cut coordinates in the form "lonmin lonmax latmin latmax" using coordinates bounding box', type=get_arg_coords)
        group_mask.add_argument("-i", "--cuts_indices", help='Cut coordinates in the form "imin imax jmin jmax" using matrix indices', type=get_arg_coords)
        group_mask.add_argument("-l", "--ldd", help='Path to LDD file in netcdf format (.nc)')
        group_mask.add_argument("-N", "--stations",
                                help='Path to stations.txt file.'
                                     'Read documentation to know about the format')

        group_filelist.add_argument("-f", "--folder", help='Directory with netCDF files to be cut')
        group_filelist.add_argument("-F", "--file", help='netCDF file to be cut')
        group_filelist.add_argument("-S", "--subdir",
                                    help='Directory containing folders '
                                         'Output files will have same directory-folders structure')

        self.add_argument("-o", "--outpath", help='path where to save cut files',
                          default='./cutmaps_out', required=True)
        self.add_argument("-W", "--overwrite", help='Set flag to overwrite existing files',
                          default=False, required=False, action='store_true')


def main(cliargs):
    parser = ParserHelpOnError(description='Cut netCDF file: {}'.format(version))
    parser.add_args()
    args = parse_and_check_args(parser, cliargs)
    mask = args.mask
    cuts = args.cuts
    cuts_indices = args.cuts_indices
    mask_nc = None

    ldd = args.ldd
    stations = args.stations

    # =========================================================================
    # Input Path Handling with Validation and Error Handling
    # =========================================================================
    # Extract input paths from arguments with proper validation
    input_folder = args.folder
    input_file = args.file
    static_data_folder = args.subdir
    overwrite = args.overwrite
    pathout = args.outpath

    # Validate input sources - ensure at least one input is provided
    # This provides early failure with clear error messages
    if not any([input_folder, input_file, static_data_folder]):
        parser.error(
            'No input source specified. Please provide one of: '
            '--folder, --file, or --subdir'
        )

    # Validate input_folder if provided - check existence and accessibility
    if input_folder:
        input_folder_path = Path(input_folder)
        if not input_folder_path.exists():
            raise FileNotFoundError(
                f"Input folder does not exist: {input_folder}"
            )
        if not input_folder_path.is_dir():
            raise NotADirectoryError(
                f"Input path is not a directory: {input_folder}"
            )
        # Resolve to absolute path for consistent handling
        input_folder = str(input_folder_path.resolve())

    # Validate input_file if provided
    if input_file:
        input_file_path = Path(input_file)
        if not input_file_path.exists():
            raise FileNotFoundError(
                f"Input file does not exist: {input_file}"
            )
        if not input_file_path.is_file():
            raise ValueError(
                f"Input path is not a file: {input_file}"
            )
        # Resolve to absolute path for consistent handling
        input_file = str(input_file_path.resolve())

    # Validate static_data_folder if provided
    if static_data_folder:
        static_data_path = Path(static_data_folder)
        if not static_data_path.exists():
            raise FileNotFoundError(
                f"Static data folder does not exist: {static_data_folder}"
            )
        if not static_data_path.is_dir():
            raise NotADirectoryError(
                f"Static data path is not a directory: {static_data_folder}"
            )
        # Resolve to absolute path for consistent handling
        static_data_folder = str(static_data_path.resolve())

    # Validate output path
    pathout_path = Path(pathout)
    if pathout_path.exists() and not pathout_path.is_dir():
        raise ValueError(
            f"Output path exists but is not a directory: {pathout}"
        )
    if not os.path.exists(pathout):
        logger.warning('\nOutput folder %s not existing. Creating it...', pathout)
        os.mkdir(pathout)
    if ldd and stations:
        logger.info('\nTry to produce a mask from LDD and stations points: %s %s', ldd, stations)
        mask, outlets_nc, mask_nc = mask_from_ldd(ldd, stations)
        # copy outlets (produced from stations txt file) and the new mask to output folder
        shutil.copy(outlets_nc, os.path.join(pathout, OUTLETS_FILENAME))
        shutil.copy(mask, os.path.join(pathout, SMALL_MASK_FILENAME))
        shutil.copy(mask_nc, os.path.join(pathout, FULL_MASK_FILENAME))

    x_min, x_max, y_min, y_max = get_cuts(cuts=cuts, cuts_indices=cuts_indices, mask=mask)
    logger.info('\n\nCutting using: %s\n Files to cut from: %s\n Output: %s\n Overwrite existing: %s\n\n',
                mask or ([x_min, x_max, y_min, y_max if cuts or cuts_indices else None]),
                input_folder or static_data_folder,
                pathout, overwrite)

    list_to_cut = get_filelist(input_folder, static_data_folder, input_file)

    # walk through list_to_cut
    for file_to_cut in list_to_cut:

        filename, ext = os.path.splitext(file_to_cut)

        # localdir used only with static_data_folder.
        # It will track folder structures in a EFAS/GloFAS like setup and replicate it in output folder
        localdir = os.path.dirname(file_to_cut)\
            .replace(os.path.dirname(static_data_folder), '')\
            .lstrip('/') if static_data_folder else ''

        fileout = os.path.join(pathout, localdir, os.path.basename(file_to_cut))
        if os.path.isdir(file_to_cut) and static_data_folder:
            # just create folder
            os.makedirs(fileout, exist_ok=True)
            continue
        if ext != '.nc':
            if static_data_folder:
                logger.warning('%s is not in netcdf format, just copying to ouput folder', file_to_cut)
                shutil.copy(file_to_cut, fileout)
            else:
                logger.warning('%s is not in netcdf format, skipping...', file_to_cut)
            continue
        elif os.path.isfile(fileout) and os.path.exists(fileout) and not overwrite:
            logger.warning('%s already existing. This file will not be overwritten', fileout)
            continue

        logger.info(f"## Processing file (cutmap): {file_to_cut} -> {fileout}")
        cutmap(file_to_cut, fileout, x_min, x_max, y_min, y_max, use_coords=(cuts_indices is None))
        logger.info(f"## Finished processing file (cutmap): {file_to_cut} -> {fileout}")
        if ldd and stations:
            logger.info(f"## Applying mask to file: {fileout}")
            mask_path = os.path.join(pathout, SMALL_MASK_FILENAME)
            mask_map_values = load_mask_file(mask_path)

            if mask_map_values is None:
                logger.error(f"Failed to load mask file: {mask_path}")
                continue

            with Dataset(fileout, 'r+', format='NETCDF4_CLASSIC') as file_out:
                for output_var_name, variable in file_out.variables.items():
                    _apply_variable_mask(
                        variable=variable,
                        mask_values=mask_map_values,
                    )
            logger.info(f"## Finished applying mask to file: {fileout}")

def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
