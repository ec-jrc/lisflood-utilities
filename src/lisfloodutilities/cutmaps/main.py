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
import logging
import multiprocessing
import os
import shutil
import sys
import time
import warnings
from pathlib import Path
from typing import Optional, List, Union, NoReturn

import xarray as xr
import numpy as np

import dask
from dask.diagnostics.progress import ProgressBar

from .. import version, logger
from .cutlib import (
    mask_from_ldd,
    get_filelist,
    get_cuts,
    cutmap,
    MASK_VALUE,
    SMALL_MASK_FILENAME,
    FULL_MASK_FILENAME,
    OUTLETS_FILENAME,
)
from .helpers import COORDINATE_NAMES, TIME_NAMES

# Maximum number of retry attempts for a frozen file before giving up
MAX_RETRIES = 3
# Timeout in seconds: if output file doesn't change for this duration, consider processing frozen
FREEZE_TIMEOUT_SECONDS = 300

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


def _is_nfs_path(path: Union[str, Path]) -> bool:
    """Check if a path resides on an NFS filesystem.

    Uses /proc/mounts (Linux) to determine the filesystem type of the mount
    point that contains the given path. Handles autofs-managed NFS mounts
    by also checking whether the mount source references /etc/auto.nfs or
    the actual underlying mount is NFS.

    Args:
        path: File or directory path to check.

    Returns:
        True if the path is on an NFS (nfs, nfs4, or autofs-backed NFS)
        filesystem, False otherwise.
    """
    try:
        resolved = Path(path).resolve()
    except (OSError, ValueError):
        return False

    try:
        with open("/proc/mounts", "r") as f:
            mounts = f.readlines()
    except (OSError, IOError):
        # /proc/mounts not available (non-Linux or restricted environment)
        return False

    # Find the longest matching mount point for the resolved path.
    # Track both the best match and any NFS mount that matches
    # (autofs may mount the NFS volume with a longer path than the autofs entry).
    best_match = ""
    best_fstype = ""
    best_source = ""
    for line in mounts:
        parts = line.split()
        if len(parts) < 3:
            continue
        source = parts[0]
        mount_point = parts[1]
        fstype = parts[2]
        # Check if the resolved path starts with this mount point
        if (str(resolved) == mount_point or
                str(resolved).startswith(mount_point.rstrip("/") + "/")):
            if len(mount_point) >= len(best_match):
                best_match = mount_point
                best_fstype = fstype
                best_source = source

    fstype_lower = best_fstype.lower()

    # Direct NFS mount
    if fstype_lower in ("nfs", "nfs4"):
        return True

    # autofs mount backed by NFS (source is /etc/auto.nfs or similar auto.nfs* file)
    if fstype_lower == "autofs" and "auto.nfs" in best_source.lower():
        return True

    return False


def _configure_nfs_workarounds(*paths: Optional[str]) -> None:
    """Enable HDF5/dask workarounds if any of the given paths reside on NFS.

    On NFS filesystems, HDF5 file locking can cause intermittent hangs because
    the NFS lock daemon may be slow or unresponsive. Additionally, dask's default
    threaded scheduler can amplify the problem by opening the same file from many
    threads simultaneously.

    This function:
    - Sets HDF5_USE_FILE_LOCKING=FALSE to disable advisory locks.
    - Switches dask to the synchronous scheduler to avoid concurrent file access.

    These are only applied when at least one path is on NFS and the user has not
    already set HDF5_USE_FILE_LOCKING explicitly.

    Args:
        *paths: Paths to check (None values are skipped).
    """
    any_nfs = any(_is_nfs_path(p) for p in paths if p is not None)
    if not any_nfs:
        return

    if "HDF5_USE_FILE_LOCKING" not in os.environ:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        logger.info(
            "NFS filesystem detected. Setting HDF5_USE_FILE_LOCKING=FALSE "
            "to prevent potential hangs."
        )

    # Use synchronous scheduler to avoid multi-threaded NFS lock contention.
    dask.config.set(scheduler="synchronous")
    logger.info(
        "NFS filesystem detected. Using dask synchronous scheduler "
        "to avoid concurrent file access."
    )


def parse_and_check_args(parser: argparse.ArgumentParser, cliargs: List[str]) -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Args:
        parser: ArgumentParser instance with configured arguments.
        cliargs: List of command-line arguments to parse.

    Returns:
        Parsed arguments namespace.

    Raises:
        SystemExit: If argument validation fails.
    """
    args = parser.parse_args(cliargs)

    # Check mutual exclusivity of mask-related arguments
    mask_args_count = sum([
        args.mask is not None,
        args.cuts is not None,
        args.cuts_indices is not None,
    ])
    if mask_args_count > 1:
        parser.error("[--mask | --cuts | --cuts_indices] arguments are mutually exclusive")

    # Check that at least one mask/cut method is provided
    has_mask_method = args.mask or args.cuts or args.cuts_indices
    has_ldd_stations = args.ldd and args.stations
    if not has_mask_method and not has_ldd_stations:
        parser.error(
            "(--mask | --cuts | --cuts_indices | [--ldd, --stations]) "
            "You need to pass mask path or cuts coordinates "
            "or a list of stations along with LDD path"
        )

    # Check mutual exclusivity between mask methods and ldd/stations
    if has_mask_method and has_ldd_stations:
        parser.error(
            "(--mask | --cuts | --cuts_indices | [--ldd, --stations]) "
            "--mask, --cuts, --cuts_indices and --ldd and --stations arguments are mutually exclusive"
        )

    return args


def get_arg_coords(value: str) -> map:
    """Parse coordinate values from command line argument.

    Args:
        value: String containing 4 space-separated values (coordinates or indices).

    Returns:
        Map object yielding the parsed values as int or float.

    Raises:
        argparse.ArgumentTypeError: If the value doesn't contain exactly 4 values.
    """
    apply_type = float if "." in value else int
    values = value.split()
    if len(values) != 4:
        raise argparse.ArgumentTypeError(
            f"Expected 4 values, got {len(values)}: '{value}'"
        )
    return map(apply_type, values)


def load_mask_file(mask_path: Union[str, Path]) -> Optional[np.ndarray]:
    """Load mask values from a NetCDF mask file.

    Args:
        mask_path: Path to the mask NetCDF file.

    Returns:
        NumPy array with mask values, or None if loading fails.
    """
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Duplicate dimension names present",
                category=UserWarning,
            )
            with xr.open_dataset(mask_path, engine="netcdf4") as mask_ds:
                mask_var = next(
                    (mask_ds[var] for var in mask_ds.data_vars if var not in COORDINATE_NAMES),
                    None,
                )
                if mask_var is not None:
                    return mask_var.values  # Load into memory once.
    except Exception as exc:  # pragma: no cover - defensive programming.
        logger.exception("Failed to read mask file %s: %s", mask_path, exc)
    return None


def _apply_variable_mask_xarray(
    dataset: xr.Dataset,
    var_name: str,
    mask_values: np.ndarray,
) -> None:
    """Apply a mask to an xarray DataArray variable, replacing masked values with fill_value.

    This function handles both 2D (raster) and 3D (time-series of rasters) data.
    It preserves the variable's scale_factor and add_offset attributes when
    computing the encoded fill value.

    Args:
        dataset: xarray Dataset containing the variable to mask.
        var_name: Name of the variable to mask.
        mask_values: 2D NumPy array with mask values (1 = keep, other = fill).
    """
    var = dataset[var_name]

    # Skip string variables and known coordinate/reference variables.
    if var.dtype == np.dtype("|S1") or var_name in EXCLUDED_VAR_NAMES:
        return

    # Skip variables with duplicate dimension names (e.g. string/label variables
    # encoded with repeated character dimensions like ('string1', 'string1')).
    # These are not spatial rasters and cannot be masked.
    if len(var.dims) != len(set(var.dims)):
        logger.debug(
            "Skipping variable '%s' with duplicate dimensions %s",
            var_name,
            var.dims,
        )
        return

    # Retrieve scaling metadata (if any).
    scale_factor = var.attrs.get("scale_factor")
    add_offset = var.attrs.get("add_offset")
    fill_value = var.attrs.get("_FillValue")

    # Handle missing_value as fallback for _FillValue
    if fill_value is None:
        fill_value = var.attrs.get("missing_value")

    # Default fill value handling
    if fill_value is None:
        fill_value = np.nan

    # Apply scale/offset to the fill value so that NaNs are encoded correctly.
    if (
        fill_value is not None
        and fill_value == fill_value  # NaN-safe check
        and scale_factor is not None
        and add_offset is not None
    ):
        fill_value = fill_value * scale_factor + add_offset

    # Load data once and apply mask in a single vectorised operation.
    data = var.values

    if data.ndim == 2:
        # 2-D case (e.g., raster).
        masked = np.where(mask_values == MASK_VALUE, data, fill_value)
        dataset[var_name].values = masked
    elif data.ndim == 3:
        # 3-D case (e.g., time-series of rasters).
        # Reshape mask to broadcast across the time axis, then apply in one call.
        mask_3d = mask_values[np.newaxis, :, :]
        masked = np.where(mask_3d == MASK_VALUE, data, fill_value)
        dataset[var_name].values = masked
    else:
        logger.warning(
            "Skipping variable '%s' with unsupported dimensionality %d",
            var_name,
            data.ndim,
        )


def _process_single_file(
    file_to_cut: str,
    pathout: str,
    static_data_folder: Optional[str],
    x_min, x_max, y_min, y_max: float,
    use_coords: bool,
    overwrite: bool,
    ldd: Optional[str],
    stations: Optional[str],
) -> Optional[str]:
    """Process a single file for cutting and masking.
    
    Args:
        file_to_cut: Path to the file to process.
        pathout: Output directory path.
        static_data_folder: Static data folder path (for directory structure).
        x_min, x_max, y_min, y_max: Cut coordinates.
        use_coords: Whether to use coordinates (True) or indices (False).
        overwrite: Whether to overwrite existing files.
        ldd: Path to LDD file (if using mask from LDD/stations).
        stations: Path to stations file (if using mask from LDD/stations).
    
    Returns:
        Path to the processed file, or None if skipped/failed.
    """
    filename, ext = os.path.splitext(file_to_cut)
    
    # localdir used only with static_data_folder.
    localdir = (
        os.path.dirname(file_to_cut)
        .replace(os.path.dirname(static_data_folder), "")
        .lstrip("/")
        if static_data_folder
        else ""
    )
    
    fileout = os.path.join(pathout, localdir, os.path.basename(file_to_cut))
    
    if os.path.isdir(file_to_cut) and static_data_folder:
        os.makedirs(fileout, exist_ok=True)
        return None
    
    if ext != ".nc":
        if static_data_folder:
            logger.warning(
                "%s is not in netcdf format, just copying to output folder",
                file_to_cut,
            )
            shutil.copy(file_to_cut, fileout)
        else:
            logger.warning("%s is not in netcdf format, skipping...", file_to_cut)
        return None
    
    if os.path.isfile(fileout) and not overwrite:
        logger.warning("%s already existing. This file will not be overwritten", fileout)
        return None
    
    logger.info(f"## Processing file (cutmap): {file_to_cut} -> {fileout}")
    cutmap(file_to_cut, fileout, x_min, x_max, y_min, y_max, use_coords=use_coords)
    logger.info(f"## Finished processing file (cutmap): {file_to_cut} -> {fileout}")
    
    if ldd and stations:
        logger.info(f"## Applying mask to file: {fileout}")
        mask_path = os.path.join(pathout, SMALL_MASK_FILENAME)
        mask_map_values = load_mask_file(mask_path)
        
        if mask_map_values is None:
            logger.error(f"Failed to load mask file: {mask_path}")
            return fileout  # Return anyway, cut was successful
        
        # Load dataset, apply mask, then write back to a temp file and replace
        # Suppress xarray warning about duplicate dimension names in string/label
        # variables — we skip those variables during mask application.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Duplicate dimension names present",
                category=UserWarning,
            )
            with xr.open_dataset(fileout, engine="netcdf4", decode_cf=False) as ds:
                for var_name in ds.data_vars:
                    _apply_variable_mask_xarray(
                        dataset=ds,
                        var_name=str(var_name),
                        mask_values=mask_map_values,
                    )
                # Write the masked dataset to a temporary file
                tmp_file = fileout + ".tmp"
                ds.to_netcdf(tmp_file)
        # Replace original with the masked version
        os.replace(tmp_file, fileout)
        logger.info(f"## Finished applying mask to file: {fileout}")
    
    return fileout


def _run_process_single_file_worker(
    result_queue: multiprocessing.Queue,
    file_to_cut: str,
    pathout: str,
    static_data_folder: Optional[str],
    x_min, x_max, y_min, y_max,
    use_coords: bool,
    overwrite: bool,
    ldd: Optional[str],
    stations: Optional[str],
) -> None:
    """Worker function that runs _process_single_file in a child process.

    Communicates the result (or exception) back through a multiprocessing Queue.
    """
    try:
        result = _process_single_file(
            file_to_cut, pathout, static_data_folder,
            x_min, x_max, y_min, y_max,
            use_coords, overwrite, ldd, stations,
        )
        result_queue.put(("ok", result))
    except Exception as exc:
        result_queue.put(("error", str(exc)))


def _get_output_file_path(
    file_to_cut: str,
    pathout: str,
    static_data_folder: Optional[str],
) -> str:
    """Compute the expected output file path for a given input file.

    This mirrors the logic in _process_single_file for determining fileout.
    """
    localdir = (
        os.path.dirname(file_to_cut)
        .replace(os.path.dirname(static_data_folder), "")
        .lstrip("/")
        if static_data_folder
        else ""
    )
    return os.path.join(pathout, localdir, os.path.basename(file_to_cut))


def _process_with_watchdog(
    file_to_cut: str,
    pathout: str,
    static_data_folder: Optional[str],
    x_min, x_max, y_min, y_max,
    use_coords: bool,
    overwrite: bool,
    ldd: Optional[str],
    stations: Optional[str],
    freeze_timeout: int = FREEZE_TIMEOUT_SECONDS,
    max_retries: int = MAX_RETRIES,
) -> Optional[str]:
    """Process a single file with watchdog monitoring for frozen operations.

    Runs _process_single_file in a child process and monitors progress by checking
    the output file's modification time. If the file doesn't change for longer than
    freeze_timeout seconds, the process is killed, the output file is deleted, and
    processing is retried. After max_retries failed attempts, raises a RuntimeError.

    Args:
        file_to_cut: Path to the file to process.
        pathout: Output directory path.
        static_data_folder: Static data folder path.
        x_min, x_max, y_min, y_max: Cut coordinates or indices.
        use_coords: Whether to use coordinates or indices.
        overwrite: Whether to overwrite existing files.
        ldd: Path to LDD file.
        stations: Path to stations file.
        freeze_timeout: Seconds without progress before considering the process frozen.
        max_retries: Maximum number of retry attempts.

    Returns:
        Path to the processed file, or None if skipped.

    Raises:
        RuntimeError: If processing fails after max_retries attempts due to freezing.
    """
    fileout = _get_output_file_path(file_to_cut, pathout, static_data_folder)

    for attempt in range(1, max_retries + 1):
        logger.info(
            "Processing file (attempt %d/%d): %s",
            attempt, max_retries, file_to_cut,
        )

        result_queue = multiprocessing.Queue()
        proc = multiprocessing.Process(
            target=_run_process_single_file_worker,
            args=(
                result_queue, file_to_cut, pathout, static_data_folder,
                x_min, x_max, y_min, y_max,
                use_coords, overwrite, ldd, stations,
            ),
        )
        proc.start()

        # Monitor the process for freezes
        last_activity_time = time.time()
        last_file_mtime = None
        frozen = False

        while proc.is_alive():
            # Check if the output file has been modified
            try:
                if os.path.exists(fileout):
                    current_mtime = os.path.getmtime(fileout)
                    if last_file_mtime is None or current_mtime != last_file_mtime:
                        last_file_mtime = current_mtime
                        last_activity_time = time.time()
                # Also check for .tmp files (used during mask application)
                tmp_file = fileout + ".tmp"
                if os.path.exists(tmp_file):
                    tmp_mtime = os.path.getmtime(tmp_file)
                    if tmp_mtime != last_file_mtime:
                        last_activity_time = time.time()
            except OSError:
                pass  # File may be in flux

            # Check if we've exceeded the freeze timeout
            elapsed_since_activity = time.time() - last_activity_time
            if elapsed_since_activity > freeze_timeout:
                frozen = True
                logger.warning(
                    "WATCHDOG: Processing of %s appears frozen "
                    "(no progress for %d seconds). Killing process (attempt %d/%d).",
                    file_to_cut, freeze_timeout, attempt, max_retries,
                )
                proc.kill()
                proc.join(timeout=10)
                if proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=5)
                break

            # Check result queue for early completion
            if not result_queue.empty():
                break

            time.sleep(2)  # Poll every 2 seconds

        # If process finished naturally
        if not frozen:
            proc.join(timeout=10)
            if not result_queue.empty():
                status, result = result_queue.get_nowait()
                if status == "ok":
                    return result
                else:
                    logger.error(
                        "Error processing %s: %s", file_to_cut, result
                    )
                    # Non-freeze error, don't retry with watchdog logic
                    raise RuntimeError(
                        f"Error processing {file_to_cut}: {result}"
                    )
            # Process ended without putting result in queue (crashed?)
            logger.warning(
                "Process for %s ended unexpectedly (exit code: %s)",
                file_to_cut, proc.exitcode,
            )
            frozen = True  # Treat as a failure worth retrying

        # Clean up the output file after a frozen/failed attempt
        for f_path in (fileout, fileout + ".tmp"):
            if os.path.exists(f_path):
                try:
                    os.remove(f_path)
                    logger.info("Deleted incomplete output file: %s", f_path)
                except OSError as e:
                    logger.warning("Could not delete %s: %s", f_path, e)

        if attempt < max_retries:
            logger.info(
                "Retrying file %s (next attempt: %d/%d)...",
                file_to_cut, attempt + 1, max_retries,
            )

    # All retries exhausted
    raise RuntimeError(
        f"FATAL: Processing of file '{file_to_cut}' froze {max_retries} times "
        f"(no progress for {freeze_timeout}s each time). Giving up. "
        f"This file may be corrupted or too large to process."
    )


class ParserHelpOnError(argparse.ArgumentParser):
    """ArgumentParser that prints help on error."""

    def error(self, message: str) -> NoReturn:
        """Print error message and help, then exit."""
        sys.stderr.write(f"Error: {message}\n")
        self.print_help()
        sys.exit(1)

    def add_args(self) -> None:
        """Add command-line arguments to the parser."""
        group_mask = self.add_argument_group(
            title="Cut with a provided mask or a bounding box or "
            "create mask cookie-cutter on-fly from stations list and ldd map"
        )
        group_filelist = self.add_mutually_exclusive_group(required=True)

        group_mask.add_argument(
            "-m", "--mask",
            help="mask file cookie-cutter in pcraster (.map) or netcdf (.nc) format"
        )
        group_mask.add_argument(
            "-c", "--cuts",
            help='Cut coordinates in the form "lonmin lonmax latmin latmax" using coordinates bounding box',
            type=get_arg_coords
        )
        group_mask.add_argument(
            "-i", "--cuts_indices",
            help='Cut coordinates in the form "imin imax jmin jmax" using matrix indices',
            type=get_arg_coords
        )
        group_mask.add_argument(
            "-l", "--ldd",
            help="Path to LDD file in netcdf format (.nc)"
        )
        group_mask.add_argument(
            "-N", "--stations",
            help="Path to stations.txt file. Read documentation to know about the format"
        )

        group_filelist.add_argument(
            "-f", "--folder",
            help="Directory with netCDF files to be cut"
        )
        group_filelist.add_argument(
            "-F", "--file",
            help="netCDF file to be cut"
        )
        group_filelist.add_argument(
            "-S", "--subdir",
            help="Directory containing folders. "
            "Output files will have same directory-folders structure"
        )

        self.add_argument(
            "-o", "--outpath",
            help="path where to save cut files",
            default="./cutmaps_out",
            required=True
        )
        self.add_argument(
            "-W", "--overwrite",
            help="Set flag to overwrite existing files",
            default=False,
            required=False,
            action="store_true"
        )
        self.add_argument(
            "-v", "--verbose",
            help="Show detailed info messages during processing",
            default=False,
            required=False,
            action="store_true"
        )
        self.add_argument(
            "-t", "--freeze-timeout",
            help="Enable freeze detection watchdog. If a file does not advance "
                 f"for more than {FREEZE_TIMEOUT_SECONDS} seconds, it is killed "
                 f"and retried up to {MAX_RETRIES} times. "
                 "Without this flag, processing runs without timeout monitoring.",
            default=False,
            required=False,
            action="store_true"
        )


def _validate_input_path(
    path: Optional[str],
    path_type: str,
    must_exist: bool = True,
    must_be_dir: bool = False,
    must_be_file: bool = False,
) -> Optional[str]:
    """Validate an input path and return the resolved absolute path.

    Args:
        path: Path to validate.
        path_type: Description of the path for error messages.
        must_exist: Whether the path must exist.
        must_be_dir: Whether the path must be a directory.
        must_be_file: Whether the path must be a file.

    Returns:
        Resolved absolute path as string, or None if path is not provided.

    Raises:
        FileNotFoundError: If path must exist but doesn't.
        NotADirectoryError: If path must be a directory but isn't.
        ValueError: If path must be a file but isn't.
    """
    if path is None:
        return None

    path_obj = Path(path)

    if must_exist and not path_obj.exists():
        raise FileNotFoundError(f"{path_type} does not exist: {path}")

    if must_be_dir and not path_obj.is_dir():
        raise NotADirectoryError(f"{path_type} is not a directory: {path}")

    if must_be_file and not path_obj.is_file():
        raise ValueError(f"{path_type} is not a file: {path}")

    return str(path_obj.resolve())


def main(cliargs: List[str]) -> None:
    """Main entry point for the cutmaps tool.

    Args:
        cliargs: Command-line arguments (excluding script name).
    """
    parser = ParserHelpOnError(description=f"Cut netCDF file: {version}")
    parser.add_args()
    args = parse_and_check_args(parser, cliargs)

    # Set logger verbosity based on --verbose flag
    if not args.verbose:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    # Extract arguments
    mask = args.mask
    cuts = args.cuts
    cuts_indices = args.cuts_indices
    ldd = args.ldd
    stations = args.stations
    input_folder = args.folder
    input_file = args.file
    static_data_folder = args.subdir
    overwrite = args.overwrite
    pathout = args.outpath

    # Validate input sources
    if not any([input_folder, input_file, static_data_folder]):
        parser.error(
            "No input source specified. Please provide one of: "
            "--folder, --file, or --subdir"
        )

    # Validate input paths
    input_folder = _validate_input_path(input_folder, "Input folder", must_be_dir=True)
    input_file = _validate_input_path(input_file, "Input file", must_be_file=True)
    static_data_folder = _validate_input_path(
        static_data_folder, "Static data folder", must_be_dir=True
    )

    # Validate and create output path
    pathout_path = Path(pathout)
    if pathout_path.exists() and not pathout_path.is_dir():
        raise ValueError(f"Output path exists but is not a directory: {pathout}")

    if not os.path.exists(pathout):
        logger.warning("\nOutput folder %s not existing. Creating it...", pathout)
        os.mkdir(pathout)

    # Apply NFS workarounds (disable HDF5 locking, use synchronous dask scheduler)
    # if any of the input/output paths reside on an NFS filesystem.
    _configure_nfs_workarounds(input_folder, input_file, static_data_folder, pathout)

    mask_nc = None
    if ldd and stations:
        logger.info("\nTry to produce a mask from LDD and stations points: %s %s", ldd, stations)
        mask, outlets_nc, mask_nc = mask_from_ldd(ldd, stations)
        # Copy outlets (produced from stations txt file) and the new mask to output folder
        shutil.copy(outlets_nc, os.path.join(pathout, OUTLETS_FILENAME))
        shutil.copy(mask, os.path.join(pathout, SMALL_MASK_FILENAME))
        shutil.copy(mask_nc, os.path.join(pathout, FULL_MASK_FILENAME))

    x_min, x_max, y_min, y_max = get_cuts(cuts=cuts, cuts_indices=cuts_indices, mask=mask)

    logger.info(
        "\n\nCutting using: %s\n Files to cut from: %s\n Output: %s\n Overwrite existing: %s\n\n",
        mask or ([x_min, x_max, y_min, y_max if cuts or cuts_indices else None]),
        input_folder or static_data_folder,
        pathout,
        overwrite,
    )

    list_to_cut = get_filelist(
        input_folder or "", static_data_folder or "", input_file or ""
    )

    # Process files sequentially
    num_files = len(list_to_cut)
    use_coords = cuts_indices is None
    freeze_timeout_enabled = args.freeze_timeout

    if freeze_timeout_enabled:
        logger.info(f"Processing {num_files} files with watchdog monitoring "
                    f"(freeze timeout: {FREEZE_TIMEOUT_SECONDS}s, max retries: {MAX_RETRIES})")
    else:
        logger.info(f"Processing {num_files} files sequentially")

    for file_to_cut in list_to_cut:
        try:
            if freeze_timeout_enabled:
                result = _process_with_watchdog(
                    str(file_to_cut),
                    pathout,
                    static_data_folder,
                    x_min, x_max, y_min, y_max,
                    use_coords,
                    overwrite,
                    ldd,
                    stations,
                )
            else:
                result = _process_single_file(
                    str(file_to_cut),
                    pathout,
                    static_data_folder,
                    x_min, x_max, y_min, y_max,
                    use_coords,
                    overwrite,
                    ldd,
                    stations,
                )
            if result:
                logger.info(f"Successfully processed: {result}")
        except RuntimeError as exc:
            logger.error(str(exc))
            sys.exit(1)
        except Exception as exc:
            logger.error(f"File {file_to_cut} generated an exception: {exc}")
            sys.exit(1)


def main_script() -> None:
    """Entry point for script execution."""
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
