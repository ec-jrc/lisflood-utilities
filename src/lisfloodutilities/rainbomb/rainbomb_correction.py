"""

Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import argparse
import os
import sys
import time
from typing import Any, Optional

import numpy as np
import pandas as pd
import xarray as xr

import climetlab as cml # type: ignore
import logging

# suppress cfgrib warnings about missing lat/lon coordinates when opening GRIB files with xarray
logging.getLogger('cfgrib').setLevel(logging.ERROR)


cur_folder = os.getcwd()

NETCDF_EXTENSIONS = {'.nc', '.nc4'}
GRIB_EXTENSIONS = {'.grib', '.grb', '.grib2', '.grb2'}

DATA_CONVERSION_FACTOR = 1000 # convert to mm

# file with the closest ERA5 neighbours for each point, used to check the rainbombs
FILENAME_CLOSEST_NEIGHBOURS = 'neighbours_era5_closest.nc'
# file with the intermediate and upper thresholds based on the SEAS5 data, used to check the rainbombs
FILENAME_THRESHOLDS = 'thresholds.csv'
# GRIB template for saving the corrected data, with the correct grid and fields but random date (to be set after with grib_set)
FILENAME_TEMPLATE = 'rain_template.grb'

DATA_KEY_REF = 'Ref'  # central point
DATA_KEY_MAX_NEIGH = 'MaxNeigh'  # max neighbour
DATA_KEY_REPLACE = 'Replace'  # replacement value if case of rainbomb
DATA_KEY_C_MAX = 'c_max'  # upper buffer threshold (Buffer2)
DATA_KEY_C_INTERM = 'c_interm'  # intermediate buffer threshold (Buffer1)
DATA_KEY_LOWER_BUFFER = 'LowerBuffer'  # lower buffer threshold (MaxNeigh + Buffer1)
DATA_KEY_RAINFALL_BIN = 'RainfallBin'  # rain bin based on MaxNeigh, used to assign the buffers
DATA_KEY_PRECIPITATION_VARIABLE = 'tp'  # variable name for precipitation in the dataset
DATA_KEY_NEIGHBOUR  = 'neighbour'
DATA_KEY_POINT_ID = 'point_id'
DATA_KEY_VALUES = 'values'

OUTPUT_FIRST_STEP_TIME = 0 # time to set in the GRIB output metadata (in hours)
OUTPUT_FIRST_STEP = 0 # first step to set in the GRIB output metadata
OUTPUT_LAST_STEP = 24 # last step to set in the GRIB output metadata


def correct_rainbomb(data: pd.Series) -> float:
    """
    Correct the rainbomb based on the max rainfall from neighbours, a replacement array
    (e.g. max/mean rain), and auxiliary information depending on the rain intensity.

    All methods analyzed are kept here (commented), but for operational purposes the cor3
    is used, which gave the best results.

    The function takes a dataframe row with the data needed for the correction, in different columns.
    The data refer to daily resolution, and the needed columns are:
    - Ref: daily precipitation of the central point
    - MaxNeigh: daily precipitation of the maximum neighbour for the central point
    - Replace: replacement value for the correction of the rainbomb. This is equal to MaxNeigh;
      a different argument for research purposes
    - c_max: upper buffer based on 99.99 percentile from long-term SEAS5 data
    - c_interm: intermediate buffer based on 99.50 percentile from long-term SEAS5 data

    The function reads the above data, and compares the Ref value with the MaxNeigh and the two buffers.
    If the rainfall excess of the Ref point is more than c_max mm higher than the MaxNeigh, all excess
    is removed as erroneous.
    If the excess is between the two buffers, the updated Ref value is linearly decreased considering
    Ref, MaxNeigh, c_max, and c_interm values.
    In all other cases Ref values stays intact.

    Parameters:
    -----------
    data: pandas.Series
        A row from a DataFrame containing the following columns:
        - Ref: daily precipitation of the central point
        - MaxNeigh: daily precipitation of the maximum neighbour
        - Replace: replacement value for the correction
        - c_max: upper buffer threshold (Buffer2)
        - c_interm: intermediate buffer threshold (Buffer1)

    Returns:
    --------
    float: The corrected rainfall value
    """
    # get values from the data
    ref_val = data[DATA_KEY_REF]  # central point
    max_val = data[DATA_KEY_MAX_NEIGH]  # max neighbour
    replacement_value = data[DATA_KEY_REPLACE]  # replacement value if case of rainbomb
    c_max = data[DATA_KEY_C_MAX]  # upper buffer threshold (Buffer2)
    c_interm = data[DATA_KEY_C_INTERM]  # intermediate buffer threshold (Buffer1)

    if ref_val > max_val + c_max:
        fore3 = replacement_value
    elif ref_val > max_val + c_interm:
        # weight based on relative difference of values
        wgt = (max_val + c_max - ref_val) / (c_max - c_interm)
        wgt = wgt.clip(0, 1)  # ensure weight is between 0-1
        fore3 = replacement_value + (ref_val - replacement_value) * wgt
    else:
        fore3 = ref_val

    return fore3


def read_era5(input_file: str, verbose: bool = False) -> tuple[xr.Dataset, str]:
    """
    Read ERA5 data from a NetCDF or GRIB file.
    Supports both formats based on file extension.
    Parameters:
    -----------
    input_file: str
        Path to the input ERA5 file (NetCDF or GRIB format)
    verbose: bool, optional
        If True, print the detected file format (default: False)
    Returns:
    --------
        tuple[xr.Dataset, str]: The loaded ERA5 dataset and its file extension
    Raises:
    -------
    RuntimeError: If the file format is unsupported or if cfgrib engine is required but not installed for GRIB files
    """
    # Detect file format based on extension
    file_ext = os.path.splitext(input_file)[1].lower()
    if verbose:
        print(f"Detected file extension: {file_ext}")

    # read daily raw field
    # Support both NetCDF and GRIB formats
    if file_ext in NETCDF_EXTENSIONS:
        if verbose:
            print("Opening as NetCDF file...")
        raw_era5 = xr.open_dataset(input_file)
    elif file_ext in GRIB_EXTENSIONS:
        if verbose:
            print("Opening as GRIB file...")

        try:
            raw_era5 = xr.open_dataset(input_file, engine='cfgrib')
        except ImportError:
            raise RuntimeError(
                "cfgrib engine required for GRIB files. Install with: pip install cfgrib"
            )
    else:
        # Try NetCDF as default (some files don't have extensions)
        if verbose:
            print("Unknown extension, attempting to open as NetCDF...")
        raw_era5 = xr.open_dataset(input_file)
    return raw_era5, file_ext


def correct_rainbomb_dataset(
    input_file: str,
    output_file: str,
    neighbours_file: Optional[str] = None,
    thresholds_file: Optional[str] = None,
    template_file: Optional[str] = None,
    parent_dir: Optional[str] = None,
    verbose: bool = False,
    set_grib_date_flag: bool = False,
) -> None:
    """
    Correct rainbomb artifacts in ERA5 daily precipitation data.

    A rainbomb is an unrealistically high rainfall value at a single grid point
    that is not supported by surrounding points. This function identifies and
    corrects such artifacts by comparing each grid point against its neighbours.

    Parameters:
    -----------
    input_file: str
        Path to the input ERA5 file (NetCDF format)
    output_file: str
        Path to the output corrected file (GRIB format)
    neighbours_file: str, optional
        Path to the neighbours data file (NetCDF). If not provided, defaults to
        'neighbours_era5_closest.nc' in parent_dir.
    thresholds_file: str, optional
        Path to the thresholds CSV file. If not provided, defaults to
        'thresholds.csv' in parent_dir.
    template_file: str, optional
        Path to the GRIB template file. If not provided, defaults to
        'template.grib' in parent_dir.
    parent_dir: str, optional
        Parent directory containing auxiliary data (used as fallback if individual
        files are not specified). Expected files:
        - 'neighbours_era5_closest.nc': closest ERA5 neighbours for each point
        - 'thresholds.csv': intermediate and upper thresholds based on SEAS5 data
        - 'template.grib': GRIB template for output formatting
    verbose: bool, optional
        If True, print progress information (default: False)
    set_grib_date_flag: bool, optional
        If True, set the correct date in the GRIB output file based on the input
        NetCDF time coordinate (default: False)

    Returns:
    --------
    None
    """
    # Validate input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(
            f"Input file not found: {input_file}. "
            "Please provide a valid path to an existing ERA5 NetCDF file."
        )

    # Resolve auxiliary file paths
    if parent_dir is None and neighbours_file is None and thresholds_file is None and template_file is None:
        raise ValueError(
            "Either parent_dir or individual auxiliary file paths must be provided"
        )

    # Set default paths based on parent_dir
    neighbours_file, thresholds_file, template_file = validate_config_file_paths(neighbours_file,
                                                                                 thresholds_file,
                                                                                 template_file,
                                                                                 parent_dir)

    if verbose:
        print(f"Reading input file: {input_file}")

    raw_era5, file_ext = read_era5(input_file, verbose)

    if verbose:
        print(f"Loading neighbours data from: {neighbours_file}")

    neighbours = xr.open_dataarray(neighbours_file) # type: ignore

    if verbose:
        print(f"Loading thresholds from: {thresholds_file}")

    thresholds_df = pd.read_csv(thresholds_file, index_col=0) # type: ignore

    # ----------- Correct daily field -----------
    if verbose:
        print("Processing rainfall data...")

    # final corrected data in the same format and units as the original dataset, to be saved as GRIB after
    corrected_ds = process_rainfall_data(raw_era5, neighbours, thresholds_df, verbose) 

    # ----------- Save as grb file -----------

    if verbose:
        print(f"Saving corrected data to: {output_file}")
        print(f"Using template file: {template_file}")

    # grb template to be used for saving the data, with the correct grid
    # and fields but random date (to be set after with grib_set)
    if file_ext in GRIB_EXTENSIONS:
        # If input is GRIB, use the same file as template to ensure correct grid and metadata
        template_file = input_file
    source = cml.load_source("file", template_file)
    template = source[0] # type: ignore

    template_step_range = template['stepRange']  # e.g., '5-6'
    template_end_step = template_step_range.split('-')[1]

    # Diagnostic: Print template info
    if verbose:
        print(f"Template missingValue: {template['missingValue']}")
        print(f"Template shape: {template.shape}")

    # Diagnostic: Check corrected data for issues
    corrected_data = corrected_ds[DATA_KEY_PRECIPITATION_VARIABLE].values.copy()
    if verbose:
        print(f"Corrected data dtype: {corrected_data.dtype}")
        print(f"Corrected data has NaN: {np.any(np.isnan(corrected_data))}")
        print(f"Corrected data has inf: {np.any(np.isinf(corrected_data))}")
        print(f"Corrected data min: {np.nanmin(corrected_data)}")
        print(f"Corrected data max: {np.nanmax(corrected_data)}")

    grib_metadata = {}

    # Extract the date from the input file to set correct date in GRIB output
    # Only run if the user specified --set-grib-date flag
    if set_grib_date_flag:
        # Handle time coordinate extraction for both NetCDF and GRIB formats
        time_val = extract_timestamp(raw_era5, file_ext)
        # Convert to string in YYYYMMDD format
        date_str = pd.to_datetime(time_val).strftime('%Y%m%d')

        grib_metadata = {
            'date': date_str,
            'time': OUTPUT_FIRST_STEP_TIME,
            'step': OUTPUT_LAST_STEP,
            'stepRange': f'{OUTPUT_FIRST_STEP}-{OUTPUT_LAST_STEP}',
        }

    # Get the corrected data and handle NaN values explicitly
    # This prevents the "failed to set key 'missingValue'" error that can occur
    # when climetlab tries to convert NaN to float32 max (3.4028235e+38)

    # Replace NaN and inf values with the template's missingValue
    template_missing_value = template['missingValue']
    if np.any(np.isnan(corrected_data)) or np.any(np.isinf(corrected_data)):
        if verbose:
            print(f"Replacing NaN/inf values with template missingValue: {template_missing_value}")
        corrected_data = np.where(np.isfinite(corrected_data), corrected_data, template_missing_value)

    with cml.new_grib_output(output_file, template=template) as output:
        # Use check_nans=False since we handled NaN values manually above
        output.write(corrected_data, check_nans=False, metadata=grib_metadata)

    if verbose:
        print("Correction complete!")


def process_rainfall_data(raw_era5: xr.Dataset, neighbours: xr.DataArray,
                          thresholds_df: pd.DataFrame, verbose: bool) -> xr.Dataset:
    """
    Process the rainfall data to identify and correct rainbombs.
    Parameters:
    -----------
    raw_era5: xr.Dataset
        The raw ERA5 dataset containing the precipitation variable
    neighbours: xr.DataArray
        DataArray containing the indices of the closest neighbours for each grid point
    thresholds_df: pd.DataFrame
        DataFrame containing the intermediate and upper thresholds based on SEAS5 data
    verbose: bool
        If True, print progress information
    Returns:
    --------    xr.Dataset
        A new Dataset with the corrected precipitation variable, in the same format and units as the original dataset
    """
    rain_mm = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE] * DATA_CONVERSION_FACTOR # convert to mm

    # max rain of neighbours (mask non-existing ones!)
    rain_neighbours_max = rain_mm.isel(values=neighbours).where(neighbours >= 0)

    # rename for compatibility among sets
    rain_neighbours_max = rain_neighbours_max.max(DATA_KEY_NEIGHBOUR).rename({DATA_KEY_POINT_ID: DATA_KEY_VALUES})

    # Mask potential rainbombs
    suspect_mask = rain_mm > rain_neighbours_max
    suspect_indices = np.where(suspect_mask.values)[0]
    if suspect_indices.size == 0:
        return raw_era5  # no postprocessing needed

    # keep another column for replace value, in case we want to use different than the max neighbour
    df = pd.DataFrame({
        DATA_KEY_REF: rain_mm.values[suspect_indices],
        DATA_KEY_MAX_NEIGH: rain_neighbours_max.values[suspect_indices],
        DATA_KEY_REPLACE: rain_neighbours_max.values[suspect_indices]
    }, index=suspect_indices)

    # keep only the cells whose value is over the max_rain + minimum buffer (Buffer2)
    # so we check the smallest number of points possible

    # assign each analyzed instance in rain bin based on MaxNeigh
    bins = np.digitize(df[DATA_KEY_MAX_NEIGH], thresholds_df.Limit, True)
    bins = np.clip(bins, 0, len(thresholds_df) - 1)
    # df[DATA_KEY_RAINFALL_BIN] = bins

    # get the buffers for each instance based on the rain bin
    df[DATA_KEY_C_INTERM] = thresholds_df.Buffer1.values[bins]
    df[DATA_KEY_C_MAX] = thresholds_df.Buffer2.values[bins]
    df[DATA_KEY_LOWER_BUFFER] = df[DATA_KEY_MAX_NEIGH] + df[DATA_KEY_C_INTERM]

    # keep only the suspicious instances where Ref - MaxNeigh > interm_buffer
    df = df[df[DATA_KEY_REF] > df[DATA_KEY_LOWER_BUFFER]]

    if verbose:
        print(f"Found {len(df)} rainbomb instances to correct")

    # final corrected dataset
    corrections = df.apply(correct_rainbomb, axis=1).values  # correct the checked instances
    corrected_np = rain_mm.values.copy()  # get a np array with the original ERA5 values
    corrected_np[df.index] = corrections  # replace the rainbombs with the corrected values

    # Convert back to original units
    corrected_values = corrected_np / DATA_CONVERSION_FACTOR

    # Validate corrected values before creating dataset
    if np.any(np.isnan(corrected_values)) or np.any(np.isinf(corrected_values)):
        # Replace invalid values with original rain_mm values (in original units)
        original_values = rain_mm.values / DATA_CONVERSION_FACTOR
        corrected_values = np.where(np.isfinite(corrected_values), corrected_values, original_values)

    # Create new dataset with same structure as original
    corrected_ds = xr.zeros_like(raw_era5)

    # Assign corrected precipitation data to the appropriate variable
    corrected_ds[DATA_KEY_PRECIPITATION_VARIABLE].values = corrected_values

    # Restore original attributes to preserve metadata (units, long_name, etc.)
    original_attrs = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE].attrs
    corrected_ds[DATA_KEY_PRECIPITATION_VARIABLE].attrs = original_attrs.copy()

    return corrected_ds


def extract_timestamp_valid_time(raw_era5: xr.Dataset, file_ext: str) -> pd.Timestamp:
    """Extract timestamp from a NetCDF or GRIB file.
    
    This function handles various formats of time coordinate extraction from xarray
    datasets, ensuring the returned value is always a pandas Timestamp.
    
    Parameters:
    -----------
    raw_era5: xr.Dataset
        The input xarray dataset (NetCDF or GRIB)
    file_ext: str
        The file extension to determine the file type
        
    Returns:
    --------
    pd.Timestamp
        The extracted timestamp
        
    Raises:
    -------
    ValueError
        If no valid timestamp can be extracted from the file
    """
    time_val: Any = None
    if file_ext in GRIB_EXTENSIONS:
        # GRIB files opened with cfgrib may have time in 'valid_time' or 'time'
        if 'valid_time' in raw_era5.attrs:
            time_val = raw_era5.attrs['valid_time']
        elif hasattr(raw_era5, 'valid_time') and raw_era5.valid_time is not None:
            time_val = raw_era5.valid_time.values if hasattr(raw_era5.valid_time, 'values') else raw_era5.valid_time
        elif 'time' in raw_era5:
            time_coord = raw_era5['time']
            time_val = time_coord.values if hasattr(time_coord, 'values') else time_coord
        else:
            # Try to get from the data variable
            if DATA_KEY_PRECIPITATION_VARIABLE in raw_era5:
                time_val = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE].attrs.get('valid_time', None)
                if time_val is None and hasattr(raw_era5[DATA_KEY_PRECIPITATION_VARIABLE], 'valid_time'):
                    time_val = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE].valid_time
                    time_val = time_val.values if hasattr(time_val, 'values') else time_val
            if time_val is None:
                raise ValueError("Could not extract time/date from input GRIB file. Please check the file structure.")
    else:
        # NetCDF files - use time coordinate
        time_coord = raw_era5.time
        if hasattr(time_coord, 'values'):
            # xarray DataArray - get the first time value
            time_values = time_coord.values
            if len(time_values) == 0:
                raise ValueError("Time coordinate is empty in input netCDF file.")
            time_val = time_values[0]
        else:
            # Already a datetime object
            time_val = time_coord
        if time_val is None:
            raise ValueError("Could not extract time/date from input netCDF file. Please check the file structure.")
    
    # Convert to pandas Timestamp if needed
    if isinstance(time_val, pd.Timestamp):
        return time_val
    
    # Handle numpy arrays - extract scalar value
    if isinstance(time_val, np.ndarray):
        if time_val.ndim == 0:
            # Scalar numpy array - extract the Python value
            time_val = time_val.item()
        else:
            # Multi-dimensional array - get first element
            time_val = time_val[0]
    
    # Handle xarray DataArray - extract value
    if isinstance(time_val, xr.DataArray):
        time_arr = time_val.values
        if isinstance(time_arr, np.ndarray) and time_arr.ndim == 0:
            time_val = time_arr.item()
        else:
            time_val = time_arr[0] if hasattr(time_arr, '__getitem__') else time_arr
    
    # Handle numpy datetime64
    if isinstance(time_val, np.datetime64):
        return pd.Timestamp(time_val)
    
    # Handle numpy scalars that might be datetime64
    if isinstance(time_val, np.generic):
        # Convert numpy scalar to Python native type
        time_val = time_val.item()
    
    # Final conversion to pd.Timestamp - handle datetime objects
    if isinstance(time_val, (pd.Timestamp, np.datetime64)):
        return pd.Timestamp(time_val)
    
    try:
        # Convert to datetime first, then to Timestamp
        return pd.Timestamp(pd.to_datetime(time_val))
    except (ValueError, TypeError, OSError) as e:
        raise ValueError(f"Could not convert time value '{time_val}' (type: {type(time_val).__name__}) to Timestamp: {e}")


def extract_timestamp(raw_era5: xr.Dataset, file_ext: str) -> pd.Timestamp:
    """Extract timestamp from a NetCDF or GRIB file.
    
    This function handles various formats of time coordinate extraction from xarray
    datasets, ensuring the returned value is always a pandas Timestamp.
    
    Parameters:
    -----------
    raw_era5: xr.Dataset
        The input xarray dataset (NetCDF or GRIB)
    file_ext: str
        The file extension to determine the file type
        
    Returns:
    --------
    pd.Timestamp
        The extracted timestamp
        
    Raises:
    -------
    ValueError
        If no valid timestamp can be extracted from the file
    """
    time_val: Any = None
    if file_ext in GRIB_EXTENSIONS:
        # GRIB files opened with cfgrib may have time in 'valid_time' or 'time'
        if 'time' in raw_era5:
            time_coord = raw_era5['time']
            time_val = time_coord.values if hasattr(time_coord, 'values') else time_coord
        elif 'valid_time' in raw_era5.attrs:
            time_val = raw_era5.attrs['valid_time']
        elif hasattr(raw_era5, 'valid_time') and raw_era5.valid_time is not None:
            time_val = raw_era5.valid_time.values if hasattr(raw_era5.valid_time, 'values') else raw_era5.valid_time
        else:
            # Try to get from the data variable
            if DATA_KEY_PRECIPITATION_VARIABLE in raw_era5:
                time_val = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE].attrs.get('valid_time', None)
                if time_val is None and hasattr(raw_era5[DATA_KEY_PRECIPITATION_VARIABLE], 'valid_time'):
                    time_val = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE].valid_time
                    time_val = time_val.values if hasattr(time_val, 'values') else time_val
            if time_val is None:
                raise ValueError("Could not extract time/date from input GRIB file. Please check the file structure.")
    else:
        # NetCDF files - use time coordinate
        time_coord = raw_era5.time
        if hasattr(time_coord, 'values'):
            # xarray DataArray - get the first time value
            time_values = time_coord.values
            if len(time_values) == 0:
                raise ValueError("Time coordinate is empty in input netCDF file.")
            time_val = time_values[0]
        else:
            # Already a datetime object
            time_val = time_coord
        if time_val is None:
            raise ValueError("Could not extract time/date from input netCDF file. Please check the file structure.")
    
    # Convert to pandas Timestamp if needed
    if isinstance(time_val, pd.Timestamp):
        return time_val
    
    # Handle numpy arrays - extract scalar value
    if isinstance(time_val, np.ndarray):
        if time_val.ndim == 0:
            # Scalar numpy array - extract the Python value
            time_val = time_val.item()
        else:
            # Multi-dimensional array - get first element
            time_val = time_val[0]
    
    # Handle xarray DataArray - extract value
    if isinstance(time_val, xr.DataArray):
        time_arr = time_val.values
        if isinstance(time_arr, np.ndarray) and time_arr.ndim == 0:
            time_val = time_arr.item()
        else:
            time_val = time_arr[0] if hasattr(time_arr, '__getitem__') else time_arr
    
    # Handle numpy datetime64
    if isinstance(time_val, np.datetime64):
        return pd.Timestamp(time_val)
    
    # Handle numpy scalars that might be datetime64
    if isinstance(time_val, np.generic):
        # Convert numpy scalar to Python native type
        time_val = time_val.item()
    
    # Final conversion to pd.Timestamp - handle datetime objects
    if isinstance(time_val, (pd.Timestamp, np.datetime64)):
        return pd.Timestamp(time_val)
    
    try:
        # Convert to datetime first, then to Timestamp
        return pd.Timestamp(pd.to_datetime(time_val))
    except (ValueError, TypeError, OSError) as e:
        raise ValueError(f"Could not convert time value '{time_val}' (type: {type(time_val).__name__}) to Timestamp: {e}")


def validate_config_file_paths(
    neighbours_file: Optional[str] = None,
    thresholds_file: Optional[str] = None,
    template_file: Optional[str] = None,
    parent_dir: Optional[str] = None,
) -> tuple[str, str, str]:
    """Validate and resolve auxiliary file paths for rainbomb correction.
    
    Parameters:
    -----------
    neighbours_file: Optional[str]
        Path to neighbours data file, or None to use default from parent_dir
    thresholds_file: Optional[str]
        Path to thresholds CSV file, or None to use default from parent_dir
    template_file: Optional[str]
        Path to GRIB template file, or None to use default from parent_dir
    parent_dir: Optional[str]
        Parent directory containing auxiliary files (used if individual paths are None)
    
    Returns:
    --------
    tuple[str, str, str]
        Tuple of (neighbours_file, thresholds_file, template_file) - all guaranteed to be str
    """
    if neighbours_file is None:
        neighbours_file = os.path.join(parent_dir, FILENAME_CLOSEST_NEIGHBOURS)  # type: ignore[arg-type]
    if thresholds_file is None:
        thresholds_file = os.path.join(parent_dir, FILENAME_THRESHOLDS)  # type: ignore[arg-type]
    if template_file is None:
        template_file = os.path.join(parent_dir, FILENAME_TEMPLATE)  # type: ignore[arg-type]

    # Validate auxiliary files exist
    if not os.path.isfile(neighbours_file):  # type: ignore[arg-type]
        raise FileNotFoundError(
            f"Neighbours file not found: {neighbours_file}. "
            "Please provide a valid path to the neighbours data file."
        )

    if not os.path.isfile(thresholds_file):  # type: ignore[arg-type]
        raise FileNotFoundError(
            f"Thresholds file not found: {thresholds_file}. "
            "Please provide a valid path to the thresholds CSV file."
        )

    if not os.path.isfile(template_file):  # type: ignore[arg-type]
        raise FileNotFoundError(
            f"Template file not found: {template_file}. "
            "Please provide a valid path to the GRIB template file."
        )
        
    return neighbours_file, thresholds_file, template_file  # type: ignore[return-value]


def getargs(argv=sys.argv) -> argparse.Namespace:
    """Get program arguments.

    Returns:
    --------
    args: argparse.Namespace
        Namespace of program arguments
    """
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Script for correcting ERA5 single-grid rainbombs for daily fields.

        A rainbomb is an unrealistically high rainfall value at a single grid point
        that is not supported by surrounding points. This utility identifies and
        corrects such artifacts by comparing each grid point against its neighbours
        and applying threshold-based corrections.
        """,
        prog=prog,
    )
    parser.add_argument(
        "-d",
        "--parent_dir",
        type=str,
        default=f'{cur_folder}/config',
        help="Parent directory where the auxiliary data for the rainbomb correction are located. "
             "Ignored if individual auxiliary files are specified.",
    )
    parser.add_argument(
        "-n",
        "--neighbours_file",
        type=str,
        default=None,
        help="Path to the neighbours data file (NetCDF). "
             f"If not provided, defaults to '{FILENAME_CLOSEST_NEIGHBOURS}' in parent_dir.",
    )
    parser.add_argument(
        "-t",
        "--thresholds_file",
        type=str,
        default=None,
        help="Path to the thresholds CSV file. "
             f"If not provided, defaults to '{FILENAME_THRESHOLDS}' in parent_dir.",
    )
    parser.add_argument(
        "-m",
        "--template_file",
        type=str,
        default=None,
        help="Path to the GRIB template file. "
             f"If not provided, defaults to '{FILENAME_TEMPLATE}' in parent_dir.",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        required=True,
        help="Input file for correction (raw ERA5 in NetCDF or GRIB format)",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Output file (corrected ERA5 in GRIB format)",
    )
    parser.add_argument(
        "--set-grib-date",
        action="store_true",
        default=False,
        help="Set the correct date in the GRIB output file based on the input NetCDF or GRIB time coordinate",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Print progress information",
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    """Main function for running the rainbomb correction script."""

    args = getargs(argv)

    # Validate arguments
    if args.parent_dir is None and args.neighbours_file is None and args.thresholds_file is None and args.template_file is None:
        raise ValueError(
            "Either --parent_dir or at least one individual auxiliary file "
            "(--neighbours_file, --thresholds_file, --template_file) must be provided"
        )

    try:
        start_time = time.perf_counter()

        if args.verbose:
            print("Starting rainbomb correction...")
            print(f"Input file: {args.input_file}")
            print(f"Output file: {args.output_file}")
            if args.parent_dir:
                print(f"Parent directory: {args.parent_dir}")
            if args.neighbours_file:
                print(f"Neighbours file: {args.neighbours_file}")
            if args.thresholds_file:
                print(f"Thresholds file: {args.thresholds_file}")
            if args.template_file:
                print(f"Template file: {args.template_file}")

        correct_rainbomb_dataset(
            input_file=args.input_file,
            output_file=args.output_file,
            neighbours_file=args.neighbours_file,
            thresholds_file=args.thresholds_file,
            template_file=args.template_file,
            parent_dir=args.parent_dir,
            verbose=args.verbose,
            set_grib_date_flag=args.set_grib_date,
        )

        elapsed_time = time.perf_counter() - start_time
        if args.verbose:
            print(f"Time elapsed: {elapsed_time:0.2f} seconds")
    except FileNotFoundError as fnfe:
        print(f'ERROR: {fnfe}')
    except Exception as e:
        raise RuntimeError(f'{e}')


def main_script():
    sys.exit(main())


if __name__ == "__main__":
    main_script()