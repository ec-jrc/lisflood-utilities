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
import subprocess
import sys
import time
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

import climetlab as cml # type: ignore


cur_folder = os.path.dirname(os.path.realpath(__file__))

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


def set_grib_date(output_file: str, init_date: str, step: int = 24, verbose: bool = False) -> None:
    """
    Set the correct date in a GRIB file using grib_set.

    The GRIB file created from a template may have random date fields.
    This function uses the ecCodes grib_set tool to update the dataDate,
    dataTime, and step in the GRIB file.

    Produces a new GRIB file named 'era5_corrected_{init_date}.grb' with the updated date fields.

    Parameters:
    -----------
    output_file: str
        Path to the GRIB file to modify
    init_date: str
        Date in YYYYMMDD format (e.g., '20200101')
    step: int, optional
        Forecast step in hours (default: 24)
    verbose: bool, optional
        If True, print the grib_set command being executed (default: False)

    Returns:
    --------
    None

    Raises:
    -------
    RuntimeError: If grib_set command fails
    """
    output_file_parent = os.path.dirname(output_file)
    corrected_file = os.path.join(output_file_parent, f'era5_corrected_{init_date}.grb')
    cmd = ['grib_set', '-s', f'step={step},dataDate={init_date},dataTime=0', output_file, corrected_file]

    if verbose:
        print(f"Executing: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        if verbose and result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"grib_set failed: {e.stderr}") from e
    except FileNotFoundError:
        raise RuntimeError(
            "grib_set not found. Please ensure ecCodes is installed and grib_set is in PATH"
        ) from None


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
    # Resolve auxiliary file paths
    if parent_dir is None and neighbours_file is None and thresholds_file is None and template_file is None:
        raise ValueError(
            "Either parent_dir or individual auxiliary file paths must be provided"
        )

    # Set default paths based on parent_dir
    if neighbours_file is None:
        neighbours_file = os.path.join(parent_dir, FILENAME_CLOSEST_NEIGHBOURS) # type: ignore
    if thresholds_file is None:
        thresholds_file = os.path.join(parent_dir, FILENAME_THRESHOLDS) # type: ignore
    if template_file is None:
        template_file = os.path.join(parent_dir, FILENAME_TEMPLATE) # type: ignore

    if verbose:
        print(f"Reading input file: {input_file}")

    # read daily raw field
    raw_era5 = xr.open_dataset(input_file)

    if verbose:
        print(f"Loading neighbours data from: {neighbours_file}")

    neighbours = xr.open_dataarray(neighbours_file) # type: ignore

    if verbose:
        print(f"Loading thresholds from: {thresholds_file}")

    thresholds_df = pd.read_csv(thresholds_file, index_col=0) # type: ignore

    # ----------- Correct daily field -----------
    if verbose:
        print("Processing rainfall data...")

    rain_mm = raw_era5[DATA_KEY_PRECIPITATION_VARIABLE] * DATA_CONVERSION_FACTOR # convert to mm

    # max rain of neighbours (mask non-existing ones!)
    rain_neighbours_max = rain_mm.isel(values=neighbours).where(neighbours >= 0)

    # rename for compatibility among sets
    rain_neighbours_max = rain_neighbours_max.max(DATA_KEY_NEIGHBOUR).rename({DATA_KEY_POINT_ID: DATA_KEY_VALUES})

    # keep another column for replace value, in case we want to use different than the max neighbour
    initial_check = pd.DataFrame({
        DATA_KEY_REF: rain_mm.isel(values=(rain_mm > rain_neighbours_max)),
        DATA_KEY_MAX_NEIGH: rain_neighbours_max.isel(values=(rain_mm > rain_neighbours_max)),
        DATA_KEY_REPLACE: rain_neighbours_max.isel(values=(rain_mm > rain_neighbours_max))
    }, index=rain_mm.isel(values=(rain_mm > rain_neighbours_max))[DATA_KEY_VALUES].values)

    # keep only the cells whose value is over the max_rain + minimum buffer (Buffer2)
    # so we check the smallest number of points possible

    # assign each analyzed instance in rain bin based on MaxNeigh
    thresholds_row = np.digitize(initial_check[DATA_KEY_MAX_NEIGH], thresholds_df.Limit, True)
    thresholds_row[thresholds_row >= len(thresholds_df)] = len(thresholds_df) - 1  # if over last row's value, go last row
    initial_check[DATA_KEY_RAINFALL_BIN] = thresholds_row

    # get the buffers for each instance based on the rain bin
    initial_check[DATA_KEY_C_INTERM] = thresholds_df.Buffer1.iloc[thresholds_row].values
    initial_check[DATA_KEY_C_MAX] = thresholds_df.Buffer2.iloc[thresholds_row].values
    initial_check[DATA_KEY_LOWER_BUFFER] = initial_check[DATA_KEY_MAX_NEIGH] + thresholds_df.Buffer1.iloc[thresholds_row].values

    # keep only the suspicious instances where Ref - MaxNeigh > interm_buffer
    initial_check = initial_check[initial_check[DATA_KEY_REF] > initial_check[DATA_KEY_LOWER_BUFFER]]

    if verbose:
        print(f"Found {len(initial_check)} rainbomb instances to correct")

    # final corrected dataset
    corrections = [i[1] for i in list(initial_check.iterrows())]  # get the data for corrections (i[0] is the index)
    corrections = [correct_rainbomb(i) for i in corrections]  # correct the checked instances
    corrected_np = rain_mm.values.copy()  # get a np array with the original ERA5 values
    corrected_np[initial_check.index] = np.array(corrections)  # replace the rainbombs with the corrected values
    corrected_ds = raw_era5 * 0 + corrected_np / DATA_CONVERSION_FACTOR  # final corrected data in the same format and units as the original data

    # ----------- Save as grb file -----------

    if verbose:
        print(f"Saving corrected data to: {output_file}")

    # grb template to be used for saving the data
    template_file = os.path.join(parent_dir, FILENAME_TEMPLATE) # type: ignore
    source = cml.load_source("file", template_file)
    template = source[0]

    with cml.new_grib_output(output_file, template=template) as output:
        output.write(corrected_ds[DATA_KEY_PRECIPITATION_VARIABLE].values)

    # Extract the date from the input file to set correct date in GRIB output
    # Only run if the user specified --set-grib-date flag
    if set_grib_date_flag:
        # The input NetCDF has a time coordinate we can use
        time_coord = raw_era5.time
        if hasattr(time_coord, 'values'):
            # xarray DataArray - get the first time value
            time_val = time_coord.values[0]
        else:
            # Already a datetime object
            time_val = time_coord

        # Convert to string in YYYYMMDD format
        date_str = pd.to_datetime(time_val).strftime('%Y%m%d')

        # Set the correct date in the GRIB file using grib_set
        set_grib_date(output_file, date_str, step=24, verbose=verbose)

    if verbose:
        print("Correction complete!")


def getargs():
    """Get program arguments.

    Returns:
    --------
    args: argparse.Namespace
        Namespace of program arguments
    """
    prog = os.path.basename(sys.argv[0])
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
        help="Input file for correction (raw ERA5 in NetCDF format)",
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
        help="Set the correct date in the GRIB output file based on the input NetCDF time coordinate",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Print progress information",
    )
    args = parser.parse_args(sys.argv[1:]) # parser.parse_args()
    return args


def main(argv=sys.argv):
    """Main function for running the rainbomb correction script."""
    
    args = getargs()

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

    except Exception as e:
        raise RuntimeError(f'{e}')


def main_script():
    sys.exit(main())


if __name__ == "__main__":
    main_script()