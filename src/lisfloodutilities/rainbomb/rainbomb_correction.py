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
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

import climetlab as cml


def correct_rainbomb(data):
    """Correct the rainbomb based on the max rainfall from neighbours, a replacement array
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
    ref_val = data['Ref']  # central point
    max_val = data['MaxNeigh']  # max neighbour
    replacement_value = data['Replace']  # replacement value if case of rainbomb
    c_max = data['c_max']  # upper buffer threshold (Buffer2)
    c_interm = data['c_interm']  # intermediate buffer threshold (Buffer1)

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


def correct_rainbomb_dataset(
    input_file: str,
    output_file: str,
    neighbours_file: Optional[str] = None,
    thresholds_file: Optional[str] = None,
    template_file: Optional[str] = None,
    parent_dir: Optional[str] = None,
    verbose: bool = False
) -> None:
    """Correct rainbomb artifacts in ERA5 daily precipitation data.

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
        'rain_template.grb' in parent_dir.
    parent_dir: str, optional
        Parent directory containing auxiliary data (used as fallback if individual
        files are not specified). Expected files:
        - neighbours_era5_closest.nc: closest ERA5 neighbours for each point
        - thresholds.csv: intermediate and upper thresholds based on SEAS5 data
        - rain_template.grb: GRIB template for output formatting
    verbose: bool, optional
        If True, print progress information (default: False)

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
        neighbours_file = os.path.join(parent_dir, 'neighbours_era5_closest.nc') # type: ignore
    if thresholds_file is None:
        thresholds_file = os.path.join(parent_dir, 'thresholds.csv') # type: ignore
    if template_file is None:
        template_file = os.path.join(parent_dir, 'rain_template.grb') # type: ignore

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

    rain_mm = raw_era5.tp * 1000  # convert to mm

    # max rain of neighbours (mask non-existing ones!)
    rain_neighbours_max = rain_mm.isel(values=neighbours).where(neighbours >= 0)

    # rename for compatibility among sets
    rain_neighbours_max = rain_neighbours_max.max('neighbour').rename({'point_id': 'values'})

    # keep another column for replace value, in case we want to use different than the max neighbour
    initial_check = pd.DataFrame({
        'Ref': rain_mm.isel(values=(rain_mm > rain_neighbours_max)),
        'MaxNeigh': rain_neighbours_max.isel(values=(rain_mm > rain_neighbours_max)),
        'Replace': rain_neighbours_max.isel(values=(rain_mm > rain_neighbours_max))
    }, index=rain_mm.isel(values=(rain_mm > rain_neighbours_max))['values'].values)

    # keep only the cells whose value is over the max_rain + minimum buffer (Buffer2)
    # so we check the smallest number of points possible

    # assign each analyzed instance in rain bin based on MaxNeigh
    thresholds_row = np.digitize(initial_check.MaxNeigh, thresholds_df.Limit, True)
    thresholds_row[thresholds_row >= len(thresholds_df)] = len(thresholds_df) - 1  # if over last row's value, go last row
    initial_check['RainfallBin'] = thresholds_row

    # get the buffers for each instance based on the rain bin
    initial_check['c_interm'] = thresholds_df.Buffer1.iloc[thresholds_row].values
    initial_check['c_max'] = thresholds_df.Buffer2.iloc[thresholds_row].values
    initial_check['LowerBuffer'] = initial_check.MaxNeigh + thresholds_df.Buffer1.iloc[thresholds_row].values

    # keep only the suspicious instances where Ref - MaxNeigh > interm_buffer
    initial_check = initial_check[initial_check.Ref > initial_check.LowerBuffer]

    if verbose:
        print(f"Found {len(initial_check)} rainbomb instances to correct")

    # final corrected dataset
    corrections = [i[1] for i in list(initial_check.iterrows())]  # get the data for corrections (i[0] is the index)
    corrections = [correct_rainbomb(i) for i in corrections]  # correct the checked instances
    corrected_np = rain_mm.values.copy()  # get a np array with the original ERA5 values
    corrected_np[initial_check.index] = np.array(corrections)  # replace the rainbombs with the corrected values
    corrected_ds = raw_era5 * 0 + corrected_np / 1000  # final corrected data in the same format as the original data

    # ----------- Save as grb file -----------

    if verbose:
        print(f"Saving corrected data to: {output_file}")

    # grb template to be used for saving the data
    template_file = os.path.join(parent_dir, 'rain_template.grb') # type: ignore
    source = cml.load_source("file", template_file)
    template = source[0]

    with cml.new_grib_output(output_file, template=template) as output:
        output.write(corrected_ds['tp'].values)

    if verbose:
        print("Correction complete!")

    # Note: the above python line uses a template with random grib fields in the date,
    # so we need to set the correct ones. This can be done after this script with
    # the grib_set function when changing the dataDate to the correct one:
    # grib_set -s step=24,dataDate=$idate,dataTime=0 $output era5_corrected_${idate}.grb


def getarg():
    """Get program arguments.

    Returns:
    --------
    args: argparse.Namespace
        Namespace of program arguments
    """
    parser = argparse.ArgumentParser(
        description="""
        Script for correcting ERA5 single-grid rainbombs for daily fields.

        A rainbomb is an unrealistically high rainfall value at a single grid point
        that is not supported by surrounding points. This utility identifies and
        corrects such artifacts by comparing each grid point against its neighbours
        and applying threshold-based corrections.
        """,
        prog=os.path.basename(sys.argv[0]),
    )
    parser.add_argument(
        "-d",
        "--parent_dir",
        type=str,
        required=False,
        default="/home/ecm8227/Repositories/glofas_forcings/version5/rainbomb_correction/",
        help="Parent directory where the auxiliary data for the rainbomb correction are located.",
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
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Print progress information",
    )
    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """Main function for running the rainbomb correction script."""
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
        default=None,
        help="Parent directory where the auxiliary data for the rainbomb correction are located. "
             "Ignored if individual auxiliary files are specified.",
    )
    parser.add_argument(
        "-n",
        "--neighbours_file",
        type=str,
        default=None,
        help="Path to the neighbours data file (NetCDF). "
             "If not provided, defaults to 'neighbours_era5_closest.nc' in parent_dir.",
    )
    parser.add_argument(
        "-t",
        "--thresholds_file",
        type=str,
        default=None,
        help="Path to the thresholds CSV file. "
             "If not provided, defaults to 'thresholds.csv' in parent_dir.",
    )
    parser.add_argument(
        "-m",
        "--template_file",
        type=str,
        default=None,
        help="Path to the GRIB template file. "
             "If not provided, defaults to 'rain_template.grb' in parent_dir.",
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
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Print progress information",
    )

    args = parser.parse_args(argv[1:])

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
            verbose=args.verbose
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