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
import pandas as pd
import os
import sys
import time
import xarray as xr
import cfgrib
from pathlib import Path
from typing import Union, Optional

def read_points(inputcsv: Union[str, Path]) -> xr.Dataset:
    """It reads a CSV file with coordinates of points: gauging stations, reservoirs...

    Parameters:
    -----------
    inputcsv: string or pathlib.Path
        a CSV indicating the points for which the time series will be extracted. It must contain three columns: the identification of the point, and its two coordinates

    Returns:
    --------
    poi_xr: xarray.Dataset
        the input data converted into a Xarray object
    """

    if not os.path.isfile(inputcsv):
        print(f'ERROR: {inputcsv} is missing!')
        sys.exit(1)

    try:
        poi_df = pd.read_csv(inputcsv)
        original_columns = poi_df.columns.copy()
        poi_df.columns = original_columns.str.lower()
        # find columns representing coordinates
        coord_1 = [col for col in poi_df.columns if col.startswith('lon') or col.startswith('x')][0]
        coord_2 = [col for col in poi_df.columns if col.startswith('lat') or col.startswith('y')][0]
        # find the column representing the point ID
        idx_col = poi_df.columns.difference([coord_1, coord_2])[0]
        # convert to xarray.Dataset
        poi_xr = poi_df.set_index(idx_col)[[coord_1, coord_2]].to_xarray()
        rename_dim = {idx_col: col for col in original_columns if col.lower() == idx_col}
        poi_xr = poi_xr.rename(rename_dim)
    except:
        print('ERROR: Please check that the CSV file is formatted correctly!')
        sys.exit(2)

    return poi_xr



def read_inputmaps(directory: Union[str, Path]) -> xr.Dataset:
    """It extract from a series of input files (NetCDF or GRIB) the time series of a set of points

    Parameters:
    -----------
   directory: string or pathlib.Path
        the directory containing the input files, which can be either in NetCDF or GRIB format

    Returns:
    --------
    ds: xarray.dataset
        containing the concatenation of all the input maps
    """

    pattern_engine = {'*.nc': 'netcdf4',
                      '*.grib': 'cfgrib'}

    if not os.path.isdir(directory):
        print(f'ERROR: {directory} is missing or not a directory!')
        sys.exit(1)
    else:
        directory = Path(directory)

    filepaths = []
    for pattern, engine in pattern_engine.items():
        filepaths = list(directory.glob(pattern))
        if len(filepaths) > 0:
            break

    if not filepaths:
        print(f'ERROR: No NetCDF/GRIB file found in {directory}')
        sys.exit(2)

    print(f'{len(filepaths)} input {engine} file(s) found in "{directory}"')

    # chunks is set to auto for general purpose processing
    # it could be optimized depending on input NetCDF
    ds = xr.open_mfdataset(filepaths, engine=engine, chunks='auto', parallel=True)

    return ds



def extract_timeseries(poi: xr.Dataset,
                       ds: xr.Dataset,
                       outputfile: Optional[Union[str, Path]] = None
                      ) -> Optional[xr.Dataset]:
    """It extract from a series of input files (NetCDF or GRIB) the time series of a set of points

    Parameters:
    -----------
    inputcsv: xarray.Dataset
        a Dataset indicating the coordinates of the points of interest. It must have only two variables (the coordinates), and the names of this variables must be dimensions in "ds"
    directory: xarray.Dataset
        the time stack of input maps from which the time series will be extracted
    ouputfile: optional, string or pathlib.Path
        the file where the results will be saved. It can be either a CSV or a NetCDF file

    Returns:
    --------
    By default, it puts out an xarray.Dataset with the extracted time series. Instead, if "outputfile" is provided, results will be saved to a file (CSV or NetCDF)
    """

    coord_1, coord_2 = list(poi)
    if not all(coord in ds.coords for coord in [coord_1, coord_2]):
        print(f'ERROR: The variables in "poi" (coordinates) are not coordinates in "ds"')
        sys.exit(1)

    # extract time series
    ds_poi = ds.sel({coord_1: poi[coord_1], coord_2: poi[coord_2]}, method='nearest')

    if outputfile is None:
        return ds_poi.compute()

    else:
        outputfile = Path(outputfile)
        if outputfile.suffix == '.nc':
            ds_poi.to_netcdf(outputfile)
        elif outputfile.suffix == '.csv':
            df = ds_poi.to_dataframe()
            df_reset = df.reset_index()
            df_reset.to_csv(outputfile, index=False)
        else:
            print('ERROR: the extension of the output file must be either ".nc" or ".csv"')
            sys.exit(2)
        print(f'Results exported as {outputfile}')



def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to extract time series of values from
        (multiple) NetCDF files at specific coordinates.
        Coordinates of points of interest must be
        included in a CSV file with at least 3 columns
        named id, lat, lon.
        """,
        prog=prog,
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV file (id, lat, lon)")
    parser.add_argument("-d", "--directory", required=True, help="Input directory with .nc files")
    parser.add_argument("-o", "--output", required=True, help="Output file. Two extensions are supported: .csv or .nc")

    args = parser.parse_args()

    try:
        start_time = time.perf_counter()

        print('Reading input CSV...')
        points = read_points(args.input)
        print('Reading input maps...')
        maps = read_inputmaps(args.directory)
        print(maps)
        print('Processing...')
        extract_timeseries(points, maps, args.output)

        elapsed_time = time.perf_counter() - start_time
        print(f"Time elapsed: {elapsed_time:0.2f} seconds")

    except Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1)

def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()
