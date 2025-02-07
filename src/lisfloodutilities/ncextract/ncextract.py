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
from datetime import datetime



def read_points(inputcsv: Union[str, Path]) -> xr.Dataset:
    """Reads a CSV file with the coordinates of the points of interest: gauging stations, reservoirs...

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
        raise FileNotFoundError(f'{inputcsv} does not exist!')

    try:
        # read input CSV
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
        poi_xr = poi_xr.rename({idx_col: 'id'})
    except:
        raise ValueError(f"Could not read CSV properly. Please check the format.\nDetails: {e}")

    return poi_xr



def read_inputmaps(
    directory: Union[str, Path],
    start: datetime = None,
    end: datetime = None
) -> xr.Dataset:
    """Reads a set of input files (NetCDF or GRIB) and selects the period of interest

    Parameters:
    -----------
    directory: string or pathlib.Path
        the directory containing the input files, which can be either in NetCDF or GRIB format
    start: datetime
        Start of the extraction period
    end: datetime
        End of the extraction period

    Returns:
    --------
    maps: xarray.dataset
        containing the concatenation of all the input maps
    """

    pattern_engine = {
        '*.nc': 'netcdf4',
        '*.grib': 'cfgrib'
    }

    if not os.path.isdir(directory):
        raise FileNotFoundError(f'{directory} is missing or not a directory!')
    
    directory = Path(directory)
    filepaths = []
    for pattern, engine in pattern_engine.items():
        filepaths = list(directory.glob(pattern))
        if filepaths:
            break

    if not filepaths:
        raise FileNotFoundError(f'No NetCDF/GRIB files found in {directory}')

    print(f'{len(filepaths)} input {engine} file(s) found in "{directory}"')
    
    try:
        # load dataset
        maps = xr.open_mfdataset(filepaths, engine=engine, chunks='auto', parallel=True)
        # Note: chunks is set to auto for general purpose processing
        #       it could be optimized depending on input NetCDF
    except Exception as e:
        raise RuntimeError(f'Failed to open datasets using engine "{engine}": {str(e)}')
    
    # Validate start and end dates
    if start and not isinstance(start, datetime):
        raise ValueError(f"The 'start' parameter must be a datetime object, got {type(start)}")
    if end and not isinstance(end, datetime):
        raise ValueError(f"The 'end' parameter must be a datetime object, got {type(end)}")
    
    if start or end:
        time_dims = [dim for dim in maps.dims if 'time' in dim.lower()]
        if not time_dims:
            print('WARNING: No time dimension found, skipping time filtering')
            return maps
        try:
            maps = maps.sel({time_dims[0]: slice(start, end)})
        except Exception as e:
            raise ValueError(f'Failed to apply time filter: {str(e)}')

    return maps



def extract_timeseries(
    maps: xr.Dataset,
    poi: xr.Dataset,
    output_dir: Optional[Union[str, Path]] = None,
    output_format: str = 'nc'
) -> Optional[xr.Dataset]:
    """Extracts time series for each point of interest and saves them separately.

    Parameters:
    -----------
    maps: xarray.Dataset
        the time stack of input maps from which the time series will be extracted
    poi: xarray.Dataset
        A Dataset indicating the coordinates of the points of interest. It must have only two variables (the coordinates), and the names of this variables must be dimensions in "maps"
    output_dir: optional, string or pathlib.Path
        The directory where the results will be saved. If not provided, returns an xarray.Dataset.
    output_format: optional, str
        The format of output files: "nc" (NetCDF) or "csv". Default is "nc".

    Returns:
    --------
    If 'output_dir' is None, returns an xarray.Dataset with extracted time series.
    Otherwise, saves results in the specified format.
    """
    
    if "id" not in poi.coords:
        raise ValueError('ERROR: "poi" must contain an "id" coordinate.')
        
    coord_1, coord_2 = list(poi)
    if not all(coord in maps.coords for coord in [coord_1, coord_2]):
        raise ValueError(f'The variables in "poi" (coordinates) are not coordinates in "maps"')

    # create output directory
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # check that the output_format is correct
        if output_format not in ['csv', 'nc']:
            raise ValueError('the extension of the output file must be either ".nc" or ".csv"')

    maps_poi = []
    for ID in poi.id.data:
        # extract time series of the point
        series = maps.sel(
            {coord_1: poi.sel(id=ID)[coord_1].item(),
             coord_2: poi.sel(id=ID)[coord_2].item()},
            method='nearest'
        )
        series = series.expand_dims(dim={'id': [ID]})

        # save time series
        if output_dir is None:
            maps_poi.append(series)
            print('{0} | Time series for point {1} extracted'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%s'),
                                                                      ID))
        else:           
            output_file = output_dir / f'{ID}.{output_format}'
            if output_format == 'nc':
                series.to_netcdf(output_file)
            elif output_format == 'csv':
                df = series.to_dataframe().reset_index().drop(['id', coord_1, coord_2], axis=1)
                df.to_csv(output_file, index=False)                
            print('{0} | Time series for point {1} saved in {2}'.format(datetime.now().strftime('%Y-%m-%d %H:%M'),
                                                                       ID,
                                                                       output_file))
    
    if output_dir is None:
        return xr.concat(maps_poi, dim='id').compute()
    
    return None



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
    parser.add_argument("-p", "--points", required=True, help="CSV file of points of interest (id, lat, lon)")
    parser.add_argument("-i", "--input", required=True, help="Input directory with .nc files")
    parser.add_argument("-o", "--output", required=True, help="Output directory for time series")
    parser.add_argument("-f", "--format", choices=["nc", "csv"], default="nc", help="Output format: 'nc' of 'csv' (default: 'nc')")
    parser.add_argument("-s", "--start", type=str, default=None, help='Start datetime (YYYY-MM-DD)')
    parser.add_argument("-e", "--end", type=str, default=None, help='End datetime (YYYY-MM-DD)')

    args = parser.parse_args()

    # parse dates
    if args.start:
        try:
            args.start = datetime.strptime(args.start, "%Y-%m-%d")
        except ValueError:
            raise ValueError("Invalid date format in the 'start' argument. Use 'YYYY-MM-DD'.")
    if args.end:
        try:
            args.end = datetime.strptime(args.end, "%Y-%m-%d")
        except ValueError:
            raise ValueError("Invalid date format in the 'end' argument. Use 'YYYY-MM-DD'.")
                    
    try:

        start_time = time.perf_counter()
        
        print('Reading input CSV...')
        points = read_points(args.points)
        
        print('Reading input maps...')
        maps = read_inputmaps(args.input, start=args.start, end=args.end)
        print(maps)
        
        print('Processing...')
        extract_timeseries(maps, points, args.output, args.format)

        elapsed_time = time.perf_counter() - start_time
        print(f"Time elapsed: {elapsed_time:0.2f} seconds")

    except Exception as e:
        raise RuntimeError(f'{e}')
        sys.exit(1)

def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()
