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
import numpy as np
import pandas as pd
import os
import sys
import time
import xarray as xr
import cfgrib
from pathlib import Path
from typing import Union, Optional, Tuple
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
        x_coord = [col for col in poi_df.columns if col.startswith('lon') or col.startswith('x')][0]
        y_coord = [col for col in poi_df.columns if col.startswith('lat') or col.startswith('y')][0]

        # find the column representing the point ID
        idx_col = poi_df.columns.difference([x_coord, y_coord])[0]

        # convert to xarray.Dataset
        poi_xr = poi_df.set_index(idx_col)[[x_coord, y_coord]].to_xarray()
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



def read_ldd(
    file: Union[str, Path],
) -> xr.DataArray:
    """Reads the local drainage direction map

    Parameters:
    -----------
    file: string or pathlib.Path
        the NetCDF file of the local drainage direction map

    Returns:
    --------
    ldd: xarray.DataArray
    """
    
    return xr.open_dataset(file)['Band1']
    

def rename_geographic_coords(
    target: Union[xr.Dataset, xr.DataArray],
    reference: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
    """Renames the geographical coordinates/variables in the target dataset to match those in the reference dataset.
    
    Parameters:
    -----------
    target: xarray.Dataset or xarray.DataArray
        Object whose geographical coordinates will be renamed
    target: xarray.Dataset or xarray.DataArray
        Reference names of the geographical coordinates
        
    Returns:
    --------
    A similar object as 'target', but with the names of the geographical coordinates in 'reference'
    """
    
    # names of geographical coordinates in the target dataset
    try:
        x_obj = [coord for coord in target.coords if coord.startswith('lon') or coord.startswith('x')][0]
    except:
        x_obj = [coord for coord in target.variables if coord.startswith('lon') or coord.startswith('x')][0]
    try:
        y_obj = [coord for coord in target.coords if coord.startswith('lat') or coord.startswith('y')][0]
    except:
        y_obj = [coord for coord in target.variables if coord.startswith('lat') or coord.startswith('y')][0]
    
    # names of geographical coordinates in the reference dataset
    try:
        x_ref = [coord for coord in reference.coords if coord.startswith('lon') or coord.startswith('x')][0]
    except:
        x_ref = [coord for coord in reference.variables if coord.startswith('lon') or coord.startswith('x')][0]
    try:
        y_ref = [coord for coord in reference.coords if coord.startswith('lat') or coord.startswith('y')][0]
    except:
        y_ref = [coord for coord in reference.variables if coord.startswith('lat') or coord.startswith('y')][0]
    
    return target.rename({x_obj: x_ref, y_obj: y_ref})



def find_inflow_points(
    lat: float,
    lon: float,
    ldd: xr.DataArray
) -> xr.Dataset:
    """This function finds the upstream coordinates of the pixels flowing into the input coordinates
    
    Parameteres:
    ------------
    lat: float
        latitude of the input point
    lon: float
        longitued of the input point
    ldd: xarray.DataArray
        map of local drainage directions
        
    Returns:
    --------
    points: xarra.Dataset
        Contains the coordinates of the pixels flowing into the point of interest.
    """
    
    # Determine coordinate system
    y_coord = [coord for coord in ldd.coords if coord.startswith('lat') or coord.startswith('y')][0]
    x_coord = [coord for coord in ldd.coords if coord.startswith('lon') or coord.startswith('x')][0]

    # spatial resolution of the input map
    resolution = np.round(np.mean(np.diff(ldd[x_coord].values)), 4)

    # Define window around the input pixel
    window = 1.5 * resolution
    ldd_window = ldd.sel({y_coord: slice(lat + window, lat - window),
                          x_coord: slice(lon - window, lon + window)})
    # 2D arrays of the coordinates of the pixels in the window
    lons, lats = np.meshgrid(ldd_window[x_coord].data, ldd_window[y_coord].data)

    # create a 1D mask of inflow pixels
    inflow = np.array([[3, 2, 1],
                       [6, 5, 4],
                       [9, 8, 7]])
    mask = (ldd_window == inflow).data.flatten()

    # apply mask
    lons, lats = lons.flatten()[mask], lats.flatten()[mask]
    
    # convert to xarray.Dataset
    points = pd.DataFrame(data={y_coord: lats, x_coord: lons})
    points.index.name = 'inflow'
    points = points.to_xarray()

    return points



def extract_timeseries(
    maps: xr.Dataset,
    poi: xr.Dataset,
    inflow: bool = False,
    ldd: Optional[xr.Dataset] = None,
    output: Optional[Union[str, Path]] = None,
    overwrite: bool = False
) -> Optional[xr.Dataset]:
    """Extracts time series for each point of interest and saves them separately.

    Parameters:
    -----------
    maps: xarray.Dataset
        the time stack of input maps from which the time series will be extracted
    poi: xarray.Dataset
        A Dataset indicating the coordinates of the points of interest. It must have only two variables (the coordinates), and the names of this variables must be dimensions in "maps"
    inflow: boolean
        Wheter to extract the value in that pixel (False) or compute the sum of the pixels flowing into it (True)
    ldd: optional, xarray.Dataset
        Map of local drainage directions. Only needed if 'inflow' is True
    output: optional, string or pathlib.Path
        The directory where the results will be saved. If not provided, returns an xarray.Dataset.
    overwrite: boolean
        whether to overwrite or skip points of interest whose output file already exists. By default is False

    Returns:
    --------
    If 'output' is None, returns an xarray.Dataset with extracted time series.
    Otherwise, saves results in a set of NetCDF files.
    """
        
    if "id" not in poi.coords:
        raise ValueError('ERROR: "poi" must contain an "id" coordinate.')
        
    coord_1, coord_2 = list(poi)
    if not all(coord in maps.coords for coord in [coord_1, coord_2]):
        raise ValueError(f'The variables in "poi" (coordinates) are not coordinates in "maps"')
    if coord_1.startswith('y') or coord_1.startswith('lat'):
        coords = [coord_1, coord_2]
    else:
        coords = [coord_2, coord_1]
        
    if inflow and ldd is None:
        raise ValueError('An "ldd" map must be provided if the option "inflow" is enabled.')

    # create output directory
    if output:
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)

    maps_poi = []
    for ID in poi.id.data:
        
        if output:
            output_file = output / f'{ID}.nc'
            if output_file.exists() and not overwrite:
                print(f'Output file {output_file} already exists. Moving forward to the next point.')
                continue

        if inflow:
            # find pixels flowing into the point of interest
            inflows = find_inflow_points(*[poi.sel(id=ID)[coord].item() for coord in coords], ldd)
            
            # extract and sum the time series of the pixels flowing into the point
            series = maps.sel({coord: inflows[coord] for coord in coords}, method='nearest').sum('inflow')
            series = series.expand_dims({'id': [ID]})            
            series = series.assign_coords({coord: ('id',
                                                   [maps[coord].sel({coord: poi.sel(id=ID)[coord].item()}, method='nearest').item()])
                                                   # [poi.sel(id=ID)[coord].item()])
                                           for coord in coords})
            series.attrs.update(maps.attrs)

        else:
            # extract time series of the point
            series = maps.sel({coord: poi.sel(id=ID)[coord].item() for coord in coords}, method='nearest')
            series = series.expand_dims(dim={'id': [ID]})

        # save time series
        if output is None:
            maps_poi.append(series)
            print('{0} | Time series for point {1} extracted'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                                     ID))
        else:           
            series.to_netcdf(output_file)               
            print('{0} | Time series for point {1} saved in {2}'.format(datetime.now().strftime('%Y-%m-%d %H:%M'),
                                                                       ID,
                                                                       output_file))
    
    if output is None:
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
    parser.add_argument("-d", "--dir", required=True, help="Input directory with NetCDF or GRIB files")
    parser.add_argument("-o", "--output", required=True, help="Output directory for time series")
    parser.add_argument("-s", "--start", type=str, default=None, help='Start datetime (YYYY-MM-DD) (default: None)')
    parser.add_argument("-e", "--end", type=str, default=None, help='End datetime (YYYY-MM-DD) (default: None)')
    parser.add_argument("-i", "--inflow", action="store_true", default=False, help='Extract the aggregation of pixels flowing into the points of interest')
    parser.add_argument("-l", "--ldd", required=False, help="Map of local drainage directions. Only neccesary if 'inflow' is True")
    parser.add_argument("-w", "--overwrite", action="store_true", default=False, help="Overwrite existing output files")

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
        
        print('Reading input maps...')
        maps = read_inputmaps(args.dir, start=args.start, end=args.end)
        print(maps)
        
        print('Reading input CSV...')
        points = read_points(args.points)
        points = rename_geographic_coords(points, maps)
        
        if args.inflow:
            if args.ldd is None:
                raise ValueError("-ldd must be provided if --inflow is enabled.")
            else:
                print('Reading the LDD map...')
                ldd = read_ldd(args.ldd)
                args.ldd = rename_geographic_coords(ldd, maps)

        print('Processing...')
        extract_timeseries(maps, points, args.inflow, args.ldd, args.output, args.overwrite)

        elapsed_time = time.perf_counter() - start_time
        print(f"Time elapsed: {elapsed_time:0.2f} seconds")

    except Exception as e:
        raise RuntimeError(f'{e}')
        sys.exit(1)

        
        
def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()