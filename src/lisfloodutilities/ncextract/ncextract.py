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
import glob
import os
import pandas as pd
import sys
import time
import xarray as xr

def extract(inputcsv, directory, outputfile, nc):
    start_time = time.perf_counter()

    if not os.path.isfile(inputcsv):
        print(f'ERROR: {inputcsv} is missing!')
        sys.exit(0)

    if not os.path.isdir(directory):
        print(f'ERROR: {directory} is missing or not a directory!')
        sys.exit(0)

    filepaths = []
    for file in glob.glob(os.path.join(directory, '*.nc')):
        print(file)
        filepaths.append(file)

    if not filepaths:
        print(f'ERROR: No NetCDF files found in {directory}')
        sys.exit(0)

    try:
        poi = pd.read_csv(inputcsv)
        poi_indexer = poi.set_index('id')[['lat', 'lon']].to_xarray()
    except:
        print('ERROR: Please check that CSV file is formatted correctly!')
        sys.exit(0)

    # chunks is set to auto for general purpose processing
    # it could be optimized depending on input NetCDF
    ds = xr.open_mfdataset(filepaths, chunks='auto', parallel=True)

    ds_poi = ds.sel(lat=poi_indexer.lat, lon=poi_indexer.lon, method='nearest')
    print(ds_poi)
    print("Processing...")

    if nc:
        ds_poi.to_netcdf(outputfile)
    else:
        df = ds_poi.to_dataframe()
        df_reset = df.reset_index()
        df_reset.to_csv(outputfile, index=False)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Time elapsed: {elapsed_time:0.2f} seconds")

def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to extract values from (multiple) NetCDF files
        at specific coordinates (over a time series).
        Coordinates of points of interest must be 
        included in a CSV file with at least 3 columns
        named id, lat, lon.
        """,
        prog=prog,
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV file (id, lat, lon)")
    parser.add_argument("-d", "--directory", required=True, help="Input directory with .nc files")
    parser.add_argument("-o", "--output", required=True, help="Output file (default is CSV, use -nc for NetCDF)")
    parser.add_argument("-nc", "--nc", action='store_true', help='Output to NetCDF')

    args = parser.parse_args()

    extract(args.input, args.directory, args.output, args.nc)

def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()