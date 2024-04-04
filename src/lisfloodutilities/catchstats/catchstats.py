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
from pathlib import Path
import pandas as pd
import sys
import time
import xarray as xr
from typing import List, Union, Optional, Literal
from tqdm.auto import tqdm


def catchment_statistics(inputmaps: Union[str, Path],
                         mask: Union[str, Path],
                         statistic: List[Literal['mean', 'sum', 'std', 'var', 'min', 'max', 'median', 'count']], 
                         output: Union[str, Path],
                         pixarea: Optional[str] = None,
                         overwrite: bool = False
                        ):
    """
    Given a set of input maps and catchment masks, it computes catchment statistics.
    
    Parameters:
    -----------
    inputmaps: str or pathlib.Path
        directory that contains the input NetCDF files whose statistics will be computed. These files can be static (withouth time dimenion) or dynamic (with time dimension)
    mask: str or pathlib.Path
        directory that contains the NetCDF files that define the catchment boundaries. These files can be the output of the `cutmaps` tool
    statistic: list of strings
        statistics to be computed. Only some statistics are available: 'mean', 'sum', 'std', 'var', 'min', 'max', 'median', 'count'
    output: str or pathlib.Path
        directory where the resulting NetCDF files will be saved.
    pixarea: optional or str
        if provided, a NetCDF file with pixel area used to compute weighted statistics. It is specifically meant for geographic projection systems where the area of a pixel varies with latitude
    overwrite: boolean
        whether to overwrite or skip catchments whose output NetCDF file already exists. By default is False, so the catchment will be skipped
    
    Returns:
    --------
    A NetCDF file will be created for each catchment in the directory "mask"
    """

    start_time = time.perf_counter()
    
    # output directory
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    
    # check statistic
    possible_stats = ['mean', 'sum', 'std', 'var', 'min', 'max', 'median', 'count']
    assert all(stat in possible_stats for stat in statistic), "All values in 'statistic' should be one of these: {0}".format(', '.join(possible_stats))
    
    # input maps
    if not os.path.isdir(inputmaps):
        print(f'ERROR: {inputmaps} is missing or not a directory!')
        sys.exit(0)
    else:
        inputmaps = Path(inputmaps)
        filepaths = list(inputmaps.glob('*.nc'))
        if not filepaths:
            print(f'ERROR: No NetCDF files found in "{inputmaps}"')
            sys.exit(0)
        else:
            print(f'{len(filepaths)} input NetCDF files found in "{inputmaps}"')
        try:
            # chunks is set to auto for general purpose processing
            # it could be optimized depending on input NetCDF
            ds = xr.open_mfdataset(filepaths, chunks='auto', parallel=True)
        except:
            # for static maps
            ds = xr.Dataset({file.stem.split('_')[0]: xr.open_dataset(file)['Band1'] for file in filepaths})
        if 'wgs_1984' in ds:
            ds = ds.drop_vars('wgs_1984')
    
    # catchment masks
    if not os.path.isdir(mask):
        print(f'ERROR: {mask} is missing or not a directory!')
        sys.exit(0)
    else:
        mask = Path(mask)
        maskpaths = list(mask.glob('*.nc'))
        if not maskpaths:
            print(f'ERROR: No NetCDF files found in "{mask}"')
            sys.exit(0)
        else:
            maskpaths = {int(file.stem): file for file in maskpaths}
            print(f'{len(maskpaths)} mask NetCDF files found in "{mask}"')
        
    # weighing map
    if pixarea is not None:
        if not os.path.isfile(pixarea):
            print(f'ERROR: {pixarea} is missing!')
            sys.exit(0)
        else:
            try:
                weight = xr.open_dataset(pixarea)['Band1']
            except:
                print(f'ERROR: The weighing map "{pixarea}" could not be loaded')
                sys.exit(0) 

    # define coordinates and variables of the resulting Dataset
    dims = dict(ds.dims)
    dimnames = [dim.lower() for dim in dims]
    if 'lat' in dimnames and 'lon' in dimnames:
        x_dim, y_dim = 'lon', 'lat'
    else:
        x_dim, y_dim = 'x', 'y'
    del dims[x_dim]
    del dims[y_dim]
    coords = {dim: ds[dim] for dim in dims}
    stats_dict = {var: statistic for var in ds}
    variables = [f'{var}_{stat}' for var, stats in stats_dict.items() for stat in stats]
    
    # compute statistics for each catchemnt
    for ID in tqdm(maskpaths.keys(), desc='processing catchments'):
        
        fileout = output / f'{ID:04}.nc'
        if fileout.exists() &  ~overwrite:
            print(f'Output file {fileout} already exists. Moving forward to the next catchment')
            continue
        
        # create empty Dataset
        coords.update({'id': [ID]})
        ds_aoi = xr.Dataset({var: xr.DataArray(coords=coords, dims=coords.keys()) for var in variables})
        
        # read mask map
        try:
            maskpath = maskpaths[ID]
            aoi = xr.open_dataset(maskpath)['Band1']
            aoi = xr.where(aoi.notnull(), 1, aoi)
        except:
            print(f'ERROR: The mask {maskpath} could not be read')
            continue
            
        # apply mask to the dataset
        masked_ds = ds.sel({x_dim: aoi[x_dim], y_dim: aoi[y_dim]}).where(aoi == 1)
        masked_ds = masked_ds.compute()

        # apply weighting by pixel area
        if pixarea is not None:
            masked_weight = weight.sel({x_dim: aoi[x_dim], y_dim: aoi[y_dim]}).where(aoi == 1)
            weighted_ds = masked_ds.weighted(masked_weight.fillna(0))  

        # compute statistics
        for var, stats in stats_dict.items(): 
            for stat in stats:
                if stat in ['mean', 'sum', 'std', 'var']:
                    if pixarea is not None:
                        x = getattr(weighted_ds, stat)(dim=[x_dim, y_dim])[var]
                    else:
                        x = getattr(masked_ds, stat)(dim=[x_dim, y_dim])[var]
                elif stat in ['min', 'max', 'median', 'count']:
                    x = getattr(masked_ds, stat)(dim=[x_dim, y_dim])[var]
                ds_aoi[f'{var}_{stat}'].loc[{'id': ID}] = x
        
        # export
        ds_aoi.to_netcdf(fileout)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Time elapsed: {elapsed_time:0.2f} seconds")

    
    
def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to compute catchment statistics from
        (multiple) NetCDF files.
        The mask map is a NetCDF file with values in the
        area of interest and NaN elsewhere.
        """,
        prog=prog,
    )
    parser.add_argument("-i", "--input", required=True, help="Directory containint the input NetCDF files")
    parser.add_argument("-m", "--mask", required=True, help="Directory containing the mask NetCDF files")
    parser.add_argument("-s", "--statistic", nargs='+', required=True, help='List of statistics to be computed')
    parser.add_argument("-o", "--output", required=True, help="Directory where the output NetCDF files will be saved")
    parser.add_argument("-a", "--area", required=False, default=None, help="NetCDF file of pixel area used to weight the statistics")
    parser.add_argument("-W", "--overwrite", action="store_true", help="Overwrite existing output files")
    
    args = parser.parse_args()
    
    catchment_statistics(args.input, args.mask, args.statistic, args.output, args.area, args.overwrite)

    
    
def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()