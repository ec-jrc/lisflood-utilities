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
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm


PRODUCT_ERA5 = 'era5'
PRODUCT_SEASONAL = 'seasonal'
PRODUCTS = {PRODUCT_ERA5, PRODUCT_SEASONAL}


def getargs(argv=sys.argv):
    """Get program arguments.

    Parameters:
    -----------
    argv : list
        Command line arguments (default: sys.argv)

    Returns:
    --------
    args: argparse.Namespace
        Namespace of program arguments
    """
    parser = argparse.ArgumentParser(
        description="""
        Script for generating neighbour indices for each grid point.

        For reducing the time of the computations for the rainbomb correction,
        this script keeps the indices of the neighbours for each grid point.
        Thus, the correction script can iterate directly over each point and its
        neighbours. The neighbours on east/west are always the closest points in
        the same latitude. For the north and south, the search is done to the next
        upper/lower latitude and there are 3 methods used:
        - Keep all points that lie within the longitudinal range of west-east based
          on the west & east neighbours
        - Keep the closest north and south neighbour (can be 1 or even 2 if same
          distance north/south-west of the point and north/south-east)
        - Keep all points that have smaller distance than the maximum distance from
          the reference point to its west/east neighbour
        """,
    )
    parser.add_argument(
        "-p",
        "--product",
        type=str,
        required=False,
        default=f"{PRODUCT_ERA5}",
        help=f"Product to be used: {','.join(PRODUCTS)}; the default is {PRODUCT_ERA5}",
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        required=True,
        help="Work directory containing the auxiliary data files and the output files. "
             "The directory should contain the sample file '<product>_sample.grb' depending on the product used.",
    )
    args = parser.parse_args(argv[1:])  # assign namespace to args
    return args


def get_rect_neighbours_ilat(lat_index, sample_df, lats_all):
    """
    Gets the neighbours of each point (west/east at same latitude, and all points in upper/lower lat within the west-east rectangle)
    The function gives the neighbours for all points in the same latitude field
    Max neighbours north/south are 3 (this was checked for all points)
    """
    lat_used = lats_all[lat_index] # get actual latitude
    points_in_same_lat = sample_df.query('latitude==@lat_used') # subset all points within the latitude field

    all_neighbours_lat = []
    # loop through all points and find their neighbours
    for i_point_index in range(len(points_in_same_lat)): 

        i_point = points_in_same_lat.iloc[i_point_index]
        
        ##### ----------- get the west and east neighbour
        same_lat_neighbours = points_in_same_lat.longitude - i_point.longitude
        
        if (same_lat_neighbours>0).sum()>0:
            east_neighbour = same_lat_neighbours.where(same_lat_neighbours>0).idxmin()
        else:
            east_neighbour = same_lat_neighbours.idxmin()
            
        if (same_lat_neighbours<0).sum()>0:
            west_neighbour = same_lat_neighbours.where(same_lat_neighbours<0).idxmax()
        else:
            west_neighbour = same_lat_neighbours.idxmax()

        ##### -----------  get the upper neighbours (they should be within the lon range of the neighbours in the same latitude)
        if lat_index!=0:
            upper_lat = lats_all[lat_index-1]
            points_in_upper_lat = sample_df.query('latitude==@upper_lat')
            
            criterio1 = points_in_upper_lat.longitude<=sample_df.iloc[east_neighbour, 1]
            criterio2 = points_in_upper_lat.longitude>=sample_df.iloc[west_neighbour, 1]

            # if point is at the east edge so the east limit is negative and the west positive, it needs OR, otherwise AND
            edge_multiplication = sample_df.iloc[east_neighbour].longitude*sample_df.iloc[west_neighbour].longitude
            if edge_multiplication<0 and sample_df.iloc[east_neighbour].longitude<0:
                criterion_final = criterio1|criterio2
            else:
                criterion_final = criterio1&criterio2

            upper_neighbour = points_in_upper_lat.where(criterion_final).dropna().sort_values('longitude').index
            upper_neighbour = upper_neighbour.tolist()
            
            # just in case not all 3 points are available, fill the missing with NaN (flagged as -1)
            upper_neighbour += [-1]*3
            upper_neighbour = upper_neighbour[:3]
        else:
            upper_neighbour = [-1]*3 # if already in the first line, we can't go further up, so all flagged with -1

        ##### ----------- get the lower neighbours (they should be within the lon range of the neighbours in the same latitude)
        if lat_index!=len(lats_all)-1:
            lower_lat = lats_all[lat_index+1]
            points_in_lower_lat = sample_df.query('latitude==@lower_lat')
            
            criterio1 = points_in_lower_lat.longitude<=sample_df.iloc[east_neighbour, 1]
            criterio2 = points_in_lower_lat.longitude>=sample_df.iloc[west_neighbour, 1]

            edge_multiplication = sample_df.iloc[east_neighbour].longitude*sample_df.iloc[west_neighbour].longitude
            if edge_multiplication<0 and sample_df.iloc[east_neighbour].longitude<0:
                criterion_final = criterio1|criterio2
            else:
                criterion_final = criterio1&criterio2

            lower_neighbour = points_in_lower_lat.where(criterion_final).dropna().sort_values('longitude').index
            lower_neighbour = lower_neighbour.tolist()
            
            # just in case not all 3 points are available, fill the missing with NaN (flagged as -1)
            lower_neighbour += [-1]*3
            lower_neighbour = lower_neighbour[:3]
        else:
            lower_neighbour = [-1]*3 # if already in the last line, we can't go further down, so all flagged with -1
            
        ##### ----------- concatenate all neighbours of the point in a dataarray
        all_neighbours_i_point = xr.DataArray(upper_neighbour+[west_neighbour]+[east_neighbour]+lower_neighbour,
                                  coords={'neighbour': ['north1', 'north2', 'north3', 'west', 'east', 'south1', 'south2', 'south3']})
        all_neighbours_i_point = all_neighbours_i_point.assign_coords(point_id=i_point.name).expand_dims('point_id') # add point id
        all_neighbours_lat.append(all_neighbours_i_point)

    all_neighbours_lat = xr.concat(all_neighbours_lat, dim='point_id')  # concatenate to 1 xarray object
    return all_neighbours_lat


def get_four_neighbours_ilat(lat_index, sample_df, lats_all):
    """
    Gets the neighbours of each point (west/east at same latitude, and north/south)
    The function gives the neighbours for all points in the same latitude field.
    Max north/south neighbours are 2 (if symmetric from the reference point)
    """

    lat_used = lats_all[lat_index] # get actual latitude
    points_in_same_lat = sample_df.query('latitude==@lat_used')  # subset all points within the latitude field

    all_neighbours_lat = []
    # loop through all points and find their neighbours
    for i_point_index in range(len(points_in_same_lat)): 

        i_point = points_in_same_lat.iloc[i_point_index]
        
        ##### ----------- get the west and east neighbour
        same_lat_neighbours = points_in_same_lat.longitude - i_point.longitude
        
        if (same_lat_neighbours>0).sum()>0:
            east_neighbour = same_lat_neighbours.where(same_lat_neighbours>0).idxmin()
        else:
            east_neighbour = same_lat_neighbours.idxmin()
            
        if (same_lat_neighbours<0).sum()>0:
            west_neighbour = same_lat_neighbours.where(same_lat_neighbours<0).idxmax()
        else:
            west_neighbour = same_lat_neighbours.idxmax()

        ##### -----------  get the north neighbours
        if lat_index!=0:
            upper_lat = lats_all[lat_index-1]
            points_in_upper_lat = sample_df.query('latitude==@upper_lat')

            north_neighbour = np.abs(points_in_upper_lat.longitude - i_point.longitude)
            north_neighbour = north_neighbour.where(north_neighbour==north_neighbour.min()).dropna()
            north_neighbour = north_neighbour.index.tolist()
            north_neighbour += [-1]
            north_neighbour = north_neighbour[:2]
        else:
            north_neighbour = [-1]*2  # if already in the first line, we can't go further up, so flagged with -1

        ##### ----------- get the south neighbours
        if lat_index!=len(lats_all)-1:
            lower_lat = lats_all[lat_index+1]
            points_in_lower_lat = sample_df.query('latitude==@lower_lat')
            
            south_neighbour = np.abs(points_in_lower_lat.longitude - i_point.longitude)
            south_neighbour = south_neighbour.where(south_neighbour==south_neighbour.min()).dropna()
            south_neighbour = south_neighbour.index.tolist()
            south_neighbour += [-1]
            south_neighbour = south_neighbour[:2]
        else:
            south_neighbour = [-1]*2  # if already in the last line, we can't go further down, so flagged with -1
            
        ##### ----------- concatenate all neighbours of the point in a dataarray
        all_neighbours_i_point = xr.DataArray(north_neighbour+[west_neighbour]+[east_neighbour]+south_neighbour,
                                      coords={'neighbour': ['north1', 'north2', 'west', 'east', 'south1', 'south2']})
        all_neighbours_i_point = all_neighbours_i_point.assign_coords(point_id=i_point.name).expand_dims('point_id') # add point id
        all_neighbours_lat.append(all_neighbours_i_point)

    all_neighbours_lat = xr.concat(all_neighbours_lat, dim='point_id') # concatenate to 1 xarray object
    return all_neighbours_lat


def main(argv=sys.argv):
    """Main function for running the generate_neighbours script."""

    args = getargs(argv)

    # Validate arguments
    product = args.product
    if product not in PRODUCTS:
        print(f'product should be one of: {PRODUCTS}')
        sys.exit(1)

    parent_dir = args.directory
    sample = xr.open_dataset(f'{parent_dir}/{product}_sample.grb')

    # keep only the values dims
    dims_to_drop = list(sample.dims)
    dims_to_drop = [i for i in dims_to_drop if i!='values']
    sample = sample.isel({i:0 for i in dims_to_drop})

    # keep only the lat, lon coordinates
    coords_to_drop = list(sample.coords)
    coords_to_drop = [i for i in coords_to_drop if i not in ['latitude', 'longitude']]
    sample = sample.drop_vars(coords_to_drop)

    # convert longitudes to [-180, 180] from [0, 360]
    lons = sample.longitude.values
    lons = np.where(lons>180, lons-360, lons)
    sample.longitude[:] = lons

    # get a dataframe with all the coordinates
    sample_df = sample.to_dataframe()[['latitude', 'longitude']]

    lats_all = sorted(sample_df.latitude.unique())[::-1]

    # Create partial function with sample_df and lats_all
    rect_func = partial(get_rect_neighbours_ilat, sample_df=sample_df, lats_all=lats_all)

    all_indices = range(len(lats_all))
    pool = Pool(processes=16)
    all_neighboursrect = list(tqdm(pool.imap(rect_func, all_indices), total=len(all_indices), position=0))
    all_neighboursrect = xr.concat(all_neighboursrect, dim='point_id').sortby('point_id')
    all_neighboursrect.name = 'neighbours'
    encoding = {'neighbours': {'zlib': True, 'complevel': 4}}
    all_neighboursrect.to_netcdf(f'{parent_dir}/neighbours_{product}_rectangle.nc', encoding=encoding)

    #### Keep only the closest points in North/South
    four_func = partial(get_four_neighbours_ilat, sample_df=sample_df, lats_all=lats_all)

    all_indices = range(len(lats_all))
    pool = Pool(processes=16)
    all_neighbours4 = list(tqdm(pool.imap(four_func, all_indices), total=len(all_indices), position=0))
    all_neighbours4 = xr.concat(all_neighbours4, dim='point_id').sortby('point_id')
    all_neighbours4.name = 'neighbours'
    encoding = {'neighbours': {'zlib': True, 'complevel': 4}}
    all_neighbours4.to_netcdf(f'{parent_dir}/neighbours_{product}_closest.nc', encoding=encoding)


def main_script():
    """Entry point for the script when run as a module."""
    sys.exit(main())


if __name__ == "__main__":
    main_script()
