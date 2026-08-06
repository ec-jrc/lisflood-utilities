#   This is a script to create a mask with mild sloping river pixels.
#   It takes LISFLOOD channels slope map (changrad.nc), the LDD (ldd.nc) and mask of the catchment/domain.
#   Pixels where river slope < threshold are added to the mask, if drainage area is large enough.
#   A minimum number of consecutive mild sloping downstream pixels is required for a pixel to be added to the mask.
#   Uses xarray and numpy instead of PCRaster.
#
#   Usage:
#   mctrivers.py -i changrad.nc -l ec_ldd.nc -m mask.nc -u upArea.nc -E y x -S 0.001 -N 5 -U 500 -O chanmct.nc

import tempfile
import os
from typing import List, Tuple
import xarray as xr
from earthkit.hydro import distance, move, upstream
from earthkit.hydro.river_network import repair as river_network_repair
from lisfloodutilities.cutmaps.helpers import get_river_network_from_map
from .helpers import copy_coordinates_and_attributes, X_COORD_NAMES, Y_COORD_NAMES


def getarg():
    """ Get program arguments.

    :return: args:  namespace of program arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--changradfile', type=str, required=True,
                        help='Changrad file with channels riverbed bed slope (chagrad.nc)')
    parser.add_argument('-l', '--LDDfile', type=str, required=True,
                        help='LISFLOOD LDD file (ldd.nc)')
    parser.add_argument('-u', '--uparea', type=str, required=True,
                        help='Upstream area file (upArea.nc)')
    parser.add_argument('-m', '--maskfile', type=str, required=False, default='',
                        help='Mask or domain file (mask.nc)')
    parser.add_argument('-S', '--slope', type=float, required=False, default=0.001,
                        help='Slope threshold to use MCT (default slp < 0.001)')
    parser.add_argument('-N', '--nloops', type=int, required=False, default=5, choices=range(0, 100),
                        help='Number of consecutive downstream MCT gridcells to be an MCT cell (default = 5)')
    parser.add_argument('-U', '--minuparea', type=float, required=False, default=0,
                        help='Minimum upstream area (same units as in the -u file) for including a cell in the final Muskingum mask (default 0)')
    parser.add_argument('-E', '--coordsnames', type=str, nargs='+', required=False, default=[],
                        help='Coords names for lat, lon (in this order with space!) from the netcdf files used')
    parser.add_argument('-O', '--outputfilename', type=str, required=False, default="chanmct.nc",
                        help='Output file (chanmct.nc)')
    args = parser.parse_args()  # assign namespace to args
    return args

    
def mct_mask(channels_slope_file: str, ldd_file: str, uparea_file: str, mask_file: str = '',
             slp_threshold: float = 0.001, nloops: int = 5, minuparea: float = 0,
             coords_names: List[str] = []) -> xr.DataArray:
    """
    
    Builds a mask of mild sloping rivers for use in LISFLOOD with MCT diffusive river routing. It takes LISFLOOD channels slope map (changrad.nc), the LDD (ldd.nc),
    the upstream drained area map (upArea.nc) and the catchment/domain mask (mask.nc), and outputs a bolean mask (chanmct.nc). Pixels where riverbed gradient < threshold
    (slp_threshold) are added to the mask if their drainage area is large enough (minuparea) and they also have at least nloops consecutive downstream pixels that meet
    the same condition for slope (drainage area will be met as downstream the area increases).

    Usage:
    The tool requires the following input arguments:
    
    channels_slope_file: LISFLOOD channels gradient map (changrad.nc)
    ldd_file: LISFLOOD local drain direction file (ldd.nc)
    uparea_file: LISFLOOD Upstream area file (upArea.nc)
    mask_file: a mask nc file; if not given (default) all cells are considered valid.
    slp_threshold: Riverbed slope threshold to use MCT diffusive wave routing (default: 0.001)
    nloops: Number of consecutive downstream grid cells that also need to comply with the slope requirement for including a grid cell in the MCT rivers mask (default: 5)
    minuparea: Minimum upstream drainage area for a pixel to be included in the MCT rivers mask (uses the same units as in the -u file) (default: 0)
    coords_names: Coordinates names for lat, lon (in this order as list) used in the the netcdf files (default: []; checks for commonly used names ['x', 'lon', 'rlon'], similar for lat names)
    outputfile: Output file containing the rivers mask where LISFLOOD can use the MCT diffusive wave routing (default: chanmct.nc)
    
    Example for generating an MCT rivers mask with pixels where riverbed slope < 0.001, drainage area > 500 kms and at least 5 downstream pixels meet the same
    two conditions, considering the units of the upArea.nc file are given in kms:
    
    mct_mask(channels_slope_file='changrad.nc', ldd_file='ldd.nc', uparea_file='upArea.nc', mask_file='mask.nc',
             slp_threshold=0.001, nloops=5, minuparea=500, coords_names=['y' , 'x'])
    """

    # ---------------- Read LDD to get coordinates and repair it before creating network
    ldd_file_to_use = ldd_file
    LD_ds = xr.open_dataset(ldd_file_to_use)
    x_proj, y_proj = extract_coords(LD_ds, coords_names)
    LD_ds.close()

    # Create river network from original LDD to use in the lddrepair function in case of failure.
    # This is needed because if there are cycles in the LDD, the downstream and path functions will fail, 
    # and without the river network we can't run the lddrepair function to fix the LDD.;
    # This will help to ensure that the repaired LDD is hydrologically consistent
    try:
        network = get_river_network_from_map(ldd_file_to_use)
    except Exception as e:
        print(f"Error occurred while creating river network the LDD might have cycles. Trying to repair it: {e}")
        repaired_ldd_file = ldd_file.replace('.nc', '_repaired_by_mctrivers.nc')
        river_network_repair(ldd_file_to_use, repaired_ldd_file, river_network_format='pcr_d8', input_source='file')
        ldd_file_to_use = repaired_ldd_file
        print(f'LDD repaired successfuly. Switch to repaired file: {repaired_ldd_file}')
        # Create river network from repaired LDD
        network = get_river_network_from_map(repaired_ldd_file, export=False)

    # ---------------- Process channels slope netcdf
    CH = xr.open_dataset(channels_slope_file)
    CH = prepare_dataset(CH, x_proj, y_proj, 'changrad')

    # get number of rows and columns
    rows, cols = CH.sizes[y_proj], CH.sizes[x_proj]

    # get coords of map
    x_all = CH.variables[x_proj].values
    y_all = CH.variables[y_proj].values

    # ---------------- Create rivers mask (Step 1: Flag all grid cells with slope lower than threshold)
    rivers_mask = (CH.changrad < slp_threshold).astype(int)
    CH.close()

    # ---------------- Read upstream area
    UA = xr.open_dataset(uparea_file)
    UA = prepare_dataset(UA, x_proj, y_proj, 'domain')
    UA = UA['domain']  # Extract as Dataarray

    # check that the area is over the minimum
    minarea_bool = (UA >= minuparea).astype(int)
    UA.close()

    # ---------------- Read domain/basin (mask) area
    # Create domain mask from LDD - valid cells are those with valid LDD (1-9)
    # First, reload LDD to create the domain mask
    LD_for_mask_ds = xr.open_dataset(ldd_file_to_use)
    LD_for_mask = prepare_dataset(LD_for_mask_ds, x_proj, y_proj, 'ldd')
    LD_for_mask = LD_for_mask['ldd']
    LD_for_mask = LD_for_mask.fillna(0)  # Fill NaN with 0
    LD_for_mask = LD_for_mask.where((LD_for_mask > 0) & (LD_for_mask < 10))  # Valid LDD values are 1-9

    if mask_file:
        try:
            MX = xr.open_dataset(mask_file)
            MX = prepare_dataset(MX, x_proj, y_proj, 'domain')
            MX = MX['domain']  # Extract as Dataarray
            # Combine with LDD mask
            MX = MX * LD_for_mask
        except:
            print(f'The given mask path {mask_file} is not a valid path. Using LDD domain from {ldd_file_to_use}.')
            MX = LD_for_mask
    else:
        print(f'No mask file provided. Using LDD domain from {ldd_file_to_use}.')
        MX = LD_for_mask

    # use the exact same coords from channel slope file, just in case there are precision differences
    MX = MX.assign_coords({x_proj: x_all, y_proj: y_all})

    # ---------------- Loop on the basin pixels to find how many MCT pixels they have downstream
    # initiate a counter with 1 in cells that fit the slope criteria and 0 elsewhere
    sum_rivers = rivers_mask.values.copy()

    # set the initial value of the 'downstream' pixels
    downstream_cells = rivers_mask.values.copy()

    # Calculates the maximum distance to all points from the river network sinks
    downstreamdist = distance.to_sink(network, path='shortest')
    # Create mask for valid downstream cells (not pits)
    # A pit has ldd value 5, downstreamdist returns 0 for pits
    downstream_actual_mask = (downstreamdist > 0).astype(int)

    # The loop is used to count how many pixels are MCT downstream, as at each loop we move the values 1 pixel upstream
    # At the end of the loop, each element of the array has the number of downstream MCT pixels for that pixel
    for loops in range(nloops):
        # get the value on the downstream cell and put it in a mask
        downstream_cells = move.upstream(network, downstream_cells, return_type='gridded')
        downstream_cells = downstream_cells * downstream_actual_mask
        sum_rivers = sum_rivers + downstream_cells

    # ---------------- Generate a new MCT rivers mask
    # Pixels with nloops downstream MCT pixels plus their self (nloops+1 in total) go to the MCT river mask
    mct_mask_np = (sum_rivers == nloops + 1).astype(int)
    
    # Keep only the cells over the minimum area
    mct_mask_np = mct_mask_np * minarea_bool.values
    
    # Use path function to include in the MCT mask all pixels downstream of an MCT pixel
    # path requires boolean 0-1. If there are NaNs then it gives wrong results!
    # mct_mask_np = path(ldd_array, mct_mask_np)
    mct_mask_np = upstream.max(network, mct_mask_np, return_type='gridded').values.astype(int)

    # ---------------- Generate the output file
    MCT = xr.DataArray(
        mct_mask_np,
        dims=[y_proj, x_proj]
        # coords={y_proj: y_all, x_proj: x_all}
    )
    MCT.name = 'mct_mask'

    # Copy coordinates and projection attributes from the source dataset to preserve CRS information
    MCT = copy_coordinates_and_attributes(source=LD_for_mask_ds, target=MCT, y_proj=y_proj, x_proj=x_proj)

    # mask final data with the mask_file
    MCT = MCT.where(MX.notnull())
    MX.close()
    LD_for_mask_ds.close()

    return MCT


def extract_coords(ds: xr.Dataset, coords_names: List[str] = []) -> Tuple[str, str]:
    x_proj, y_proj = '', ''
    # Determine coordinate names: use provided names or auto-detect from dataset
    if len(coords_names) == 0:
        # Auto-detect coordinates by matching against known patterns
        available_coords = set(ds.coords)
        x_matches = available_coords.intersection(X_COORD_NAMES)
        y_matches = available_coords.intersection(Y_COORD_NAMES)
        
        if len(x_matches) != 1 or len(y_matches) != 1:
            raise ValueError(
                f"Unable to uniquely identify coordinates. "
                f"Found {len(x_matches)} x-coords: {x_matches}, "
                f"{len(y_matches)} y-coords: {y_matches}. "
                f"Available: {list(available_coords)}. "
                f"Please specify coordinates explicitly using the -E argument "
                f"(format: 'lat lon')."
            )
        
        x_proj = str(next(iter(x_matches)))
        y_proj = str(next(iter(y_matches)))
    elif len(coords_names) == 2:
        # Use provided coordinates (y, x order as per documentation)
        y_proj, x_proj = coords_names
    else:
        raise ValueError(
            f"coords_names must have 0 or 2 elements, got {len(coords_names)}: {coords_names}"
        )
        
    return x_proj, y_proj


def prepare_dataset(ds: xr.Dataset, x_proj: str, y_proj: str, new_name: str) -> xr.Dataset:
    """
    Change the name of the variable in the dataset that has dimensions matching x_proj and y_proj to new_name.
    Make sure to transpose the variable to have the correct order of dimensions (y_proj, x_proj) for consistency with downstream processing.
    This is needed for the pcraster conversion.
    """
    # Use set comparison instead of sorted() to avoid type issues with Hashable
    # and to handle dimension order differences
    old_name = [i for i in list(ds.data_vars) if set(ds[i].dims) == {x_proj, y_proj}]
    if not old_name:
        raise ValueError("No variable found with dimensions matching x and y coordinates")
    if len(old_name) > 1:
        import warnings
        warnings.warn(f"Multiple variables match the dimension check: {old_name}. Using first one: {old_name[0]}")
    ds = ds.rename({old_name[0]: new_name}) # only 1 variable complies with the above check
    ds[new_name] = ds[new_name].transpose(y_proj, x_proj)  # Ensure correct order of dimensions
    return ds


def main():
    'function for running from command line'
    # ---------------- Read settingss
    args = getarg()
    channels_slope_file_arg = args.changradfile
    ldd_file_arg = args.LDDfile
    uparea_file_arg = args.uparea
    mask_file_arg = args.maskfile
    slp_threshold_arg = args.slope
    nloops_arg = args.nloops
    minuparea_arg = args.minuparea
    coords_names_arg = args.coordsnames
    outputfile_arg = args.outputfilename

    mct_final = mct_mask(channels_slope_file=channels_slope_file_arg, ldd_file=ldd_file_arg, uparea_file=uparea_file_arg, mask_file=mask_file_arg, 
                         slp_threshold=slp_threshold_arg, nloops=nloops_arg, minuparea=minuparea_arg, coords_names=coords_names_arg)
    # lisflood does not read NaNs so the data are saved as boolean 0-1, with 100 being flagged as NaN for python reading
    mct_final.to_netcdf(outputfile_arg, encoding={"mct_mask": {'_FillValue': 100, 'dtype': 'int8'}})
    
    
if __name__ == "__main__":
    main()