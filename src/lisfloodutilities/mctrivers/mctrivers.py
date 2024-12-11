#   This is a script to create a mask with mild sloping river pixels.
#   It takes LISFLOOD channels slope map (changrad.nc), the LDD (ldd.nc) and mask of the catchment/domain.
#   Pixels where river slope < threshold are added to the mask, if drainage area is large enough.
#   A minimum number of consecutive mild sloping downstream pixels is required for a pixel to be added to the mask.
#   It requires PCRaster.
#   
#   Usage:
#   mctrivers.py -i changrad.nc -l ec_ldd.nc -m mask.nc -u upArea.nc -E y x -S 0.001 -N 5 -U 500 -O chanmct.nc


import xarray as xr
import pcraster as pcr
import numpy as np


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
    parser.add_argument('-E', '--coordsnames', type=str, nargs='+', required=False, default="None",
                        help='Coords names for lat, lon (in this order with space!) from the netcdf files used')
    parser.add_argument('-O', '--outputfilename', type=str, required=False, default="chanmct.nc",
                        help='Output file (chanmct.nc)')
    args = parser.parse_args()  # assign namespace to args
    return args

    
def mct_mask(channels_slope_file, ldd_file, uparea_file, mask_file='', 
             slp_threshold=0.001, nloops=5, minuparea=0, coords_names='None', 
             outputfile='chanmct.nc'):
    """
    
    Builds a mask of mild sloping rivers for use in LISFLOOD with MCT diffusive river routing. It takes LISFLOOD channels slope map (changrad.nc), the LDD (ldd.nc), 
    the upstream drained area map (upArea.nc) and the catchment/domain mask (mask.nc), and outputs a bolean mask (chanmct.nc). Pixels where riverbed gradient < threshold 
    (slp_threshold) are added to the mask if their drainage area is large enough (minuparea) and they also have at least nloops consecutive downstream pixels that meet
    the same condition for slope (drainage area will be met as downstream the area increases).

    Usage:
    The tool requires the following input arguments:
    
    channels_slope_file: LISFLOOD channels gradient map (changrad.nc)
    ldd_file: LISFLOOD local drain direction file (ldd.nc)
    uparea_file: LISFLOOD Uustream area file (upArea.nc)
    mask_file: a mask nc file; if not given (default) all cells are considered valid.
    slp_threshold: Riverbed slope threshold to use MCT diffusive wave routing (default: 0.001)
    nloops: Number of consecutive downstream grid cells that also need to comply with the slope requirement for including a grid cell in the MCT rivers mask (default: 5)
    minuparea: Minimum upstream drainage area for a pixel to be included in the MCT rivers mask (uses the same units as in the -u file) (default: 0)
    coords_names: Coordinates names for lat, lon (in this order as list) used in the the netcdf files (default: 'None': checks for commonly used names)
    outputfile: Output file containing the rivers mask where LISFLOOD can use the MCT diffusive wave routing (default: chanmct.nc)
    
    Example for generating an MCT rivers mask with pixels where riverbed slope < 0.001, drainage area > 500 kms and at least 5 downstream pixels meet the same 
    two conditions, considering the units of the upArea.nc file are given in kms:
    
    mct_mask(channels_slope_file='changrad.nc', ldd_file='ldd.nc', uparea_file='upArea.nc', mask_file='mask.nc', 
             slp_threshold=0.001, nloops=5, minuparea=0, coords_names=['y' , 'x'], 
             outputfile='chanmct.nc')
    """
    # ---------------- Read LDD (Note that for EFAS5 there is small shift of values for CH)
    LD = xr.open_dataset(ldd_file)
    
    # ---------------- Auxiliary variables
    x_checks = ['lon', 'x', 'rlon']
    y_checks = ['lat', 'y', 'rlat']
    if coords_names[0] == "None":
        x_proj = set(list(LD.coords)) & set(x_checks)
        y_proj = set(list(LD.coords)) & set(y_checks)
    
        if len(x_proj)!=1 or len(y_proj)!=1:
            print('Input dataset coords names for lat/lon are not as expected.')
            print(f'The available coords are: {list(LD.coords)}')
            print(f'The checked names are {y_checks} and {x_checks} for lat and lon respectively.')
            print('Please use -E argument and give the coords names !with space in between! in order: lan lon')
            exit(1)
        
        x_proj = list(x_proj)[0]
        y_proj = list(y_proj)[0]
    else:
        y_proj, x_proj = coords_names
    
    # ---------------- Process channels slope netcdf
    # proprocess CH dataset for having correct format
    CH = xr.open_dataset(channels_slope_file)
    old_name = [i for i in list(CH.data_vars) if sorted(CH[i].dims)==sorted([x_proj, y_proj])]
    CH = CH.rename({old_name[0]: "changrad"})  # only 1 variable complies with above check
    CH['changrad'] = CH['changrad'].transpose(y_proj, x_proj)  # make sure dims order is as pcraster needs

    # ---------------- Set clone map for pcraster
    # get number of rows and columns
    rows, cols = CH.sizes[y_proj], CH.sizes[x_proj]

    # get coords of map corners
    x_all = CH.variables[x_proj]
    y_all = CH.variables[y_proj]
    x1 = x_all[0]
    x2 = x_all[-1]
    y1 = y_all[0]

    # calc cell size
    cell_size = np.abs(x2 - x1) / (cols - 1)
    # calc coords
    x = x1 - cell_size / 2
    y = y1 + cell_size / 2

    # set the clone map for pcraster
    pcr.setclone(rows, cols, cell_size, x, y)

    # ---------------- Create rivers mask
    rivers_mask = (CH.changrad < slp_threshold)*1
    CH.close()

    # convert the xarray to pcraster
    rivers_mask_pcr = pcr.numpy2pcr(pcr.Scalar, rivers_mask.values, 0) 

    # ---------------- Process LDD
    old_name = [i for i in list(LD.data_vars) if sorted(LD[i].dims)==sorted([x_proj, y_proj])]
    LD = LD.rename({old_name[0]: "ldd"})['ldd']  # only 1 variable complies with above if

    # sometimes the masked value is flagged not with NaN (e.g., with cutmaps it was flagged with 0)
    # pcr.Ldd takes only integer values 1-9, so any other value needs to be masked
    LD = LD.fillna(-1)  # fill NaN so it can be converted to integer with no issues
    LD = LD.astype('int')
    LD = LD.where((LD>0) & (LD<10)).fillna(-1)
    
    # convert the xarray to pcraster
    LD = LD.transpose(y_proj, x_proj)  # make sure dims order is as pcraster needs
    ldd_pcr = pcr.numpy2pcr(pcr.Ldd, LD.values, -1)  # missing values in the np are flagged as -1

    # repair the ldd; needed in case ldd is created from cutmaps, so outlet is not flagged with 5 (pit) 
    ldd_pcr = pcr.lddrepair(ldd_pcr)
    

    # ---------------- Read upstream area
    UA = xr.open_dataset(uparea_file)
    old_name = [i for i in list(UA.data_vars) if sorted(UA[i].dims)==sorted([x_proj, y_proj])]
    UA = UA.rename({old_name[0]: "domain"})['domain']  # only 1 variable complies with above if
    
    # convert the xarray to pcraster
    UA = UA>=minuparea # check that the area is over the minimum
    minarea_bool_pcr = pcr.numpy2pcr(pcr.Boolean, UA.values, 0)
    UA.close()
    
    # ---------------- Read domain/basin (mask) area
    try:
        MX = xr.open_dataset(mask_file)
        old_name = [i for i in list(MX.data_vars) if sorted(MX[i].dims)==sorted([x_proj, y_proj])]
        MX = MX.rename({old_name[0]: "domain"})['domain']  # only 1 variable complies with above if
    except:
        print(f'The given mask path {mask_file} is not a valid path. All domain read from LDD file {ldd_file} is considered vaid.')
        MX = LD.copy(deep=True)
        MX = MX.fillna(0)*0+1

    # use the exact same coords from channel slope file, just in case there are precision differences
    MX = MX.assign_coords(x_proj=x_all, y_proj=y_all)
    LD.close()  # close ther LD file, after the check of mask availability
    
    # ---------------- Loop on the basin pixels to find how many MCT pixels they have downstream
    # initiate a counter with 1 in cells that fit the slope criteria and 0 elsewhere
    sum_rivers_pcr = rivers_mask_pcr
    
    # set the initial value of the 'downstream' pixels
    downstream_cells_pcr = rivers_mask_pcr
    
    # Loop nloops times and use downstream function to find out if each cell has nloop MCT cells downstream
    # Downstream function gives the value in the downstream pixel in a map:
    # here it gives 1 if the downstream pixel is Muskingum, zero otherwise.
    # The loop is used to count how many pixels are MCT downstream, as at each loop we move the values 1 pixel upstream
    # At the end of the loop, each element of the array has the number of downstream MCT pixels for that pixel
    for loops in range(0, nloops):
        # get the value on the downstream cell and put it in a mask
        downstream_cells_pcr = pcr.downstream(ldd_pcr, downstream_cells_pcr)
        sum_rivers_pcr = sum_rivers_pcr + downstream_cells_pcr
        
    # ---------------- Generate a new MCT rivers mask
    # Pixels with nloops downstream MCT pixels plus their self (nloops+1 in total) go to the MCT river mask
    mct_mask_pcr = pcr.ifthenelse(sum_rivers_pcr == nloops+1, pcr.boolean(1), pcr.boolean(0))
    
    # Keep only the cells over the minimum area
    mct_mask_pcr = pcr.ifthenelse(minarea_bool_pcr, mct_mask_pcr, pcr.boolean(0))
    
    # Use path function to include in the MCT mask all pixels downstream of an MCT pixel
    # path requires boolean 0-1. If there are NaNs then it gives wrong results!
    mct_mask_pcr = pcr.pcr2numpy(mct_mask_pcr, 0)  # get the numpy from pcr
    mct_mask_pcr = pcr.numpy2pcr(pcr.Boolean, mct_mask_pcr, -1)  # convert to Boolean with no NaN (-1 is not possible)
    mct_mask_pcr = pcr.path(ldd_pcr, mct_mask_pcr)  # get the actual paths
    
    # ---------------- Generate the output file
    mct_mask_np = pcr.pcr2numpy(mct_mask_pcr, 0)
    MCT = MX.fillna(0)*0+mct_mask_np
    MX.close()
    MCT.name = 'mct_mask'
    
    # mask final data with the mask_file
    MCT = MCT.where(MX==1)
    
    # lisflood does not read NaNs so the data are saved as boolean 0-1, with 0 being flagged as NaN for python reading
    MCT.to_netcdf(outputfile, encoding={"mct_mask": {'_FillValue': 0, 'dtype': 'int8'}})
    return MCT


def main():
    'function for runnign from command line'
    # ---------------- Read settings
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

    mct_mask(channels_slope_file=channels_slope_file_arg, ldd_file=ldd_file_arg, uparea_file=uparea_file_arg, mask_file=mask_file_arg, 
             slp_threshold=slp_threshold_arg, nloops=nloops_arg, minuparea=minuparea_arg, coords_names=coords_names_arg, 
             outputfile=outputfile_arg)

    
if __name__ == "__main__":
    main()