#   This is a script to retrieve from a selected point and variable the values for all downstream points, ordered from upstream to downstream
#   It takes the netcdf of the variable of interest, the LDD (ldd.nc) and the coordinates of the point of interest
#   It requires PCRaster.
#   
#   Usage:
#   profiler.py -i changrad.nc -l ldd.nc -X 10.01 -Y 24.01 -E y x -O profiler.csv


import pcraster as pcr
import xarray as xr
import pandas as pd
import numpy as np

def getarg():
    """ Get program arguments.

    :return: args:  namespace of program arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--Inputfile', type=str, required=True,
                        help='Input nc file from where the values of the downstream areas will be extracted')
    parser.add_argument('-l', '--LDDfile', type=str, required=True,
                        help='LISFLOOD LDDD file (ldd.nc)')
    parser.add_argument('-X', '--Location_x', type=str, required=True,
                        help='Longitude of the initial point [use the same projection system as the input nc files!]')
    parser.add_argument('-Y', '--Location_y', type=str, required=True,
                        help='Latitude of the initial point [use the same projection system as the input nc files!]')
    parser.add_argument('-E', '--coordsnames', type=str, nargs='+', required=False, default="None",
                        help='Coords names for lat, lon (in this order with space!) from the netcdf files used')
    parser.add_argument('-O', '--outputfilename', type=str, required=False, default="profiler.csv",
                        help='Output file (profiler.csv)')
    args = parser.parse_args()  # assign namespace to args
    return args

    
def profiler(input_file, ldd_file, x_coord, y_coord, coords_names='None'):
    """
    
    Creates a dataframe with the values for all grid cells downstream of a user-defined point from a netcdf file specified by the user

    Usage:
    The tool requires the following input arguments:
    
    Inputfile: Input nc file from where the values of the downstream areas will be extracted (e.g., chagrad.nc)
    LDDfile: LISFLOOD local drain direction file (ldd.nc)
    Location_x: LISFLOOD Uustream area file (upArea.nc)
    Location_y: a mask nc file; if not given (default) all cells are considered valid.
    coords_names: Coordinates names for lat, lon (in this order as list) used in the the netcdf files (default: 'None'; checks for commonly used names ['x', 'lon', 'rlon'], similar for lat names)
    outputfile: Output file containing the extracted values of the downsteam points, their coordinates, and their order in the river network (default: profiler.csv)
    
    Example:
    profiler(input_file='changrad.nc', ldd_file='ldd.nc', Location_x=15.21, Location_y=25.14, coords_names=['y' , 'x'])
    """
    # ---------------- Read LDD (Note that for EFAS5 there is small shift of values for CH)
    LDD = xr.open_dataset(ldd_file)
    input_file = xr.open_dataset(input_file)
    
    # ---------------- Auxiliary variables
    x_checks = ['lon', 'x', 'rlon']
    y_checks = ['lat', 'y', 'rlat']
    if coords_names == "None":
        x_proj = set(list(LDD.coords)) & set(x_checks)
        y_proj = set(list(LDD.coords)) & set(y_checks)
    
        if len(x_proj)!=1 or len(y_proj)!=1:
            print('Input dataset coords names for lat/lon are not as expected.')
            print(f'The available coords are: {list(LDD.coords)}')
            print(f'The checked names are {y_checks} and {x_checks} for lat and lon respectively.')
            print('Please use -E argument and give the coords names !with space in between! in order: lan lon')
            exit(1)
        
        x_proj = list(x_proj)[0]
        y_proj = list(y_proj)[0]
    else:
        y_proj, x_proj = coords_names

    # assign values of coordinates in both dataset based on LDD, in case of small precision issues
    input_file = input_file.assign_coords({y_proj: LDD[y_proj].values, x_proj: LDD[x_proj].values})
    
    # ---------------- Process LDD
    old_name = [i for i in list(LDD.data_vars) if sorted(LDD[i].dims)==sorted([x_proj, y_proj])]
    LDD = LDD.rename({old_name[0]: "ldd"})['ldd']  # only 1 variable complies with above if
    
    # sometimes the masked value is flagged not with NaN (e.g., with cutmaps it was flagged with 0)
    # pcr.Ldd takes only integer values 1-9, so any other value needs to be masked
    LDD = LDD.fillna(-1)  # fill NaN so it can be converted to integer with no issues
    LDD = LDD.astype('int')
    LDD = LDD.where((LDD>0) & (LDD<10)).fillna(-1)
    
    # convert the xarray to pcraster
    LDD = LDD.transpose(y_proj, x_proj)  # make sure dims order is as pcraster needs

    # ---------------- Set clone map for pcraster
    rows, cols = len(LDD[y_proj]), len(LDD[x_proj])
    pcr.setclone(rows, cols, 1, 0, 0)

    ldd_pcr = pcr.numpy2pcr(pcr.Ldd, LDD.values, -1)  # missing values in the np are flagged as -1

    # repair the ldd; needed in case ldd is created from cutmaps, so outlet is not flagged with 5 (pit) 
    ldd_pcr = pcr.lddrepair(ldd_pcr)

    # ---------------- Get coordinates index for the point of interest    
    try:
        final_point = LDD.sel({y_proj: y_coord, x_proj: x_coord} , method='nearest', tolerance=0.1)
    except:
        print('The provided coordinates are not valid, please check again!')
        exit()
        
    Y_index = np.argmax(LDD[y_proj].values==final_point[y_proj].values)
    X_index = np.argmax(LDD[x_proj].values==final_point[x_proj].values)

    # ---------------- Create a mask with the downstream points 
    profile_mask_pcr = LDD.fillna(0)*0
    profile_mask_pcr[Y_index, X_index] = 1
    profile_mask_pcr = profile_mask_pcr.astype('int')

    profile_mask_pcr = pcr.numpy2pcr(pcr.Boolean, profile_mask_pcr.values, -1)  # convert to Boolean with no NaN (-1 is not possible)
    profile_mask_pcr = pcr.path(ldd_pcr, profile_mask_pcr)  # get the actual paths

    profile_mask_np = pcr.pcr2numpy(profile_mask_pcr, 0)
    ProfPath = LDD.fillna(0)*0+profile_mask_np
    ProfPath = ProfPath.where(ProfPath==1)
    ProfPath.name = 'Profile'
    total_points = int(ProfPath.sum().values)
    print(f'There are {total_points} points from the point of interest until the end of the river network')

    # ---------------- Subset data for keeping only the domain of interest, and speeding up the analysis 
    lats_used = ProfPath.sum(x_proj).where(ProfPath.sum(x_proj)!=0).dropna(y_proj)
    lons_used = ProfPath.sum(y_proj).where(ProfPath.sum(y_proj)!=0).dropna(x_proj)
    
    ProfPath = ProfPath.sel({y_proj:lats_used[y_proj].values, x_proj:lons_used[x_proj].values})

    input_file = input_file.sel({y_proj: ProfPath[y_proj].values, x_proj: ProfPath[x_proj].values})
    LDD = LDD.sel({y_proj: ProfPath[y_proj].values, x_proj: ProfPath[x_proj].values})
    
    # define clone for pcraster for cropped data
    rows, cols = len(ProfPath[y_proj]), len(ProfPath[x_proj])
    pcr.setclone(rows, cols, 1, 0, 0)

    rivers_mask_pcr = pcr.numpy2pcr(pcr.Scalar, ProfPath.values, 0)
    ldd_pcr = pcr.numpy2pcr(pcr.Ldd, LDD.values, -1)
    ldd_pcr = pcr.lddrepair(ldd_pcr)  # need to repair again, because the outlet is (possible) lost due to subsetting
    
    # ---------------- Find the order of the points in the river network 
    # generate initial data with the river mask
    downstream_cells_pcr = rivers_mask_pcr
    sum_rivers_pcr = rivers_mask_pcr

    # modify data, so that the most downstream point is masked out (otherwise the below loop gives for all points the same value)
    downstream_actual_mask_pcr = pcr.downstreamdist(ldd_pcr)
    downstream_actual_mask_pcr = pcr.ifthenelse(downstream_actual_mask_pcr == 0, pcr.boolean(0), pcr.boolean(1))
    downstream_actual_mask_pcr = pcr.scalar(downstream_actual_mask_pcr)
    
    # Loop {number of cells -1} times and use downstream function to find out the order of the cells on the river flow
    print('Calculating the order of the identified points')
    for loops in range(total_points-1):
        downstream_cells_pcr = pcr.downstream(ldd_pcr, downstream_cells_pcr)
        downstream_cells_pcr = downstream_cells_pcr*downstream_actual_mask_pcr
        sum_rivers_pcr = sum_rivers_pcr + downstream_cells_pcr

    sum_rivers = (ProfPath.fillna(0)*0+pcr.pcr2numpy(sum_rivers_pcr, 0))
    sum_rivers.name = 'order'

    all_data = xr.merge([input_file, LDD, ProfPath, sum_rivers])
    all_data_df = all_data.to_dataframe().reset_index()
    all_data_df = all_data_df[~all_data_df.Profile.isna()]
    all_data_df = all_data_df.sort_values('order')

    return all_data_df


def main():
    'function for running from command line'
    # ---------------- Read settings
    args = getarg()
    input_file_arg = args.Inputfile
    ldd_file_arg = args.LDDfile
    x_coord_arg = args.Location_x
    y_coord_arg = args.Location_y
    coords_names_arg = args.coordsnames
    outputfile_arg = args.outputfilename
        
    data_df = profiler(input_file=input_file_arg, ldd_file=ldd_file_arg, x_coord=x_coord_arg, y_coord=y_coord_arg, coords_names=coords_names_arg)
    data_df.to_csv(outputfile_arg)


if __name__ == "__main__":
    main()
