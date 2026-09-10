import xarray as xr

def getarg():
    """
    Get program arguments.

    :return: args:  namespace of program arguments
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="Path of input netcdf file.",
    )
    parser.add_argument(
        "-m",
        "--mask",
        type=str,
        required=False,
        default="/ec/ws4/tc/emos/work/cems/floods/glofas/assets/4.0/maps/elv_Global_03min.nc",
        help="Path of mask file. Default is the elevation mask used for ta/td.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="Path of output file.",
    )
    args = parser.parse_args()  # assign namespace to args
    return args


def main():
    """Function for running the whole script as main"""
    
    # ----------- Read arguments -----------
    args = getarg()
    infile = args.infile
    mask = args.mask
    outfile = args.outfile

    # ----------- Read data ----------- 
    indata = xr.open_mfdataset(infile)
    mask = xr.open_dataset(mask)
    # get the variable of interest for the mask, and convert it to boolean
    mask_var = [i for i in list(mask.data_vars) if len(mask[i].dims)>=2][0]
    mask = mask[mask_var]
    mask = mask.notnull()*1
    mask = mask.where(mask==1)

    # ----------- Mask data ----------- 
    outdata = indata/mask
    outdata = outdata.astype('float32')  # convert to float32, no need to be in float64

    # ----------- Save data ----------- 
    # get current encoding (remove all keys not available for netcdf encoding)
    kept_keys = ['dtype', 'zlib', 'shuffle', 'complevel', 'fletcher32', 'chunksizes', 'original_shape']
    kept_keys += ['missing_value', '_FillValue', 'scale_factor', 'add_offset']
    encoding = {i:indata[i].encoding for i in list(indata.data_vars.keys())}
    for var in list(indata.data_vars.keys()):
        for i in list(encoding[var].keys()):
            if i not in kept_keys:
                encoding[var].pop(i, None)   
    
    outdata.to_netcdf(outfile, encoding=encoding)  # save data


if __name__ == "__main__":
    main()