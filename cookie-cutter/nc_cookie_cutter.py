import argparse
import os
import sys

import xarray as xr
import numpy as np
from pcraster.numpy_operations import pcr2numpy
from pcraster import pcraster
from netCDF4 import Dataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mask", help='mask file (cookie-cutter), .map if pcraster, .nc if netcdf', required=True)
    parser.add_argument("-l", "--list",
                        help='list of files to be cut, in a text file, one file per line, '
                             'main variable i.e. pr/e0/tx/tavg must be last variable in the input file', required=True)
    parser.add_argument("-o", "--outpath", help='path where to save cut files', required=True)
    args = parser.parse_args()

    filelist = args.list
    mask = args.mask
    pathout = args.outpath

    list_to_cut = open(filelist).readlines()
    maskname, ext = os.path.splitext(mask)

    if ext == '.map':
        if os.path.isfile(mask):
            print('maskmap', mask)
            maskmap = pcraster.setclone(mask)
            maskmap = pcraster.readmap(mask)
            masknp = pcr2numpy(maskmap, False)
        else:
            print('wrong pcraster input mask file')
            sys.exit(1)
    elif ext == '.nc':
        masknp = Dataset(mask, 'r')

    else:
        print('mask map format not recognized')
        sys.exit(1)

    mask_filter = np.where(masknp)
    x_min = np.min(mask_filter[0])
    x_max = np.max(mask_filter[0])
    y_min = np.min(mask_filter[1])
    y_max = np.max(mask_filter[1])

    # walk through list_to_cut
    for f in list_to_cut:
        fileout = os.path.join(pathout, f)
        if os.path.isfile(fileout):
            print(fileout, 'Alreay existing, are you sure you are not swiping original file with cut version? this file will not be overwrote ')
            continue
        filename, ext = os.path.splitext(f)
        if ext != '.nc':
            print('{} is not in netcdf format, skipping...'.format(f))
            continue
        try:
            nc = xr.open_dataset(f, chunks={'time': 100})
        except:  # file has no time component
            nc = xr.open_dataset(f)

        var = nc.variables.items()[-1][0]
        if 'time' in nc.variables:
            sliced_var = nc[var][:, x_min:x_max+1, y_min:y_max+1]
        else:
            sliced_var = nc[var][x_min:x_max+1, y_min:y_max+1]

        sliced_var.to_netcdf(fileout)
        nc.close()


if __name__ == '__main__':
    sys.exit(main())
