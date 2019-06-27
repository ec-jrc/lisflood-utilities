import argparse
import os
import sys
import glob

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
                             'main variable i.e. pr/e0/tx/tavg must be last variable in the input file', required=False)
    parser.add_argument("-f", "--folder", help='Directory with netCDF files to be cut', required=False)
    parser.add_argument("-o", "--outpath", help='path where to save cut files', required=True)
    args = parser.parse_args()

    filelist = args.list
    mask = args.mask
    pathout = args.outpath
    input_folder = args.folder
    print('File list{} Mask: {} Output: {} Input: {}'.format(filelist, mask, pathout, input_folder))
    if filelist:
        list_to_cut = open(filelist).readlines()
    elif input_folder:
        list_to_cut = glob.glob(os.path.join(input_folder, '*.nc'))
    else:
        raise ValueError('You must issue either list or folder input arguments.')

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
    print('{} {} {} {}'.format(x_min, x_max, y_min, y_max))
    # walk through list_to_cut
    for f in list_to_cut:
        fileout = os.path.join(pathout, os.path.basename(f))
        if os.path.isfile(fileout) and os.path.exists(fileout):
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

        var = list(nc.variables.items())[-1][0]
        if 'lat' in nc.variables:
            nc.variables[var].attrs['esri_pe_string'] ='GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.0174532925199433]]"'
        else:
            nc.variables[var].attrs['esri_pe_string'] ='PROJCS["ETRS_1989_LAEA",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]'
            nc.variables[var].attrs['proj4_params'] = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
        if 'time' in nc.variables:
            sliced_var = nc[var][:, x_min:x_max+1, y_min:y_max+1]
        else:
            sliced_var = nc[var][x_min:x_max+1, y_min:y_max+1]

        print('Creating ', fileout)
        sliced_var.to_netcdf(fileout)
        nc.close()


if __name__ == '__main__':
    sys.exit(main())
