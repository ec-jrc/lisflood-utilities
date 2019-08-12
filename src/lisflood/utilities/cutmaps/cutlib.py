import os
import sys
import logging
from pathlib import Path

import xarray as xr
import numpy as np
from pcraster.numpy_operations import pcr2numpy
from pcraster import pcraster
from netCDF4 import Dataset

logger = logging.getLogger()


def cutmap(f, fileout, x_max, x_min, y_max, y_min):
    try:
        nc = xr.open_dataset(f, chunks={'time': 100})
    except:  # file has no time component
        nc = xr.open_dataset(f)
    # FIXME assuming last variable is what we have to slice...
    var = list(nc.variables.items())[-1][0]
    if 'lat' in nc.variables:
        nc.variables[var].attrs[
            'esri_pe_string'] = 'GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.0174532925199433]]"'
    else:
        nc.variables[var].attrs[
            'esri_pe_string'] = 'PROJCS["ETRS_1989_LAEA",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]'
        nc.variables[var].attrs[
            'proj4_params'] = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
    if 'time' in nc.variables:
        sliced_var = nc[var][:, x_min:x_max + 1, y_min:y_max + 1]
    else:
        sliced_var = nc[var][x_min:x_max + 1, y_min:y_max + 1]
    logger.info('Creating: ', fileout)
    sliced_var.to_netcdf(fileout)
    nc.close()


def get_filelist(filelist, input_folder):
    list_to_cut = []
    if filelist:
        list_to_cut = open(filelist).readlines()
    elif input_folder:
        list_to_cut = [f for f in Path(input_folder).glob('**/*.nc')]
        # list_to_cut = glob.glob(os.path.join(input_folder, '*.nc'))
    logger.info('==================> Going to cut %d files', len(list_to_cut))
    return list_to_cut


def get_cuts(cuts, mask):
    if mask:
        maskname, ext = os.path.splitext(mask)
        if ext == '.map':
            if os.path.isfile(mask):
                maskmap = pcraster.setclone(mask)
                maskmap = pcraster.readmap(mask)
                masknp = pcr2numpy(maskmap, False)
            else:
                logger.error('Wrong pcraster input mask: %s not a file', mask)
                sys.exit(1)
        elif ext == '.nc':
            masknp = Dataset(mask, 'r')
        else:
            logger.error('Mask map format not recognized. Must be either .map or .nc. Found %s', ext)
            sys.exit(1)

        mask_filter = np.asarray(masknp).nonzero()
        x_min = np.min(mask_filter[0])
        x_max = np.max(mask_filter[0])
        y_min = np.min(mask_filter[1])
        y_max = np.max(mask_filter[1])
    elif cuts:
        x, y = cuts.split(':')
        x_min, x_max = x.split('-')
        y_min, y_max = x.split('-')
    return x_max, x_min, y_max, y_min
