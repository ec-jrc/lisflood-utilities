import os
import sys
import logging
from pathlib import Path

import xarray as xr
import numpy as np

from lisfloodutilities.readers import PCRasterMap

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


def cutmap(f, fileout, x_min, x_max, y_min, y_max):
    try:
        nc = xr.open_dataset(f, chunks={'time': 100}, decode_cf=False)
        num_dims = 3
    except Exception:  # file has no time component
        num_dims = 2
        nc = xr.open_dataset(f)
    var = [v for v in nc.variables if len(nc.variables[v].dims) == num_dims][0]
    logger.info('Variable: %s', var)

    if 'lat' in nc.variables:
        nc.variables[var].attrs['esri_pe_string'] = 'GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.0174532925199433]]"'
    else:
        nc.variables[var].attrs['esri_pe_string'] = 'PROJCS["ETRS_1989_LAEA",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]'
        nc.variables[var].attrs['proj4_params'] = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
    if isinstance(x_min, float):
        # bounding box input from user
        sliced_var = _cut_from_coords(nc, var, x_min, x_max, y_min, y_max)
    else:
        # user provides with indices directly (not coordinates)
        sliced_var = _cut_from_indices(nc, var, x_min, x_max, y_min, y_max)

    if sliced_var is not None:
        logger.info('Creating: %s', fileout)
        sliced_var.to_netcdf(fileout)
    nc.close()


def _cut_from_indices(nc, var, x_min, x_max, y_min, y_max):
    # note: netcdf has lats on first dimension e.g. y_min:y_max are Y/lat dimension indices
    # that in nc file are stored on first dimension: ta(time, lat, lon)
    # you can always adjust indices in input in order to match your nc files structure
    if 'time' in nc.variables:
        sliced_var = nc[var][:, y_min:y_max + 1, x_min:x_max + 1]
    else:
        sliced_var = nc[var][y_min:y_max + 1, x_min:x_max + 1]
    return sliced_var


def _cut_from_coords(nc, var, x_min, x_max, y_min, y_max):
    # we have coordinates bounds and not indices yet
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    # find indices
    ys = np.where((lats >= y_min) & (lats <= y_max))[0]
    xs = np.where((lons >= x_min) & (lons <= x_max))[0]
    if 'time' in nc.variables:
        sliced_var = nc[var][:, ys, xs]
    else:
        try:
            sliced_var = nc[var][ys, xs]
        except IndexError:
            # it happens when cutting lat or lon netcdf files (not 2D array)
            if var not in ('lats', 'lat', 'latitude', 'latitudes', 'x', 'y', 'lons', 'longitude', 'longitudes', 'lon'):
                raise ValueError('Cannot find main variable to cut')

            sliced_var = nc[var][ys] if var in ('lats', 'lat', 'latitude', 'latitudes', 'y') else nc[var][xs]
    return sliced_var


def get_filelist(filelist=None, input_folder=None, glofas_folder=None):
    list_to_cut = []
    if filelist:
        list_to_cut = open(filelist).readlines()
        list_to_cut = [Path(l.strip()) for l in list_to_cut]
    elif input_folder:
        list_to_cut = [f for f in Path(input_folder).glob('**/*.nc')]
    elif glofas_folder:
        list_to_cut = [f for f in Path(glofas_folder).glob('**/*')]
    logger.info('==================> Going to cut %d files', len(list_to_cut))
    return list_to_cut


def get_cuts(cuts=None, mask=None):
    if mask:
        if not os.path.isfile(mask):
            raise FileNotFoundError('Wrong input mask: %s not a file' % mask)
        maskname, ext = os.path.splitext(mask)
        if ext == '.map':
            mask = PCRasterMap(mask)
            lats = mask.lats
            lons = mask.lons
            mask.close()
            x_min, x_max = float(np.min(lons)), float(np.max(lons))
            y_min, y_max = float(np.min(lats)), float(np.max(lats))
        elif ext == '.nc':
            maskmap = xr.open_dataset(mask)
            latitudes = maskmap['lat'][:]
            longitudes = maskmap['lon'][:]
            y_min, y_max = float(np.min(latitudes)), float(np.max(latitudes))
            x_min, x_max = float(np.min(longitudes)), float(np.max(longitudes))
        else:
            logger.error('Mask map format not recognized. Must be either .map or .nc. Found %s', ext)
            sys.exit(1)
    elif cuts:
        # user provided coordinates bounds
        x, y = cuts.split(':')
        apply = float if '.' in x or '.' in y else int  # user can provide coords (float) or matrix indices bbox (int)
        x_min, x_max = list(map(apply, x.split('_')))
        y_min, y_max = list(map(apply, y.split('_')))
    else:
        logger.error('You must provide either cuts (in the format minlon_maxlon:minlat_maxlat) or a mask map')
        sys.exit(1)
    logger.info('CUTS: \nmin x: %s \nmax x: %s \nmin y: %s \nmax y: %s', x_min, x_max, y_min, y_max)
    return x_min, x_max, y_min, y_max
