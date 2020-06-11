"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

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
    nc, num_dims = open_dataset(f)
    var = [v for v in nc.variables if len(nc.variables[v].dims) == num_dims][0]
    logger.info('Variable: %s', var)

    if isinstance(x_min, float):
        # bounding box input from user
        sliced_var = _cut_from_coords(nc, var, x_min, x_max, y_min, y_max)
    else:
        # user provides with indices directly (not coordinates)
        sliced_var = _cut_from_indices(nc, var, x_min, x_max, y_min, y_max)

    if sliced_var is not None:
        logger.info('Creating: %s', fileout)
        if 'missing_value' in sliced_var.encoding:
            sliced_var.encoding['_FillValue'] = sliced_var.encoding['missing_value']
        sliced_var.to_netcdf(fileout)
    if 'laea' in nc.variables or 'lambert_azimuthal_equal_area' in nc.variables:
        var = nc.variables['laea'] if 'laea' in nc.variables else nc.variables['lambert_azimuthal_equal_area']
        xr.DataArray(name='laea', data=var.data, dims=var.dims, attrs=var.attrs).to_netcdf(fileout, mode='a')
    # TODO add global attrs
    """
    // global attributes:
    :history = "Created Tue Feb 03 18:06:52 2015";
    :conventions = "CF-1.6";
    :source_software = "Python netCDF3_Classic";
    :title = "Lisflood maps for European setting Januar 2015";
    :keywords = "Lisflood, Europe";
    :source = "Lisflood European maps - pb2015";
    :institition = "JRC H01";
    """
    nc_out, _ = open_dataset(fileout)
    nc_out.attrs = nc.attrs
    nc_out.attrs['conventions'] = 'CF-1.6'
    nc_out.attrs['institution'] = 'JRC E1'
    nc_out.attrs['Source_Software'] = 'lisfloodutilities cutmaps 0.12.12'
    nc_out.attrs['source_software'] = 'lisfloodutilities cutmaps 0.12.12'
    nc_out.close()
    nc_out.to_netcdf(fileout, 'a')
    nc.close()


def open_dataset(f):
    try:
        nc = xr.open_dataset(f, chunks={'time': 100}, decode_cf=False)
        num_dims = 3
    except Exception:  # file has no time component
        num_dims = 2
        nc = xr.open_dataset(f)
    return nc, num_dims


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
    if 'lats' in nc.variables:
        lats = nc.variables['lat'][:]
        lons = nc.variables['lon'][:]
    else:
        lats = nc.variables['y'][:]
        lons = nc.variables['x'][:]
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
        list_to_cut = [f for f in Path(glofas_folder).glob('**/*') if '/.git/' not in f.as_posix()]
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
