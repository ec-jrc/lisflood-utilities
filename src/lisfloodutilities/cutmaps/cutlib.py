"""

Copyright 2019-2020 European Union

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
import datetime
from pathlib import Path

import xarray as xr
import numpy as np

from dask.diagnostics import ProgressBar

from .helpers import pcraster_command
from ..readers.pcr import PCRasterMap
from ..pcr2nc import convert
from .. import version, logger

encoding_netcdf_vars = {'zlib': False}


def cutmap(f, fileout, x_min, x_max, y_min, y_max):
    nc, num_dims = open_dataset(f)
    var = [v for v in nc.variables if len(nc.variables[v].dims) == num_dims][0]
    logger.info('Variable: %s', var)

    if isinstance(x_min, float):
        # bounding box input from user  # FIXME weak isinstance test
        sliced_var = cut_from_coords(nc, var, x_min, x_max, y_min, y_max)
    else:
        # user provides with indices directly (not coordinates)
        sliced_var = cut_from_indices(nc, var, x_min, x_max, y_min, y_max)

    if sliced_var is not None:
        if 'missing_value' in sliced_var.encoding:
            sliced_var.encoding['_FillValue'] = sliced_var.encoding['missing_value']
        logger.info('Creating: %s', fileout)
        delayed_obj = sliced_var.to_netcdf(fileout, compute=False, encoding={var: encoding_netcdf_vars})
        with ProgressBar(dt=0.1):
            _ = delayed_obj.compute()

    grid_mapping = sliced_var.attrs.get('grid_mapping')
    if grid_mapping in nc.variables:
        varname = grid_mapping
        varproj = nc.variables[varname]
        logger.info('Writing projection variable: %s - %s', varname, varproj.attrs)
        del_res = xr.DataArray(name=varname, data=varproj.data, dims=varproj.dims, attrs=varproj.attrs).to_netcdf(fileout, mode='a', compute=False)
        with ProgressBar(dt=0.1):
            _ = del_res.compute()

    # adding global attributes
    nc_out, _ = open_dataset(fileout)
    nc_out.attrs = nc.attrs
    nc_out.attrs['history'] = 'lisfloodutilities cutmaps {} {} \n {}'.format(
        version, datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), nc_out.attrs.get('history', '')
    )
    nc_out.attrs['conventions'] = 'CF-1.6'
    nc_out.attrs['institution'] = 'JRC E1'
    nc_out.attrs['source_software'] = 'lisfloodutilities cutmaps {}'.format(version)
    nc_out.attrs.pop('Source_Software', None)
    nc_out.attrs.pop('Institution', None)
    nc_out.attrs.pop('Conventions', None)
    nc_out.close()
    try:
        logger.info('Writing additional attrs to: %s - %s', fileout, nc_out.attrs)
        del_res = nc_out.to_netcdf(fileout, 'a', compute=False)
        with ProgressBar(dt=0.1):
            _ = del_res.compute()
    except ValueError as e:
        logger.warning('Cannot add global attributes to %s - %s', fileout, e)
    finally:
        nc.close()


def open_dataset(f):
    try:
        nc = xr.open_dataset(f, chunks={'time': 250}, decode_cf=False)
        num_dims = 3
    except Exception:  # file has no time component
        num_dims = 2
        nc = xr.open_dataset(f, decode_cf=False)
    return nc, num_dims


def cut_from_indices(nc, var, x_min, x_max, y_min, y_max):
    # note: netcdf has lats on first dimension e.g. y_min:y_max are Y/lat dimension indices
    # that in nc file are stored on first dimension: ta(time, lat, lon)
    # you can always adjust indices in input in order to match your nc files structure
    if 'time' in nc.variables:
        sliced_var = nc[var][:, y_min:y_max + 1, x_min:x_max + 1]
    else:
        sliced_var = nc[var][y_min:y_max + 1, x_min:x_max + 1]
    return sliced_var


def cut_from_coords(nc, var, x_min, x_max, y_min, y_max):
    # we have coordinates bounds and not indices yet
    if 'lat' in nc.variables:
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


def get_filelist(input_folder=None, static_data_folder=None):
    list_to_cut = []
    if input_folder:
        list_to_cut = [f for f in Path(input_folder).glob('**/*.nc')]
    elif static_data_folder:
        list_to_cut = [f for f in Path(static_data_folder).glob('**/*') if '/.git/' not in f.as_posix()]
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


def mask_from_ldd(ldd_map, stations):
    """

    """
    try:
        from pcraster import accuflux
    except ImportError as e:
        logger.error('PCRaster not installed. Try to install PCRaster >= 4.3.0 in a conda environment '
                     'and then execute `pip install lisflood-utilities`')
        raise e

    path = os.path.dirname(stations)

    station_map = os.path.join(path, 'outlets.map')
    pcraster_command(cmd='col2map F0 F1 -N --clone F2 --large', files=dict(F0=stations, F1=station_map, F2=ldd_map))
    pcr2nc_metadata = {'variable': {'description': 'outlets points', 'longname': 'outlets', 'units': '',
                                    'shortname': 'outlets', 'mv': '0'},
                       'source': 'JRC E.1 Space, Security, Migration',
                       'reference': 'JRC E.1 Space, Security, Migration',
                       'geographical': {'datum': ''}}
    outlets_nc = os.path.join(path, 'outlets.nc')
    convert(station_map, outlets_nc, pcr2nc_metadata)

    tmp_txt = os.path.join(path, 'tmp.txt')
    tmp_map = os.path.join(path, 'tmp.map')
    with open(stations) as f:
        for line in f.readlines():
            x, y, idx = line.split()
            with open(tmp_txt, "w") as f1:
                f1.write("%s %s %s\n" % (x, y, 1))
            catchment_map = os.path.join(path, 'catchmask%05d.map' % int(idx))
            pcraster_command(cmd='col2map F0 F1 -N --clone F2 --large', files=dict(F0=tmp_txt, F1=tmp_map, F2=ldd_map))
            pcraster_command(cmd="pcrcalc 'F0 = boolean(catchment(F1, F2))'", files=dict(F0=catchment_map, F1=ldd_map, F2=tmp_map))
            pcraster_command(cmd="pcrcalc 'F0 = if((scalar(F0) gt (scalar(F0) * 0)) then F0)'", files=dict(F0=catchment_map))
    os.unlink(tmp_txt)
    os.unlink(tmp_map)

    regions_map = os.path.join(path, 'area_mask_regions.map')
    smallmask_map = os.path.join(path, 'mask.map')
    tempmask_map = os.path.join(path, 'tempmask.map')
    # init area map
    pcraster_command(cmd="pcrcalc 'F0 = scalar(F1) * 0 - 1'", files=dict(F0=regions_map, F1=ldd_map))
    with open(stations) as f:
        for line in f.readlines():
            x, y, idx = line.split()
            catchment_map = os.path.join(path, "catchmask%05d.map" % int(idx))

            pcraster_command(cmd="pcrcalc 'F0 = F0 * (1-scalar(cover(F1,0)))'",
                             files=dict(F0=regions_map, F1=catchment_map))
            pcraster_command(cmd="pcrcalc 'F0 = F0 + scalar(cover(F1, 0)) * %d'" % int(idx),
                             files=dict(F0=regions_map, F1=catchment_map))
            os.unlink(catchment_map)
    pcraster_command(cmd="pcrcalc 'F0 = boolean(if(scalar(F1) != -1, scalar(1)))'",
                     files=dict(F0=tempmask_map, F1=regions_map))
    pcraster_command(cmd='resample -c 0 F0 F1', files=dict(F0=tempmask_map, F1=smallmask_map))
    os.unlink(tempmask_map)
    os.unlink(regions_map)
    # create mask map in netCDF
    pcr2nc_metadata = {'variable': {'description': 'Mask Area', 'longname': 'mask', 'units': '',
                                    'shortname': 'mask', 'mv': '0'},
                       'source': 'JRC E.1 Space, Security, Migration',
                       'reference': 'JRC E.1 Space, Security, Migration',
                       'geographical': {'datum': ''}}
    maskmap_nc = os.path.join(path, 'my_mask.nc')
    convert(smallmask_map, maskmap_nc, pcr2nc_metadata)
    return smallmask_map, outlets_nc
