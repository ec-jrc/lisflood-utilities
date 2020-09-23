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

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities import IS_PYTHON2
from lisfloodutilities.cutmaps.cutlib import get_filelist, get_cuts, cutmap, mask_from_ldd
from lisfloodutilities.nc2pcr import convert

if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path


class TestCutlib:

    def test_getfiles_to_cut_folder(self):
        res = sorted(get_filelist(input_folder='tests/data/folder_a'))
        assert sorted(list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc']))) == res

    def test_getfiles_to_cut_glofas_setup(self):
        res = sorted(get_filelist(glofas_folder='tests/data/folder_a'))
        assert sorted(list(
            map(Path, ['tests/data/folder_a/1.map', 'tests/data/folder_a/2.map',
                       'tests/data/folder_a/3.map', 'tests/data/folder_a/4.map',
                       'tests/data/folder_a/5.map', 'tests/data/folder_a/qLakeIn.tss',
                       'tests/data/folder_a/test_tol_fail_1.tss', 'tests/data/folder_a/test_tol_ok_1.tss',
                       'tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_a/ta.nc',
                       'tests/data/folder_a/tp.nc']))) == res

    def test_get_cuts_withcoords(self):
        # lonmin_lonmax:latmin_latmax
        cuts = '-127.0_-126.5:53.2_53.4'
        x_min, x_max, y_min, y_max = get_cuts(cuts=cuts)
        assert (x_min, x_max, y_min, y_max) == (-127.0, -126.5, 53.2, 53.4)
        fin = 'tests/data/folder_a/ta.nc'
        fout = 'tests/data/folder_a/ta_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x = np.round(np.min(lons), 2)
            res_y = np.round(np.min(lats), 2)
        os.unlink(fout)
        assert res_x == -126.95
        assert res_y == 53.25

    def test_get_cuts_indices(self):
        # minxi_maxxi:minyi_maxyi
        cuts = '2_5:0_2'
        x_min, x_max, y_min, y_max = get_cuts(cuts=cuts)
        assert (x_min, x_max, y_min, y_max) == (2, 5, 0, 2)
        fin = 'tests/data/folder_a/ta.nc'
        fout = 'tests/data/folder_a/ta_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x = np.round(np.min(lons), 2)
            res_y = np.round(np.min(lats), 2)
        os.unlink(fout)
        assert res_x == -127.05
        assert res_y == 53.25

    def test_get_cuts_withmaskfile(self):
        maskfile = 'tests/data/area.nc'
        x_min, x_max, y_min, y_max = get_cuts(mask=maskfile)
        x_minr, x_maxr, y_minr, y_maxr = np.round(x_min, 2), np.round(x_max, 2), np.round(y_min, 2), np.round(y_max, 2)
        assert (x_minr, x_maxr, y_minr, y_maxr) == (-127.25, -126.15, 53.05, 53.45)
        fin = 'tests/data/area_global.nc'
        fout = 'tests/data/area_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x = np.round(np.min(lons), 2)
            res_y = np.round(np.min(lats), 2)
        os.unlink(fout)
        assert res_x == -127.25
        assert res_y == 53.05

    def test_get_cuts_withmaskpcr(self):
        maskfile = 'tests/data/asia.map'
        x_min, x_max, y_min, y_max = get_cuts(mask=maskfile)
        x_minr, x_maxr, y_minr, y_maxr = np.round(x_min, 2), np.round(x_max, 2), np.round(y_min, 2), np.round(y_max, 2)
        assert (x_minr, x_maxr, y_minr, y_maxr) == (58.65, 179.95, 0.65, 81.25)
        fin = 'tests/data/area_global.nc'
        fout = 'tests/data/area_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x = np.round(np.min(lons), 2)
            res_y = np.round(np.min(lats), 2)
        os.unlink(fout)
        assert res_x == 58.75
        assert res_y == 0.65

    def test_get_cuts_ldd(self):
        ldd = 'tests/data/cutmaps/ldd_eu.nc'
        clonemap = 'tests/data/cutmaps/area_eu.map'
        stations = 'tests/data/cutmaps/stations.txt'

        ldd_pcr = convert(ldd, clonemap, 'tests/data/cutmaps/ldd_eu_test.map', is_ldd=True)[0]
        mask = mask_from_ldd(ldd_pcr, stations)
        x_min, x_max, y_min, y_max = get_cuts(mask=mask)
        # 4067500.0 4397500.0 1282500.0 1577500.0
        assert (x_min, x_max, y_min, y_max) == (4067500, 4397500, 1282500, 1577500)

        fin = 'tests/data/cutmaps/ldd_eu.nc'
        fout = 'tests/data/area_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['x'][:]
            lats = nc.variables['y'][:]
            res_x = np.min(lons)
            res_y = np.min(lats)
        os.unlink(fout)
        os.unlink(ldd_pcr)
        assert res_x == 4067500
        assert res_y == 1282500
