import os

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities import IS_PYTHON2
from lisfloodutilities.cutmaps.cutlib import get_filelist, get_cuts, cutmap

if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path


class TestCutlib:

    def test_getfiles_to_cut_filelist(self):
        res = sorted(get_filelist(filelist='tests/data/cutlist.txt'))
        assert sorted(list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc',
                                      'tests/data/folder_b/ta.nc', 'tests/data/folder_d/a.nc',
                                      'tests/data/folder_d/b.nc']))) == res

    def test_getfiles_to_cut_folder(self):
        res = sorted(get_filelist(input_folder='tests/data/folder_a'))
        assert sorted(list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc']))) == res

    def test_getfiles_to_cut_glofas_setup(self):
        res = sorted(get_filelist(glofas_folder='tests/data/folder_a'))
        assert sorted(list(
            map(Path, ['tests/data/folder_a/1.map', 'tests/data/folder_a/2.map',
                       'tests/data/folder_a/3.map', 'tests/data/folder_a/4.map',
                       'tests/data/folder_a/5.map', 'tests/data/folder_a/qLakeIn.tss',
                       'tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_a/ta.nc',
                       'tests/data/folder_a/tp.nc']))) == res

    def test_get_cuts_withcuts(self):
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
            assert np.round(np.min(lons), 2) == -126.95
            assert np.round(np.min(lats), 2) == 53.25
        os.unlink(fout)

    def test_get_cuts_indices(self):
        # minxi_maxxi:minyi_maxyi
        cuts = '2_5:0_3'
        x_min, x_max, y_min, y_max = get_cuts(cuts=cuts)
        assert (x_min, x_max, y_min, y_max) == (2, 5, 0, 3)
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
        assert res_y == 53.15

    def test_get_cuts_withmaskfile(self):
        pass
