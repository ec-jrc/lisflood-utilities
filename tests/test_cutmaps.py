"""
Copyright 2019-2026 European Union

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
import logging
import tempfile

from lisfloodutilities.cutmaps.main import get_arg_coords

logging.basicConfig(format="%(threadName)s:%(message)s")
from lisfloodutilities.cutmaps.cutlib import get_filelist, get_cuts, cutmap, mask_from_ldd
from lisfloodutilities.nc2pcr import convert as nc2pcr_convert
from lisfloodutilities.compare.nc import NetCDFComparator
from lisfloodutilities.cutmaps.main import main

from . import TestWithCleaner

from pathlib import Path


class TestCutlib(TestWithCleaner):

    def test_getfiles_to_cut_file(self):
        res = get_filelist(input_file='tests/data/folder_a/ta.nc')
        assert [Path('tests/data/folder_a/ta.nc')] == res

    def test_getfiles_to_cut_folder(self):
        res = sorted(get_filelist(input_folder='tests/data/folder_a'))
        assert sorted(list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc']))) == res

    def test_getfiles_to_cut_static_setup(self):
        res = sorted(get_filelist(static_data_folder='tests/data/folder_a'))
        assert sorted(list(
            map(Path, ['tests/data/folder_a/1.map', 'tests/data/folder_a/2.map',
                       'tests/data/folder_a/3.map', 'tests/data/folder_a/4.map',
                       'tests/data/folder_a/5.map', 'tests/data/folder_a/qLakeIn.tss',
                       'tests/data/folder_a/test_tol_fail_1.tss', 'tests/data/folder_a/test_tol_ok_1.tss',
                       'tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_a/ta.nc',
                       'tests/data/folder_a/tp.nc']))) == res

    def test_get_cuts_withcoords(self):
        # "lonmin lonmax latmin latmax"
        cuts = '-127.0 -126.5 53.2 53.4'
        cuts = get_arg_coords(cuts)
        x_min, x_max, y_min, y_max = get_cuts(cuts=cuts)
        assert (x_min, x_max, y_min, y_max) == (-127.0, -126.5, 53.2, 53.4)
        fin = 'tests/data/folder_a/ta.nc'
        fout = 'tests/data/folder_a/ta_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x_min = round(np.min(lons), 2)
            res_y_min = round(np.min(lats), 2)
            res_x_max = round(np.max(lons), 2)
            res_y_max = round(np.max(lats), 2)
        # ta.nc input file has not exact coordinates (e.g. there is no -127.0 in lons)
        assert (-126.95, -126.55, 53.25, 53.35) == (res_x_min, res_x_max, res_y_min, res_y_max)

    def test_get_cuts_indices(self):
        # "minxi maxxi minyi maxyi"
        cuts_indices = '3 7 1 2'
        cuts_indices = get_arg_coords(cuts_indices)
        ix_min, ix_max, iy_min, iy_max = get_cuts(cuts_indices=cuts_indices)
        assert (ix_min, ix_max, iy_min, iy_max) == (3, 7, 1, 2)
        fin = 'tests/data/folder_a/ta.nc'
        fout = 'tests/data/folder_a/ta_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(fin, fout, ix_min, ix_max, iy_min, iy_max, use_coords=False)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x_min = round(np.min(lons), 2)
            res_y_min = round(np.min(lats), 2)
            res_x_max = round(np.max(lons), 2)
            res_y_max = round(np.max(lats), 2)

        assert (-126.95, -126.55, 53.25, 53.35) == (res_x_min, res_x_max, res_y_min, res_y_max)

    def test_get_cuts_withmaskfile(self):
        maskfile = 'tests/data/masks/area.nc'
        x_min, x_max, y_min, y_max = get_cuts(mask=maskfile)
        x_minr, x_maxr, y_minr, y_maxr = np.round(x_min, 2), np.round(x_max, 2), np.round(y_min, 2), np.round(y_max, 2)
        assert (x_minr, x_maxr, y_minr, y_maxr) == (-127.25, -126.15, 53.05, 53.45)
        fin = 'tests/data/masks/world.nc'
        fout = 'tests/data/area_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x_min = np.min(lons)
            res_y_min = np.min(lats)
            res_x_max = np.max(lons)
            res_y_max = np.max(lats)

        assert (x_min, x_max, y_min, y_max) == (res_x_min, res_x_max, res_y_min, res_y_max)

    def test_get_cuts_withmaskfile_compare(self):
        maskfile = 'tests/data/submask/subcatchment_mask.nc'
        x_min, x_max, y_min, y_max = get_cuts(mask=maskfile)
        assert (x_min, x_max, y_min, y_max) == (4052500.0, 4232500.0, 2332500.0, 2542500.0)
        fin = 'tests/data/submask/dis.nc'
        fout = 'tests/data/submask/dis_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['x'][:]
            lats = nc.variables['y'][:]
            res_x_min = np.min(lons)
            res_y_min = np.min(lats)
            res_x_max = np.max(lons)
            res_y_max = np.max(lats)
        assert (x_min, x_max, y_min, y_max) == (res_x_min, res_x_max, res_y_min, res_y_max)
        comparator = NetCDFComparator(array_equal=True)
        comparator.compare_files(fout, 'tests/data/submask/dis_subdomain.nc')
        comparator = NetCDFComparator(mask=maskfile, array_equal=True)
        comparator.compare_files(fin, fout)


    def test_get_cuts_withmaskpcr(self):
        maskfile = 'tests/data/masks/asia.map'
        x_min, x_max, y_min, y_max = get_cuts(mask=maskfile)
        x_minr, x_maxr, y_minr, y_maxr = np.round(x_min, 3), np.round(x_max, 3), np.round(y_min, 3), np.round(y_max, 3)
        assert (x_minr, x_maxr, y_minr, y_maxr) == (58.65, 179.95, 0.65, 81.25)
        fin = 'tests/data/masks/world.nc'
        fout = 'tests/data/area_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        with Dataset(fout) as nc:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            res_x_min = np.round(np.min(lons), 3)
            res_y_min = np.round(np.min(lats), 3)
            res_x_max = np.round(np.max(lons), 3)
            res_y_max = np.round(np.max(lats), 3)
        assert (x_minr, x_maxr, y_minr, y_maxr) == (res_x_min, res_x_max, res_y_min, res_y_max)

    def test_get_cuts_ldd(self):
        ldd = 'tests/data/cutmaps/ldd_eu.nc'
        stations = 'tests/data/cutmaps/stations.txt'

        mask, outlets_points, mask_nc = mask_from_ldd(ldd, stations)
        self.cleanups.append((os.unlink, (mask,)))
        self.cleanups.append((os.unlink, (outlets_points,)))
        self.cleanups.append((os.unlink, (mask_nc,)))

        x_min, x_max, y_min, y_max = get_cuts(mask=mask)

        fin = 'tests/data/cutmaps/ldd_eu.nc'
        fout = 'tests/data/cutmaps/area_cut.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        self.cleanups.append((os.unlink, (fout,)))
        with Dataset(fout) as nc:
            lons = nc.variables['x'][:]
            lats = nc.variables['y'][:]
            res_x_min = np.min(lons)
            res_y_min = np.min(lats)
            res_x_max = np.max(lons)
            res_y_max = np.max(lats)
        assert (x_min, x_max, y_min, y_max) == (res_x_min, res_x_max, res_y_min, res_y_max)

    def test_get_cuts_ldd_3arcmin(self):
        ldd = 'tests/data/cutmaps/ldd_3arcmin.nc'
        stations = 'tests/data/cutmaps/stations_3arcmin.txt'

        mask, outlets_points, mask_nc = mask_from_ldd(ldd, stations)
        self.cleanups.append((os.unlink, (mask,)))
        self.cleanups.append((os.unlink, (outlets_points,)))
        self.cleanups.append((os.unlink, (mask_nc,)))

        x_min, x_max, y_min, y_max = get_cuts(mask=mask)
        x_minr, x_maxr, y_minr, y_maxr = np.round(x_min, 3), np.round(x_max, 3), np.round(y_min, 3), np.round(y_max, 3)

        assert_msg = 'Unexpected cuts from LDD with 3 arcmin resolution'
        assert (x_minr, x_maxr, y_minr, y_maxr) == (10.225, 12.175, 42.975, 44.275), assert_msg

        fin = 'tests/data/cutmaps/ldd_3arcmin.nc'
        fout = 'tests/data/cutmaps/area_cut_3arcmin.nc'
        cutmap(fin, fout, x_min, x_max, y_min, y_max)
        self.cleanups.append((os.unlink, (fout,)))
        with Dataset(fout) as nc:
            lons = nc.variables['x'][:]
            lats = nc.variables['y'][:]
            res_x_min = np.round(np.min(lons), 3)
            res_y_min = np.round(np.min(lats), 3)
            res_x_max = np.round(np.max(lons), 3)
            res_y_max = np.round(np.max(lats), 3)

        assert_msg = 'Unexpected cuts from cutmap with LDD with 3 arcmin resolution'
        assert (x_minr, x_maxr, y_minr, y_maxr) == (res_x_min, res_x_max, res_y_min, res_y_max), assert_msg

    def test_get_cuts_ldd_onestation(self):
        # this tests the case when LDD is in netCDF format
        ldd = 'tests/data/cutmaps/ldd_eu.nc'
        stations = 'tests/data/cutmaps/stations2.txt'
    
        mask, outlets_points, mask_nc = mask_from_ldd(ldd, stations)
        self.cleanups.append((os.unlink, (mask,)))  # produced by mask_from_ldd
        self.cleanups.append((os.unlink, (outlets_points,)))  # produced by mask_from_ldd
        self.cleanups.append((os.unlink, (mask_nc,)))  # produced by mask_from_ldd
        x_min, x_max, y_min, y_max = get_cuts(mask=mask)

        assert_msg = 'Unexpected cuts from LDD with ETRS89 projection and one station'
        assert (x_min, x_max, y_min, y_max) == (4347500.0, 4372500.0, 1282500.0, 1307500.0), assert_msg

        fout = 'tests/data/cutmaps/ldd_eu_cut.nc'
        self.cleanups.append((os.unlink, (fout,)))
        cutmap(ldd, fout, x_min, x_max, y_min, y_max)

        with Dataset(fout) as nc:
            lons = nc.variables['x'][:]
            lats = nc.variables['y'][:]
            res_x_min = np.min(lons)
            res_y_min = np.min(lats)
            res_x_max = np.max(lons)
            res_y_max = np.max(lats)

        assert_msg = 'Unexpected cuts from cutmap with LDD with ETRS89 projection and one station'
        assert (x_min, x_max, y_min, y_max) == (res_x_min, res_x_max, res_y_min, res_y_max), assert_msg


    def test_main_with_ldd_and_stations(self):
        # Paths to test data
        ldd_path = Path('tests/data/cutmaps/ldd_eu.nc')
        stations_path = Path('tests/data/cutmaps/stations.txt')
        reference_folder = Path('tests/data/cutmaps/reference')
        temp_output = Path(tempfile.gettempdir(), 'temp_output_cutmaps')
        temp_output.mkdir(exist_ok=True)
        # Input file to cut (use the same LDD file for simplicity)
        input_file = ldd_path

        # Temprorary files that will be created by the main function based on the input file location
        path_input_file = input_file.parent
        mask_file = path_input_file / 'mask_full.nc'
        outlets_points_file = path_input_file / 'outlets.nc'
        mask_nc_file = path_input_file / 'mask_small.nc'

        # Build CLI arguments
        args = [
            '-l', str(ldd_path),
            '-N', str(stations_path),
            '-F', str(input_file),
            '-o', str(temp_output),
            '-W'  # allow overwrite if file already exists
        ]

        # Run the main function
        main(args)

        # Register cleanup for temporary files created by the main function
        self.cleanups.append((os.unlink, (mask_file,)))  # produced by mask_from_ldd
        self.cleanups.append((os.unlink, (outlets_points_file,)))  # produced by mask_from_ldd
        self.cleanups.append((os.unlink, (mask_nc_file,)))  # produced by mask_from_ldd

        # Expected output file name matches the input file basename
        output_file = temp_output / input_file.name
        assert output_file.is_file(), f"Output file {output_file} not created"

        # Verify that mask files have been copied to the output directory
        small_mask_out = temp_output / 'mask_small.nc'
        full_mask_out = temp_output / 'mask_full.nc'
        outlets_out = temp_output / 'outlets.nc'
        for f in (small_mask_out, full_mask_out, outlets_out):
            assert f.is_file(), f"Expected mask file {f} not found"

        # Basic sanity check: open the output NetCDF and ensure at least one variable exists
        with Dataset(output_file) as ds:
            assert len(ds.variables) > 0, "No variables found in output NetCDF"
        
        # refrence files
        small_mask_ref = reference_folder / 'mask_small.nc'
        full_mask_ref = reference_folder / 'mask_full.nc'
        outlets_ref = reference_folder / 'outlets.nc'
        output_ref = reference_folder / input_file.name

        # Compare the output files with the reference files
        comparator = NetCDFComparator(array_equal=True, for_testing=True)
        for out, ref in ((small_mask_out, small_mask_ref), (full_mask_out, full_mask_ref), (outlets_out, outlets_ref), (output_file, output_ref)):
            comparator.compare_files(out, ref)

