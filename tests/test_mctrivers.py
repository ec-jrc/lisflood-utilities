import os
import shutil
import xarray as xr

from lisfloodutilities.compare.nc import NetCDFComparator

from lisfloodutilities.mctrivers.mctrivers import mct_mask


def mk_path_out(p):
    path_out = os.path.join(os.path.dirname(__file__), p)
    if os.path.exists(path_out):
        shutil.rmtree(path_out)
    os.mkdir(path_out)
    return path_out

class TestMctMask():

    case_dir = os.path.join(os.path.dirname(__file__), 'data', 'mctrivers')

    def run(self, slp_threshold, nloops, minuparea, coords_names, type):

        # setting
        self.out_path_ref = os.path.join(self.case_dir, type, 'reference')
        self.out_path_run = os.path.join(self.case_dir, type, 'out')
        mk_path_out(self.out_path_run)

        channels_slope_file = os.path.join(self.case_dir, type, 'changrad.nc')
        ldd_file = os.path.join(self.case_dir, type, 'ldd.nc')
        uparea_file = os.path.join(self.case_dir, type, 'upArea.nc')
        mask_file = os.path.join(self.case_dir, type, 'mask.nc')
        outputfile = os.path.join(self.out_path_run, 'chanmct.nc')

        # generate the mct river mask
        mct_mask(channels_slope_file, ldd_file, uparea_file, mask_file, slp_threshold, nloops, minuparea,outputfile)

        # compare results with reference
        nc_comparator = NetCDFComparator(mask_file, array_equal=True)
        nc_comparator.compare_dirs(self.out_path_run, self.out_path_ref)

    def teardown_method(self):
        print('Cleaning directories')
        out_path = os.path.join(self.case_dir, type, 'out')
        if os.path.exists(out_path) and os.path.isdir(out_path):
            shutil.rmtree(out_path, ignore_errors=True)


class TestMctrivers(TestMctMask):

    # slp_threshold = 0.001
    # nloops = 5
    # minuparea = 0
    # coords_names = 'None'

    def test_mctrivers_etrs89(self):
        self.run(0.001, 5, 0, 'None', 'LF_ETRS89_UseCase')

    def test_mctrivers_latlon(self):
        self.run( 0.001, 5, 0, 'x' 'y', 'LF_lat_lon_UseCase')

    def cleaning(self,):
        self.teardown_method()