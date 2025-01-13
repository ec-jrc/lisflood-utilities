import os
import shutil
import xarray as xr

from lisfloodutilities.mctrivers.mctrivers import mct_mask


def mk_path_out(p):
    'make folder for storing output data'
    path_out = os.path.join(os.path.dirname(__file__), p)
    if os.path.exists(path_out):
        shutil.rmtree(path_out, ignore_errors=True)
    os.mkdir(path_out)


class TestMctMask():

    case_dir = os.path.join(os.path.dirname(__file__), 'data', 'mctrivers')

    def run(self, slp_threshold, nloops, minuparea, coords_names, case_name):
        'generate the MCT mask'
        # setting
        self.out_path_ref = os.path.join(self.case_dir, case_name, 'reference')
        self.out_path_run = os.path.join(self.case_dir, case_name, 'out')
        mk_path_out(self.out_path_run)

        channels_slope_file = os.path.join(self.case_dir, case_name, 'changrad.nc')
        ldd_file = os.path.join(self.case_dir, case_name, 'ldd.nc')
        uparea_file = os.path.join(self.case_dir, case_name, 'upArea.nc')
        mask_file = os.path.join(self.case_dir, case_name, 'mask.nc')
        outputfile = os.path.join(self.out_path_run, 'mctmask.nc')

        # generate the mct river mask
        mct_mask(channels_slope_file, ldd_file, uparea_file, mask_file, slp_threshold, nloops, minuparea, coords_names, outputfile)

        # compare the generated mask with the reference one
        ref_file = self.out_path_ref+'/mctmask.nc'
        reference = xr.open_dataset(ref_file)
        out_file = self.out_path_run+'/mctmask.nc'
        generated = xr.open_dataset(out_file)
        # check if same based on https://docs.xarray.dev/en/stable/generated/xarray.DataArray.equals.html
        all_equal = reference.equals(generated)
        generated.close()  # needs to be closed otherwise the out folder can't be deleted
        
        if not all_equal:
            fail_message = f'Test for mct river mask generation for {case_name} failed. Please check differences between' 
            fail_message += f' the generated mask "{out_file}" and the expected mask "{ref_file}".'
            assert all_equal, fail_message
        else:
            # if equal just delete the out folder
            shutil.rmtree(self.out_path_run, ignore_errors=True)


class TestMctrivers(TestMctMask):

    # slp_threshold = 0.001
    # nloops = 5
    # minuparea = 500*10**6
    # coords_names = 'None'

    def test_mctrivers_etrs89(self):
        self.run(0.001, 5, 500*10**6, 'None', 'LF_ETRS89_UseCase')

    def test_mctrivers_latlon(self):
        self.run(0.001, 5, 500*10**6, ['lat', 'lon'], 'LF_lat_lon_UseCase')