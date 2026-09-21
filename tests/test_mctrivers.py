import os
import shutil
import xarray as xr


from earthkit.hydro.river_network import repair as river_network_repair
from lisfloodutilities.mctrivers.mctrivers import extract_coords, mct_mask, prepare_dataset


def mk_path_out(p):
    'make folder for storing output data'
    path_out = os.path.join(os.path.dirname(__file__), p)
    if os.path.exists(path_out):
        shutil.rmtree(path_out, ignore_errors=True)
    os.mkdir(path_out)


class TestMctMask():
    '''
    Test the mct mask generation
    '''

    case_dir = os.path.join(os.path.dirname(__file__), 'data', 'mctrivers')

    def run(self, slp_threshold, nloops, minuparea, coords_names, case_name):
        'generate the MCT mask'
        # setting
        self.out_path_ref = os.path.join(self.case_dir, case_name, 'reference')
        self.out_path_run = os.path.join(self.case_dir, case_name, 'out')
        mk_path_out(self.out_path_run)

        channels_slope_file = os.path.join(self.case_dir, case_name, 'changrad.nc')
        original_ldd_file = os.path.join(self.case_dir, case_name, 'ldd.nc')
        ldd_file = os.path.join(self.out_path_run, 'ldd.nc')
        # Copy the ldd file to the out folder so the function can create the
        # repaired ldd in case it needs repair and later clean the output folder
        shutil.copyfile(original_ldd_file, ldd_file)
        uparea_file = os.path.join(self.case_dir, case_name, 'upArea.nc')
        mask_file = os.path.join(self.case_dir, case_name, 'mask.nc')
        outputfile = os.path.join(self.out_path_run, 'mctmask.nc')

        # generate the mct river mask
        mct_final = mct_mask(channels_slope_file, ldd_file, uparea_file, mask_file, slp_threshold, nloops, minuparea, coords_names)
        mct_final.to_netcdf(outputfile, encoding={"mct_mask": {'_FillValue': 0, 'dtype': 'int8'}})

        # compare the generated mask with the reference one
        ref_file = os.path.join(self.out_path_ref, 'mctmask.nc')
        reference = xr.open_dataset(ref_file)
        out_file = os.path.join(self.out_path_run, 'mctmask.nc')
        generated = xr.open_dataset(out_file)
        # check if same based on https://docs.xarray.dev/en/stable/generated/xarray.DataArray.equals.html
        all_equal = reference.equals(generated)
        generated.close()  # needs to be closed otherwise the out folder can't be deleted

        # Delete the out folder
        shutil.rmtree(self.out_path_run, ignore_errors=True)

        fail_message = f'Test for mct river mask generation for {case_name} failed. Please check differences between' 
        fail_message += f' the generated mask "{out_file}" and the expected mask "{ref_file}".'
        assert all_equal, fail_message


class TestMctrivers(TestMctMask):
    '''
    Test the mctrivers functionality
    '''

    def test_mctrivers_etrs89(self):
        '''Test mctrivers with ETRS89 coordinates'''
        self.run(0.001, 5, 500*10**6, [], 'LF_ETRS89_UseCase')

    def test_mctrivers_latlon(self):
        '''Test mctrivers with lat/lon coordinates'''
        self.run(0.001, 5, 500*10**6, ['lat', 'lon'], 'LF_lat_lon_UseCase')
    
    def test_earthkit_hydro_lddrepair(self):
        '''Test that the lddrepair function from earthkit.hydro.river_network works as expected'''
        case_name = 'LF_lat_lon_UseCase'
        ldd_path = os.path.join(self.case_dir, case_name, 'ldd.nc')
        ldd_repaired_path = os.path.join(self.case_dir, case_name, 'ldd_repaired.nc')
        output_path = os.path.join(self.case_dir, case_name, 'ldd_repaired_by_earthkit.nc')
        river_network_repair(ldd_path, output_path, river_network_format='pcr_d8', input_source='file')

        coords_names = ['lat', 'lon']
        LD_ds = xr.open_dataset(output_path)
        x_proj, y_proj = extract_coords(LD_ds, coords_names)
        LD = prepare_dataset(LD_ds, x_proj, y_proj, 'ldd')
        LD = LD['ldd']
        LD = LD.fillna(-1)
        LD = LD.astype('int')
        LD = LD.where((LD > 0) & (LD < 10)).fillna(-1)
        ldd_array = LD.values.astype(int)

        LD_REPAIRED_ds = xr.open_dataset(ldd_repaired_path)
        LD_REPAIRED = prepare_dataset(LD_REPAIRED_ds, x_proj, y_proj, 'ldd_repaired')
        LD_REPAIRED = LD_REPAIRED['ldd_repaired']
        LD_REPAIRED = LD_REPAIRED.fillna(-1)
        LD_REPAIRED = LD_REPAIRED.astype('int')
        ldd_repaired_array = LD_REPAIRED.values.astype(int)

        LD_ds.close()  # needs to be closed otherwise the out folder can't be deleted
        LD_REPAIRED_ds.close()  # needs to be closed otherwise the out folder can't be deleted
        try:
            os.remove(output_path)
        except FileNotFoundError:
            pass

        assert (ldd_array == ldd_repaired_array).all(), 'LDD repair did not produce the expected output. Please check the generated repaired LDD and compare it with the expected repaired LDD.'
