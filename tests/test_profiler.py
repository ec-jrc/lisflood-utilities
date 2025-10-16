import os
import shutil
import pandas as pd
import xarray as xr

from lisfloodutilities.profiler.profiler import profiler


def mk_path_out(p):
    'make folder for storing output data'
    path_out = os.path.join(os.path.dirname(__file__), p)
    if os.path.exists(path_out):
        shutil.rmtree(path_out, ignore_errors=True)
    os.mkdir(path_out)


class TestProfiler():

    case_dir = os.path.join(os.path.dirname(__file__), 'data', 'profiler')

    def run(self, x_coord, y_coord, coords_names, case_name):
        'generate the profiler data'
        # setting
        self.out_path_ref = os.path.join(self.case_dir, case_name, 'reference')
        self.out_path_run = os.path.join(self.case_dir, case_name, 'out')
        mk_path_out(self.out_path_run)

        input_file = os.path.join(self.case_dir, case_name, 'changrad.nc')
        ldd_file = os.path.join(self.case_dir, case_name, 'ldd.nc')
        outputfile = os.path.join(self.out_path_run, 'profiler.csv')

        profiler(input_file, ldd_file, x_coord, y_coord, coords_names='None')
        
        # generate the profiler
        profiler_data = profiler(input_file, ldd_file, x_coord, y_coord, coords_names)
        profiler_data.to_csv(outputfile)
        
        # compare the generated mask with the reference one
        ref_file = self.out_path_ref+'/profiler.csv'
        reference = pd.read_csv(ref_file)
        out_file = self.out_path_run+'/profiler.csv'
        generated = pd.read_csv(out_file)
        all_equal = reference.compare(generated)
        
        if len(all_equal)!=0:
            fail_message = f'Test profiler generation for {case_name} failed. Please check differences between' 
            fail_message += f' the generated mask "{out_file}" and the expected mask "{ref_file}".'
            assert all_equal, fail_message
        else:
            # if equal just delete the out folder
            shutil.rmtree(self.out_path_run, ignore_errors=True)


class TestProfilers(TestProfiler):

    def test_mctrivers_etrs89(self):
        self.run(4127500, 2412500, 'None', 'LF_ETRS89_UseCase')

    def test_mctrivers_latlon(self):
        self.run(-123.45, 53.95, ['lat', 'lon'], 'LF_lat_lon_UseCase')