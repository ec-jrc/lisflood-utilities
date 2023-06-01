import os
from shutil import copy2, rmtree
from datetime import datetime
from pathlib import Path
import numpy as np
from osgeo import gdal
from netCDF4 import Dataset
from lisfloodutilities.gridding.generate_grids import run, print_msg
from lisfloodutilities.gridding.lib.utils import FileUtils


cur_folder = os.path.dirname(os.path.realpath(__file__))


class TestGridding:

    def compare_netcdfs(self, nc1: Dataset, nc2: Dataset):
        assert nc1.variables['time'][:].size == nc2.variables['time'][:].size, 'Time axis do not have same size as reference.'
        coords_nc1_size = (nc1.variables['lat'][:].size, nc1.variables['lon'][:].size)
        coords_nc2_size = (nc2.variables['lat'][:].size, nc2.variables['lon'][:].size)
        assert coords_nc1_size == coords_nc2_size, 'Lat/Long axis do not have same size as reference.'
        max_timesteps = nc1.variables['time'][:].size
        for i in range(max_timesteps):
            assert np.allclose(nc1.variables['pr6'][i, :, :], nc2.variables['pr6'][i, :, :], atol=0.0000001), f'Timestep {i} is not equal to reference.'

    def compare_tiffs(self, tiff_folder_path1: Path, tiff_folder_path2: Path):
        for tiff_file2 in sorted(tiff_folder_path2.rglob('*.tiff')):
            ds2 = gdal.Open(str(tiff_file2))
            values2 = ds2.GetRasterBand(1).ReadAsArray()
            ds2 = None
            values1 = np.zeros(values2.shape)
            print(f'Testing {tiff_file2.name}')
            for tiff_file1 in tiff_folder_path1.rglob(tiff_file2.name):
                ds1 = gdal.Open(str(tiff_file1))
                values1 = ds1.GetRasterBand(1).ReadAsArray()
                ds1 = None
            assert values2.size == values2.size, 'Grid {tiff_file2.name} do not have same size as reference.'
            assert np.allclose(values1, values2, atol=0.0000001), f'File {tiff_file2.name} is not equal to reference.'

    def test_generate_netcdf(self):
        input_folder = 'tests/data/gridding/meteo_in/test1'
        reference_output = 'tests/data/gridding/reference/test1/pr6.nc'
        output_netcdf = 'tests/data/gridding/meteo_out/test1/pr6.nc'

        Path('tests/data/gridding/meteo_out/test1').mkdir(parents=True, exist_ok=True)

        if os.path.exists(output_netcdf):
            os.remove(output_netcdf)

        quiet_mode = True
        variable_code = 'pr6'
        config_type = '1arcmin'
        start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = False
        infolder = input_folder
        overwrite_output = False
        outfolder_or_file = output_netcdf
        processing_dates_file = None


        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, overwrite_output, start_date, end_date)

        reference = Dataset(reference_output)
        out = Dataset(output_netcdf)

        self.compare_netcdfs(out, reference)

        out = None
        reference = None

    def test_update_netcdf(self):
        input_folder = 'tests/data/gridding/meteo_in/test2'
        reference_output = 'tests/data/gridding/reference/test2/pr6.nc'
        output_netcdf = 'tests/data/gridding/meteo_out/test2/pr6.nc'
        processing_dates_file = 'tests/data/gridding/meteo_in/test2/dates_to_process.txt'

        Path('tests/data/gridding/meteo_out/test2').mkdir(parents=True, exist_ok=True)

        # This is needed so we can test updating 1 timestep on an existing netCDF file.
        copy2(os.path.join(input_folder, 'pr6.nc'), output_netcdf)

        quiet_mode = True
        variable_code = 'pr6'
        config_type = '1arcmin'
        start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = False
        infolder = input_folder
        overwrite_output = True
        outfolder_or_file = output_netcdf

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, overwrite_output, start_date, end_date)

        reference = Dataset(reference_output)
        out = Dataset(output_netcdf)

        self.compare_netcdfs(out, reference)

        out = None
        reference = None

    def test_generate_tiff(self):
        input_folder = 'tests/data/gridding/meteo_in/test1'
        reference_output = 'tests/data/gridding/reference/test3'
        output_folder = 'tests/data/gridding/meteo_out/test3'

        Path(output_folder).mkdir(parents=True, exist_ok=True)

        for filename in os.listdir(output_folder):
            file_path = os.path.join(output_folder, filename)
            if os.path.exists(file_path):
                rmtree(file_path)

        quiet_mode = True
        variable_code = 'pr6'
        config_type = '1arcmin'
        start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = True
        infolder = input_folder
        overwrite_output = False
        outfolder_or_file = output_folder
        processing_dates_file = None

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, overwrite_output, start_date, end_date)

        print('compare_tiffs')
        self.compare_tiffs(Path(output_folder), Path(reference_output))

        out = None
        reference = None
