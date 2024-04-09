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

    def compare_netcdfs(self, nc1: Dataset, nc2: Dataset, variable_name: str='pr6'):
        time_list1 = nc1.variables['time'][:].tolist()
        time_list2 = nc2.variables['time'][:].tolist()
        assert len(time_list1) == len(time_list2), 'Time axis do not have same size as reference.'
        nc1_idx = time_list1[0]
        nc2_idx = time_list2[0]
        assert nc1_idx == nc2_idx, f'Initial timestep is different. nc1.time[0]={nc1_idx} nc2.time[0]={nc2_idx}'
        coords_nc1_size = (nc1.variables['lat'][:].size, nc1.variables['lon'][:].size)
        coords_nc2_size = (nc2.variables['lat'][:].size, nc2.variables['lon'][:].size)
        assert coords_nc1_size == coords_nc2_size, 'Lat/Long axis do not have same size as reference.'
        max_timesteps = len(time_list1)
        for i in range(max_timesteps):
            assert np.allclose(nc1.variables[variable_name][i, :, :], nc2.variables[variable_name][i, :, :], atol=0.0000001), f'Timestep {i} is not equal to reference.'

    def compare_tiffs(self, tiff_folder_path1: Path, tiff_folder_path2: Path):
        for tiff_file2 in sorted(tiff_folder_path2.rglob('*.tiff')):
            ds2 = gdal.Open(str(tiff_file2))
            values2 = ds2.GetRasterBand(1).ReadAsArray()
            ds2 = None
            values1 = None
            for tiff_file1 in sorted(tiff_folder_path1.rglob(tiff_file2.name)):
                ds1 = gdal.Open(str(tiff_file1))
                values1 = ds1.GetRasterBand(1).ReadAsArray()
                ds1 = None
            assert values1 is not None, f'Grid {tiff_file2.name} was not generated.'
            assert values1.size == values2.size, f'Grid {tiff_file2.name} do not have same size as reference.'
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
        out_netcdf = True
        infolder = input_folder
        overwrite_output = False
        use_existing_file = True
        get_existing_tiff = False
        outfolder_or_file = output_netcdf
        processing_dates_file = None
        interpolation_mode = 'adw'
        use_broadcasting = True
        memory_save_mode = '5'

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, out_netcdf, overwrite_output, use_existing_file, get_existing_tiff, start_date, end_date,
            interpolation_mode=interpolation_mode, use_broadcasting=use_broadcasting, memory_save_mode=memory_save_mode)

        reference = Dataset(reference_output)
        out = Dataset(output_netcdf)

        self.compare_netcdfs(out, reference)

        out = None
        reference = None

    def test_generate_netcdf_without_start_date(self):
        input_folder = 'tests/data/gridding/meteo_in/test1'
        reference_output = 'tests/data/gridding/reference/test1/pr6_without_start_date.nc'
        output_netcdf = 'tests/data/gridding/meteo_out/test1/pr6_without_start_date.nc'

        Path('tests/data/gridding/meteo_out/test1').mkdir(parents=True, exist_ok=True)

        if os.path.exists(output_netcdf):
            os.remove(output_netcdf)

        quiet_mode = True
        variable_code = 'pr6'
        config_type = '1arcmin'
        # Not defining start_date it will use the netCDF timeunits
        # start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = False
        out_netcdf = True
        infolder = input_folder
        overwrite_output = False
        use_existing_file = True
        get_existing_tiff = True
        outfolder_or_file = output_netcdf
        processing_dates_file = None
        interpolation_mode = 'adw'
        use_broadcasting = True
        memory_save_mode = '5'

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = None # datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, out_netcdf, overwrite_output, use_existing_file, get_existing_tiff, start_date, end_date,
            interpolation_mode=interpolation_mode, use_broadcasting=use_broadcasting, memory_save_mode=memory_save_mode)

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
        out_netcdf = True
        infolder = input_folder
        overwrite_output = True
        use_existing_file = True
        get_existing_tiff = False
        outfolder_or_file = output_netcdf
        interpolation_mode = 'adw'
        use_broadcasting = False
        memory_save_mode = '5'

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, out_netcdf, overwrite_output, use_existing_file, get_existing_tiff, start_date, end_date,
            interpolation_mode=interpolation_mode, use_broadcasting=use_broadcasting, memory_save_mode=memory_save_mode)

        reference = Dataset(reference_output)
        out = Dataset(output_netcdf)

        self.compare_netcdfs(out, reference)

        out = None
        reference = None

    def test_generate_tiff(self):
        input_folder = 'tests/data/gridding/meteo_in/test1'
        reference_output = 'tests/data/gridding/reference/test3'
        output_folder = 'tests/data/gridding/meteo_out/test3'
        
        output_folder_path = Path(output_folder)

        output_folder_path.mkdir(parents=True, exist_ok=True)
        
        # Clean output folder
        for output_file in sorted(output_folder_path.rglob('*.*')):
            os.remove(output_file)

        quiet_mode = True
        variable_code = 'pr6'
        config_type = '1arcmin'
        start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = True
        out_netcdf = False
        infolder = input_folder
        overwrite_output = True
        use_existing_file = True
        get_existing_tiff = False
        outfolder_or_file = str(Path(output_folder, 'output.nc'))
        processing_dates_file = None
        interpolation_mode = 'adw'
        use_broadcasting = True
        memory_save_mode = '5'

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, out_netcdf, overwrite_output, use_existing_file, get_existing_tiff, start_date, end_date,
            interpolation_mode=interpolation_mode, use_broadcasting=use_broadcasting, memory_save_mode=memory_save_mode)

        # Remove netcdf output file since we are only interested in the output tiffs
        if os.path.exists(outfolder_or_file):
            os.remove(outfolder_or_file)
        
        # Move the output tiffs from the input folder into the output folder
        for tiff_file in sorted(Path(input_folder).rglob('*.tiff')):
            out_tiff_file = Path(output_folder, tiff_file.name)
            tiff_file.rename(out_tiff_file)

        self.compare_tiffs(output_folder_path, Path(reference_output))
        
        out = None
        reference = None

    def test_generate_netcdf_hourly(self):
        input_folder = 'tests/data/gridding/meteo_in/test4'
        reference_output = 'tests/data/gridding/reference/test4/pr1.nc'
        output_netcdf = 'tests/data/gridding/meteo_out/test4/pr1.nc'

        Path('tests/data/gridding/meteo_out/test4').mkdir(parents=True, exist_ok=True)

        if os.path.exists(output_netcdf):
            os.remove(output_netcdf)

        quiet_mode = True
        variable_code = 'pr1'
        config_type = '1arcmin'
        # Not defining start_date it will use the netCDF timeunits
        # start_date_str = '202303131200'
        end_date_str = '202303160000'
        out_tiff = False
        out_netcdf = True
        infolder = input_folder
        overwrite_output = False
        use_existing_file = True
        get_existing_tiff = False
        outfolder_or_file = output_netcdf
        processing_dates_file = None
        interpolation_mode = 'adw'
        use_broadcasting = True
        memory_save_mode = '5'

        configuration_base_folder = os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        start_date = None # datetime.strptime(start_date_str, FileUtils.DATE_PATTERN_CONDENSED)
        end_date = datetime.strptime(end_date_str, FileUtils.DATE_PATTERN_CONDENSED)

        run(config_filename, infolder, outfolder_or_file, processing_dates_file,
            file_utils, out_tiff, out_netcdf, overwrite_output, use_existing_file, get_existing_tiff, start_date, end_date,
            interpolation_mode=interpolation_mode, use_broadcasting=use_broadcasting, memory_save_mode=memory_save_mode)

        reference = Dataset(reference_output)
        out = Dataset(output_netcdf)

        self.compare_netcdfs(out, reference, variable_name=variable_code)

        out = None
        reference = None
