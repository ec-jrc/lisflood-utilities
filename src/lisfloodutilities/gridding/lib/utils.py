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

import sys
import os
import csv
import pickle
from pathlib import Path
# from pyproj import Transformer
from argparse import ArgumentTypeError
import configparser as ConfigParser
import numpy as np
import numpy.ma as ma
import pandas as pd
from decimal import *
from datetime import datetime, timedelta
from collections import OrderedDict
from scipy.spatial import cKDTree
from pyg2p import Loggable
from pyg2p.main.readers.netcdf import NetCDFReader
from pyg2p.main.interpolation.scipy_interpolation_lib import ScipyInterpolation
from numpy import delete
import importlib
from lisfloodutilities.gridding.lib.filters import KiwisFilter

from .filters import *


__DECIMAL_CASES = 20
__DECIMAL_FORMAT = '{:.20f}'

getcontext().prec = __DECIMAL_CASES


class StackFIFO:
    def __init__(self):
        self.stack = []

    def push(self, item):
        self.stack.append(item)

    def pop(self):
        if self.is_empty():
            return None
        else:
            return self.stack.pop(0)

    def as_list(self) -> list:
        return self.stack

    def is_empty(self) -> bool:
        return len(self.stack) == 0


class Printable(Loggable):

    def __init__(self, quiet_mode: bool = True):
        super().__init__()
        self.quiet_mode = quiet_mode

    def print_msg(self, msg: str = ''):
        if not self.quiet_mode:
            print(msg)
        self._log(msg)


class Dem(Printable):
    def __init__(self, dem_map: Path, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self._dem_map = dem_map
        self.print_msg(f'Reading grid settings and altitude values from: {dem_map}')
        reader = NetCDFReader(self._dem_map)
        self.nrows = reader._rows
        self.ncols = reader._cols
        self.mv = reader.mv.astype(np.float32)
        self.values = reader.values.astype(np.float32)
        self.lats, self.lons = reader.get_lat_lon_values()
        self.lats = self.lats.astype(np.float32)
        self.lons = self.lons.astype(np.float32) 
        self.lat_values = reader.get_lat_values()
        self.lon_values = reader.get_lon_values()
        self.cell_size_x = reader._pxlW
        self.cell_size_y = reader._pxlH
        self.lat_min = reader.lat_min
        self.lon_min = reader.lon_min
        self.lat_max = reader.lat_max
        self.lon_max = reader.lon_max
        self.count_nan = np.count_nonzero(np.isnan(self.values))


class FileUtils(Printable):
    DATE_PATTERN_CONDENSED = '%Y%m%d%H%M%S'
    DATE_PATTERN_CONDENSED_SHORT = '%Y%m%d%H%M'
    DATE_PATTERN_SEPARATED = '%Y-%m-%d %H:%M:%S'
    CSV_DELIMITER = '\t'
    FILES_WILDCARD = '??????????00_??????????????.txt'

    def __init__(self, var_code: str, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self.var_code = var_code
        self.var_size = len(self.var_code)

    def get_timestamp_from_filename(self, filename: Path) -> datetime:
        file_timestamp = filename.name[self.var_size:12+self.var_size] + '00'
        return datetime.strptime(file_timestamp, FileUtils.DATE_PATTERN_CONDENSED)

    def processable_file(self, file_timestamp: datetime, dates_to_process: dict = {},
                         start_date: datetime = None, end_date: datetime = None) -> bool:
        is_processable_date = True
        if len(dates_to_process) > 0:
            is_processable_date = file_timestamp in dates_to_process
        return (start_date is not None and start_date <= file_timestamp and
                end_date is not None and file_timestamp <= end_date and is_processable_date)

    def read_processing_dates_file(self, file_path: str = None) -> dict:
        dates_list = []
        if file_path is not None:
            with open(file_path) as in_file:
                for line in in_file:
                    line = line.strip()
                    try:
                        cur_date = self.get_timestamp_from_filename(Path(line))
                        dates_list.append(cur_date)
                    except Exception as e:
                        raise ArgumentTypeError(f'Inside file {file_path}. The file name [ {line} ] is not valid')
        if len(dates_list) > 0:
            return dict.fromkeys(dates_list)
        return {}

    def get_config_type_path(self, base_path: str, config_type: str) -> str:
        search_path = os.path.join(base_path, config_type)
        self.print_msg(f'Config base path: {search_path}')
        if not config_type:
            raise ArgumentTypeError('You must insert a valid grid configuration.')
        if not os.path.isdir(search_path):
            raise ArgumentTypeError(f'{config_type} is not a valid grid configuration')
        return search_path

    def get_config_file(self, base_path: str) -> str:
        if not self.var_code:
            raise ArgumentTypeError('You must insert a valid variable.')
        search_path = os.path.join(base_path, f'config_{self.var_code}.txt')
        if not os.path.isfile(search_path):
            raise ArgumentTypeError(f'{self.var_code} is not a valid variable code.')
        return search_path

    def file_type(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError('You must insert a path.')
        p = Path(fname)
        if not (p.parent.is_dir() and os.access(p.parent, os.W_OK)):
            raise ArgumentTypeError(f'Folder is not writable: {p.parent}')
        if not FileUtils.is_file_path(p):
            raise ArgumentTypeError(f'Invalid file path: {p}')
        return fname

    def folder_type(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError('You must insert a path.')
        if not os.path.isdir(fname):
            raise ArgumentTypeError(f'{fname} is not a valid path')
        return fname

    def is_file_path(file_path: Path) -> bool:
        return len(file_path.suffix) > 0

    def file_or_folder(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError('You must insert a path.')
        p = Path(fname)
        if not (p.is_dir() or (FileUtils.is_file_path(p) and p.parent.is_dir())):
            raise ArgumentTypeError(f'{fname} is not a valid path')
        return fname


class Config(Printable):
    # Each interpolation is paired with the number of neighbours used by the interpolation algorithm.
    INTERPOLATION_MODES = {'nearest': 1, 'invdist': 20, 'adw': 11, 'cdd': 11,
                           'bilinear': 4, 'triangulation': 3, 'bilinear_delaunay': 4}
    # Allows splitting the grid to interpolate into several parts and run isolated interpolations for each
    MEMORY_SAVE_MODES = {'0': None, '1': 2, '2': 4, '3': 10, '4': 20, '5': 50}
    def __init__(self, config_filename: str, start_date: datetime = None, end_date: datetime = None,
                 quiet_mode: bool = False, interpolation_mode: str = 'adw', memory_save_mode: str = '0'):
        super().__init__(quiet_mode)
        self.COLUMN_HEIGHT = "height"
        self.COLUMN_VALUE = "value"
        self.COLUMN_LON = "lon"
        self.COLUMN_LAT = "lat"
        self.VALUE_NAN = -9999.0
        self.interpolation_mode = interpolation_mode
        self.memory_save_mode = memory_save_mode
        self.__setup_variable_config(config_filename)
        self.__HEIGHT_CORRECTION_FACTOR = {"tx": 0.006, "tn": 0.006, "ta": 0.006, "ta6": 0.006, "pd": 0.00025}

        self._dem_file = Path(self.config_path).joinpath('dem.nc')
        self._dem_interpolation_file = Path(self.config_path).joinpath('dem.wgs')
        self._dem_height_correction_file = Path(self.config_path).joinpath('dem.npy')
        self._dem_height_correction_kdt_file = Path(self.config_path).joinpath('dem.kdt')

        self.__setup_dem()
        self.__setup_interpolation_parameters()
        self.__setup_data_structures()
        self.__setup_dates(start_date, end_date)

    def __setup_interpolation_parameters(self):
        self.scipy_modes_nnear = Config.INTERPOLATION_MODES
        self.grid_details = {'gridType': 'NOT rotated', 'Nj': 1, 'radius': 6367470.0}
        self.min_upper_bound = 100000000.0 # max search distance in meters
        cdd_map_path = Path(self.config_path).joinpath(f'CDDmap_{self.var_code}.nc')
        self.cdd_map = f'{cdd_map_path}'
        if self.interpolation_mode == 'cdd' and not os.path.isfile(self.cdd_map):
            raise ArgumentTypeError(f'CDD map was not found: {self.cdd_map}')
        self.cdd_mode = 'MixHofstraShepard'
        # 1) m=4; r=1/3; min 4 stations
        # 2) m=8, r=1/3; min 4 stations
        # 3) m=8; r=2/3; min 4 stations
        # 4) m=4; r=1/3; min 4 stations + CDD map without Euro4m, CarpatClim
        self.cdd_options = {
            'm_const': 4,
            'min_num_of_station': 4,
            'radius_ratio': 1/3,
            'weights_mode': 'All' # 'OnlyTOP10'
        }

    def __setup_dates(self, start_date: datetime = None, end_date: datetime = None):
        netcdf_var_time_unit_pattern = self.get_config_field('VAR_TIME','UNIT_PATTERN')
        if start_date is None:
            time_units = self.get_config_field('VAR_TIME','UNIT')
            self.start_date = datetime.strptime(time_units, netcdf_var_time_unit_pattern)
        else:
            self.start_date = start_date

        if end_date is None:
            self.end_date = datetime.strptime(datetime.now().strftime('%Y%m%d060000'), FileUtils.DATE_PATTERN_CONDENSED)
        else:
            self.end_date = end_date

    def __setup_variable_config(self, config_filename: str):
        self.config_path = os.path.dirname(config_filename)
        self.print_msg(f'Loading configuration from: {self.config_path}')
        self.__configFile = ConfigParser.ConfigParser()
        default_config = os.path.join(self.config_path, 'default.txt')
        if os.path.exists(default_config):
            self.__configFile.read(default_config)
        self.__configFile.read(config_filename)

    def __setup_dem(self):
        self.dem = Dem(self._dem_file, self.quiet_mode)
        self.dem_nrows = self.dem.nrows
        self.dem_ncols = self.dem.ncols
        self.dem_lons = self.dem.lons.flatten()
        self.dem_lats = self.dem.lats.flatten()
        self.dem_heights = self.dem.values.flatten()
        self.dem_cell_size_x = self.dem.cell_size_x
        self.dem_cell_size_y = self.dem.cell_size_y
        self.dem_min_x = self.dem.lon_min
        self.dem_max_x = self.dem.lon_max
        self.dem_min_y = self.dem.lat_min
        self.dem_max_y = self.dem.lat_max
        self.CELL_SIZE_X = Decimal(self.dem_cell_size_x)
        self.CELL_SIZE_Y = Decimal(self.dem_cell_size_y)
        self.COORDINATE_MIN_X = Decimal(self.dem_min_x) - self.CELL_SIZE_X/Decimal(2.0)
        self.COORDINATE_MAX_X = Decimal(self.dem_max_x)
        self.COORDINATE_MIN_Y = Decimal(self.dem_min_y)
        self.COORDINATE_MAX_Y = Decimal(self.dem_max_y) - self.CELL_SIZE_Y/Decimal(2.0)
        self.dem_mask = ma.getmask(ma.masked_where(np.isnan(self.dem.values), self.dem.values, copy=True))

    def __setup_data_structures(self):
        if not self._dem_interpolation_file.exists():
            self.print_msg('Creating data structures')
            df_dem = pd.DataFrame(
                {
                    self.COLUMN_LON: self.dem_lons,
                    self.COLUMN_LAT: self.dem_lats,
                    self.COLUMN_HEIGHT: self.dem_heights,
                }
            )
            df_dem = df_dem.dropna()
            df_dem = df_dem[df_dem[self.COLUMN_HEIGHT] >= -9000]
            self.print_msg('Create Height from DEM')
            self.DEM_HEIGHT_CORRECTION = np.ascontiguousarray(df_dem[self.COLUMN_HEIGHT].values)
            np.save(self._dem_height_correction_file, self.DEM_HEIGHT_CORRECTION)
            self.print_msg(f'Saved height from DEM: {self._dem_height_correction_file}')
            self.print_msg('Create KDtree DEM')
            transformed_coordinates = (df_dem[self.COLUMN_LON].values, df_dem[self.COLUMN_LAT].values)
            self.DEM_HEIGHT_CORRECTION_QUERY = cKDTree(data=np.vstack(transformed_coordinates).T, copy_data=True)
            pickle.dump(self.DEM_HEIGHT_CORRECTION_QUERY, self._dem_height_correction_kdt_file.open("wb"))
            self.print_msg(f'Saved creating KDtree DEM: {self._dem_height_correction_kdt_file}')
            self.print_msg(f'Writing file: {self._dem_interpolation_file}')
            df_dem.to_csv(
                self._dem_interpolation_file,
                index=False,
                header=False,
                float_format="%.8f",
                columns=[self.COLUMN_LON, self.COLUMN_LAT, self.COLUMN_HEIGHT],
                sep="\t",
            )
            self.print_msg('Finished data structures')
        else:
            self.print_msg('Loading data structures')
            self.print_msg(f'Reading file: {self._dem_interpolation_file}')
            df_dem = pd.read_csv(self._dem_interpolation_file, sep="\t", names=[self.COLUMN_LON, self.COLUMN_LAT, self.COLUMN_HEIGHT])
            self.print_msg(f'Reading file: {self._dem_height_correction_file}')
            self.DEM_HEIGHT_CORRECTION = np.load(self._dem_height_correction_file)
            self.print_msg(f'Reading file: {self._dem_height_correction_kdt_file}')
            self.DEM_HEIGHT_CORRECTION_QUERY = pickle.load(self._dem_height_correction_kdt_file.open("rb"))
            self.print_msg('Finish loading data structures')

    def get_config_field(self, config_group: str = '', config_property: str = '') -> str:
        return self.__configFile.get(config_group, config_property)

    @property
    def scale_factor(self) -> float:
        return float(self.get_config_field('PROPERTIES', 'VALUE_SCALE'))

    @property
    def add_offset(self) -> float:
        return float(self.get_config_field('PROPERTIES', 'VALUE_OFFSET'))

    @property
    def value_min(self) -> int:
        return int(self.get_config_field('PROPERTIES', 'VALUE_MIN'))

    @property
    def value_max(self) -> int:
        return int(self.get_config_field('PROPERTIES', 'VALUE_MAX'))

    @property
    def value_min_packed(self) -> int:
        return int((self.value_min - self.add_offset) / self.scale_factor)

    @property
    def value_max_packed(self) -> int:
        return int((self.value_max - self.add_offset) / self.scale_factor)

    @property
    def value_nan_packed(self) -> float:
        return float((self.VALUE_NAN - self.add_offset) / self.scale_factor)

    @property
    def var_code(self) -> str:
        return self.get_config_field('PROPERTIES', 'VAR_CODE')

    @property
    def do_height_correction(self) -> bool:
        return self.var_code in self.__HEIGHT_CORRECTION_FACTOR

    @property
    def height_correction_factor(self) -> float:
        return self.__HEIGHT_CORRECTION_FACTOR[self.var_code]

    @property
    def neighbours_near(self) -> int:
        return self.scipy_modes_nnear[self.interpolation_mode]

    @property
    def num_of_splits(self) -> int:
        return Config.MEMORY_SAVE_MODES[self.memory_save_mode]

    @property
    def input_wildcard(self) -> str:
        default_wildcard = self.get_config_field('GENERIC', 'INPUT_WILDCARD')
        return f'{self.var_code}{default_wildcard}'

    @property
    def input_timestamp_pattern(self) -> str:
        default_pattern = self.get_config_field('GENERIC', 'INPUT_TIMESTAMP_PATTERN')
        return f'{self.var_code}{default_pattern}'


class GriddingUtils(Printable):

    def __init__(self, conf: Config, quiet_mode: bool = False, use_broadcasting: bool = False):
        super().__init__(quiet_mode)
        self.conf = conf
        self.use_broadcasting = use_broadcasting
        self.unit_conversion = float(self.conf.get_config_field('PROPERTIES', 'UNIT_CONVERSION'))
        # self.compressed_NaN = (self.conf.VALUE_NAN - self.conf.add_offset) / self.conf.scale_factor

    def correct_height(self, df: pd.DataFrame) -> pd.DataFrame:
        if self.conf.do_height_correction:
            self.print_msg('Start height correction')
            if "hs" not in df:
                stations = np.array([df[self.conf.COLUMN_LON], df[self.conf.COLUMN_LAT]]).T
                _, idx = self.conf.DEM_HEIGHT_CORRECTION_QUERY.query(stations)
                df["hs"] = self.conf.DEM_HEIGHT_CORRECTION[idx]
            df[self.conf.COLUMN_VALUE] += df["hs"] * self.conf.height_correction_factor
            self.print_msg('Finish height correction')
        return df

    def reset_height(self, result: ma.MaskedArray) -> np.ndarray:
        if self.conf.do_height_correction:
            self.print_msg('Reset height correction')
            grid_cells_height = ma.masked_array(data=self.conf.dem_heights,
                                                mask=self.conf.dem_mask,
                                                fill_value=self.conf.VALUE_NAN)
            result -= grid_cells_height * self.conf.height_correction_factor
            if self.conf.var_code == "pd":
                # Eliminating the eventual negative values 
                # from Vapor Pressure resulting from height correction
                result[np.where(result<0)] = 0
            self.print_msg('Finish reseting height correction')
        result = result.filled()
        return result

    def prepare_grid(self, result: np.ndarray, grid_shape: np.ndarray) -> np.ndarray:
        # Compress data
        if self.unit_conversion != 1.0:
            result = result * self.unit_conversion
        result[np.where(result == self.conf.VALUE_NAN)] = np.nan
        result[np.where(result < self.conf.value_min)] = np.nan
        result[np.where(result > self.conf.value_max)] = np.nan
        result = np.round(result.astype(float), 1)
        result = (result - self.conf.add_offset) / self.conf.scale_factor
        result[np.isnan(result)] = self.conf.VALUE_NAN
        # Reshape grid
        grid_data = result.reshape(grid_shape)
        return grid_data

    def check_grid_nan(self, filename: Path, result: np.ndarray):
        # Verify if grid contains NaN different than the ones on the DEM.
        # Which means there are NaN values that shouldn't be generated.
        count_dem_nan = self.conf.dem.count_nan
        current_grid = result.copy()
        current_grid[np.where(current_grid == self.conf.VALUE_NAN)] = np.nan
        current_grid[np.where(current_grid < self.conf.value_min)] = np.nan
        current_grid[np.where(current_grid > self.conf.value_max)] = np.nan
        count_grid_nan = np.count_nonzero(np.isnan(current_grid))
        if count_dem_nan != count_grid_nan:
            print(f"WARNING: The grid interpolated from file {filename.name} contains different NaN values ({count_grid_nan}) than the DEM ({count_dem_nan}). diff: {count_grid_nan - count_dem_nan}")

    def generate_grid(self, filename: Path) -> np.ndarray:
        mv_target = self.conf.dem.mv
        mv_source = None
        grid_x, grid_y = self.conf.dem.lons, self.conf.dem.lats
        df = pd.read_csv(filename, header=None,
                         names=[self.conf.COLUMN_LON, self.conf.COLUMN_LAT, self.conf.COLUMN_VALUE],
                         na_values=self.conf.VALUE_NAN, sep=FileUtils.CSV_DELIMITER)
        df = self.correct_height(df)
        self.print_msg('Starting interpolation')
        x = df[self.conf.COLUMN_LON].values
        y = df[self.conf.COLUMN_LAT].values
        z = df[self.conf.COLUMN_VALUE].values
        xp = np.array(x).astype(np.float32)
        yp = np.array(y).astype(np.float32)
        values = np.array(z).astype(np.float32)
        df = None
        if self.conf.interpolation_mode == 'cdd':
            scipy_interpolation = ScipyInterpolation(xp, yp, self.conf.grid_details, values,
                                                     self.conf.neighbours_near, mv_target, mv_source,
                                                     target_is_rotated=False, parallel=False,
                                                     mode=self.conf.interpolation_mode,
                                                     cdd_map=self.conf.cdd_map,
                                                     cdd_mode=self.conf.cdd_mode,
                                                     cdd_options=self.conf.cdd_options,
                                                     use_broadcasting=self.use_broadcasting,
                                                     num_of_splits=self.conf.num_of_splits)
        else:
            scipy_interpolation = ScipyInterpolation(xp, yp, self.conf.grid_details, values,
                                                     self.conf.neighbours_near, mv_target, mv_source,
                                                     target_is_rotated=False, parallel=False,
                                                     mode=self.conf.interpolation_mode,
                                                     use_broadcasting=self.use_broadcasting,
                                                     num_of_splits=self.conf.num_of_splits)

        # reset search radius to a fixed value
        scipy_interpolation.min_upper_bound = self.conf.min_upper_bound
        results, weights, indexes = scipy_interpolation.interpolate(grid_x, grid_y)
        result = np.ma.masked_array(data=results, mask=self.conf.dem_mask, fill_value=self.conf.VALUE_NAN)
        self.print_msg('Finish interpolation')
        result = self.reset_height(result)
        self.check_grid_nan(filename, result)
        grid_data = self.prepare_grid(result, grid_x.shape)
        return grid_data


class KiwisLoader(Printable):

    def __init__(self, conf: Config, infolder: Path, dates_to_process: dict = {}, overwrite_file: bool = False, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self.conf = conf
        self.overwrite_file = overwrite_file
        self.var_code = self.conf.var_code
        self.var_size = len(self.var_code)
        self.inwildcard = self.conf.input_wildcard
        self.infolder = infolder
        self.dates_to_process = dates_to_process
        # Frequency between timesteps in hours
        self.time_frequency = int(self.conf.get_config_field('VAR_TIME','FREQUENCY'))
        self.is_daily_var = (self.time_frequency == 1)
        # Number of files to be read/processed simultaneously. for non daily vars time frequency is in hours
        self.read_files_step = 1 if self.is_daily_var else int(24 / self.time_frequency)
        self.files_list_dic = OrderedDict()
        self.files_list = None
        self.files_cache = StackFIFO()
        self.filter_classes = self.__get_filter_classes()

    def __iter__(self):
        self.__load_kiwis_paths()
        self.files_list = iter(self.files_list_dic)
        return self

    def __next__(self):
        # Loads cache if empty
        if self.files_cache.is_empty():
            self.__load_next_batch_of_files()
            self.__process_next_batch_of_files()
        # Get next file
        filepath_kiwis, filepath_points, kiwis_timestamps = self.files_cache.pop()
        return filepath_points

    def __load_next_batch_of_files(self):
        """
        Load next batch of files to be processed.
        If there are sub daily files, they will be processed all at the same time grouped by 24h.
        e.g.: 6 hourly var will process 4 consecutive files covering 24h period. 
        """
        next_key = next(self.files_list)
        self.files_cache.push(self.files_list_dic[next_key])
        if not self.is_daily_var:
            # Read remaining 6hourly steps
            try:
                for i in range(self.read_files_step - 1):
                    next_key = next(self.files_list)
                    self.files_cache.push(self.files_list_dic[next_key])
            except StopIteration:
                pass
        # At this point after loading the cache if it is empty means there are no more files to process
        if self.files_cache.is_empty():
            raise StopIteration

    def __process_next_batch_of_files(self):
        """
        Create the point files by filtering the Kiwis using the filter classes.
        Processing all files from the same group, to allow the filtering classes to filter
        the same stations on all the files in the same 24h group.
        """
        files_group = self.files_cache.as_list()
        filepath_kiwis = [pair[0] for pair in files_group]
        filepath_points = [pair[1] for pair in files_group]
        kiwis_timestamps = [pair[2] for pair in files_group]
        df_kiwis_array = []
        # Allow us to get the output columns out of the dataframe
        last_filter_class = KiwisFilter()
        # At this point we need to process all the files nevertheless because even
        # if 1 single file needs to be generated it might affect the results of the
        # other existing ones, since one station that is filtered in one file might
        # be filtered on the remaining files simultaneously.
        for filter_class in self.filter_classes:
            self.print_msg(f'{filter_class.get_class_description()}')
            last_filter_class = filter_class
            df_kiwis_array = filter_class.filter(filepath_kiwis, kiwis_timestamps, df_kiwis_array)
        i = 0
        for df in df_kiwis_array:
            df_output = last_filter_class.get_dataframe_output_columns(df)
            filepath = filepath_points[i]
            if self.overwrite_file or not filepath.is_file():
                df_output.to_csv(
                    filepath,
                    index=False,
                    header=False,
                    sep="\t",
                )
            i += 1
        

    def __get_filter_classes(self) -> list:
        '''
        Create and instantiate the Filter classes to process the kiwis files
        '''
        plugins_array = []
        # Load the class dynamically
        plugins_config = self.conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_PLUGIN_CLASSES')
        plugins_dic = eval(plugins_config)
        plugins = OrderedDict(plugins_dic)
        plugins_columns_def = self.conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
        plugins_columns_dic = eval(plugins_columns_def)
        plugins_columns = OrderedDict(plugins_columns_dic)
        module_name = 'lisfloodutilities.gridding.lib.filters'
        try:
            for plugin in plugins:
                class_name = plugin
                class_args = plugins[plugin]
                module = importlib.import_module(module_name)
                class_instance = getattr(module, class_name)(plugins_columns, class_args)
                plugins_array.append(class_instance)
        except ImportError:
            print(f"Error: Could not import module '{module_name}'")
        except AttributeError:
            print(f"Error: Could not find class '{class_name}' in module '{module_name}'")
        return plugins_array

    def __get_points_filename(self, kiwis_timestamp: str, filename_kiwis: Path) -> Path:
        '''
        Returns the points file path.
        If the mode is overwrite tries to get the first pointfile path it finds and if it does not find generates a new file path.
        Otherwise generates a new file path.
        '''
        if self.overwrite_file:
            for points_path in sorted(filename_kiwis.parent.rglob(f'{self.var_code}{kiwis_timestamp}_??????????????.txt')):
                if points_path.is_file():
                    return points_path
        pointfile_timestamp = datetime.now().strftime(FileUtils.DATE_PATTERN_CONDENSED)
        return Path(filename_kiwis.parent, f'{self.var_code}{kiwis_timestamp}_{pointfile_timestamp}.txt')

    def __load_kiwis_paths(self):
        netcdf_offset_file_date = int(self.conf.get_config_field('VAR_TIME','OFFSET_FILE_DATE'))
        for filename_kiwis in sorted(self.infolder.rglob(self.inwildcard)):
            kiwis_timestamp = self.__get_timestamp_from_filename(filename_kiwis)
            file_timestamp = kiwis_timestamp + timedelta(days=netcdf_offset_file_date)
            if self.__processable_file(file_timestamp, self.conf.start_date, self.conf.end_date):
                kiwis_timestamp_str = kiwis_timestamp.strftime(FileUtils.DATE_PATTERN_CONDENSED_SHORT)
                filename_points = self.__get_points_filename(kiwis_timestamp_str, filename_kiwis)
                self.files_list_dic[kiwis_timestamp_str] = (filename_kiwis, filename_points, kiwis_timestamp_str)

    def __get_timestamp_from_filename(self, filename: Path) -> datetime:
        return datetime.strptime(filename.name, self.conf.input_timestamp_pattern)

    def __processable_file(self, file_timestamp: datetime, start_date: datetime = None, end_date: datetime = None) -> bool:
        is_processable_date = True
        if len(self.dates_to_process) > 0:
            is_processable_date = file_timestamp in self.dates_to_process
        return (start_date is not None and start_date <= file_timestamp and
                end_date is not None and file_timestamp <= end_date and is_processable_date)

