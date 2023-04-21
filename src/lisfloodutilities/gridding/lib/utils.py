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
from scipy.spatial import cKDTree
from pyg2p import Loggable
from pyg2p.main.readers.netcdf import NetCDFReader
from pyg2p.main.interpolation.scipy_interpolation_lib import ScipyInterpolation

__DECIMAL_CASES = 20
__DECIMAL_FORMAT = '{:.20f}'

getcontext().prec = __DECIMAL_CASES


class Printable(Loggable):

    def __init__(self, quiet_mode: bool = True):
        super().__init__()
        self.quiet_mode = quiet_mode

    def print_msg(self, msg: str = ''):
        if not self.quiet_mode:
            print(msg)
        self._log(msg)


class Dem(Printable):
    def __init__(self, dem_map: Path):
        super().__init__()
        self._dem_map = dem_map
        self.print_msg(f'Reading grid settings and altitude values from: {dem_map}')
        reader = NetCDFReader(self._dem_map)
        self.nrows = reader._rows
        self.ncols = reader._cols
        self.mv = reader.mv
        self.values = reader.values
        self.lats, self.lons = reader.get_lat_lon_values()
        self.lat_values = reader.get_lat_values()
        self.lon_values = reader.get_lon_values()
        self.cell_size_x = reader._pxlW
        self.cell_size_y = reader._pxlH
        self.lat_min = reader.lat_min
        self.lon_min = reader.lon_min
        self.lat_max = reader.lat_max
        self.lon_max = reader.lon_max


class FileUtils(Printable):
    DATE_PATTERN_CONDENSED = '%Y%m%d%H%M%S'
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
                        raise ArgumentTypeError(f"Inside file {file_path}. The file name [ {line} ] is not valid")
        if len(dates_list) > 0:
            return dict.fromkeys(dates_list)
        return {}

    def get_config_type_path(self, base_path: str, config_type: str) -> str:
        search_path = os.path.join(base_path, config_type)
        self.print_msg(f'Config base path: {search_path}')
        if not config_type:
            raise ArgumentTypeError("You must insert a valid grid configuration.")
        if not os.path.isdir(search_path):
            raise ArgumentTypeError("%s is not a valid grid configuration" % (config_type))
        return search_path

    def get_config_file(self, base_path: str) -> str:
        if not self.var_code:
            raise ArgumentTypeError("You must insert a valid variable.")
        search_path = os.path.join(base_path, 'config_%s.txt' % self.var_code)
        if not os.path.isfile(search_path):
            raise ArgumentTypeError("%s is not a valid variable code." % (self.var_code))
        return search_path

    def file_type(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError("You must insert a path.")
        if not os.path.isfile(fname):
            raise ArgumentTypeError("%s is not a valid path" % (fname))
        return fname

    def folder_type(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError("You must insert a path.")
        if not os.path.isdir(fname):
            raise ArgumentTypeError("%s is not a valid path" % (fname))
        return fname

    def is_file_path(file_path: Path) -> bool:
        return len(file_path.suffix) > 0

    def file_or_folder(fname: str) -> str:
        if not fname:
            raise ArgumentTypeError("You must insert a path.")
        p = Path(fname)
        if not (p.is_dir() or (FileUtils.is_file_path(p) and p.parent.is_dir())):
            raise ArgumentTypeError("%s is not a valid path" % (fname))
        return fname


class Config(Printable):
    def __init__(self, config_filename: str, start_date: datetime = None,
                 end_date: datetime = None, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self.COLUMN_HEIGHT = "height"
        self.COLUMN_VALUE = "value"
        self.COLUMN_LON = "lon"
        self.COLUMN_LAT = "lat"
        self.VALUE_NAN = -9999.0
        self.__setup_variable_config(config_filename)
        self.__HEIGHT_CORRECTION_FACTOR = {"tx": 0.006, "tn": 0.006, "ta": 0.006, "ta6": 0.006, "pd": 0.00025}

        self._dem_file = Path(self.config_path).joinpath('dem.nc')
        self._dem_interpolation_file = Path(self.config_path).joinpath('dem.wgs')
        self._dem_height_correction_file = Path(self.config_path).joinpath('dem.npy')
        self._dem_height_correction_kdt_file = Path(self.config_path).joinpath('dem.kdt')

        self.__setup_dem()

        # self.interpolation_mode = 'nearest'
        self.interpolation_mode = 'invdist'
        # self.interpolation_mode = 'bilinear'
        # self.interpolation_mode = 'triangulation'
        # self.interpolation_mode = 'bilinear_delaunay'
        self.scipy_modes_nnear = {'nearest': 1, 'invdist': 4, 'bilinear': 4, 'triangulation': 3, 'bilinear_delaunay': 4}
        self.grid_details = {'gridType': 'NOT rotated', 'Nj': 1, 'radius': 1.0}

        self.__setup_data_structures()
        self.__setup_dates(start_date, end_date)

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
        self.dem = Dem(self._dem_file)
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

    def get_config_field(self, config_group: str = '', config_property: str = ''):
        return self.__configFile.get(config_group, config_property)

    @property
    def scale_factor(self):
        return float(self.get_config_field('PROPERTIES', 'VALUE_SCALE'))

    @property
    def add_offset(self):
        return float(self.get_config_field('PROPERTIES', 'VALUE_OFFSET'))

    @property
    def var_code(self):
        return self.get_config_field('PROPERTIES', 'VAR_CODE')

    @property
    def do_height_correction(self):
        return self.var_code in self.__HEIGHT_CORRECTION_FACTOR

    @property
    def height_correction_factor(self):
        return self.__HEIGHT_CORRECTION_FACTOR[self.var_code]

    @property
    def neighbours_near(self):
        return self.scipy_modes_nnear[self.interpolation_mode]


class GriddingUtils(Printable):

    def __init__(self, conf: Config, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self.conf = conf

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
        result[np.where(result!=self.conf.VALUE_NAN)] /= self.conf.scale_factor
        result[np.where(result!=self.conf.VALUE_NAN)] += self.conf.add_offset
        # Reshape grid
        grid_data = result.reshape(grid_shape)
        return grid_data

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
        xp = np.array(x)
        yp = np.array(y)
        values = np.array(z)
        scipy_interpolation = ScipyInterpolation(xp, yp, self.conf.grid_details, values,
                                                 self.conf.neighbours_near, mv_target, mv_source,
                                                 target_is_rotated=False, parallel=False,
                                                 mode=self.conf.interpolation_mode)
        results, weights, indexes = scipy_interpolation.interpolate(grid_x, grid_y)
        result = np.ma.masked_array(data=results, mask=self.conf.dem_mask, fill_value=self.conf.VALUE_NAN)
        self.print_msg('Finish interpolation')
        result = self.reset_height(result)
        grid_data = self.prepare_grid(result, grid_x.shape)
        return grid_data
