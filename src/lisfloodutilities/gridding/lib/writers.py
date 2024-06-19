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

import os
import time as timex
import warnings
import numpy as np
import copy
from argparse import ArgumentTypeError
from pathlib import Path
from datetime import datetime, timedelta
from netCDF4 import Dataset, default_fillvals, date2num, num2date
from netCDF4._netCDF4 import Variable as NetCDF4Variable
from osgeo import osr, gdal
from lisfloodutilities.gridding.lib.utils import Printable, Config, FileUtils
from lisfloodutilities import version


class OutputWriter(Printable):

    def __init__(self, conf: Config, overwrite_file: bool = False, quiet_mode: bool = False):
        super().__init__(quiet_mode)
        self.filepath = None
        self.overwrite_file = overwrite_file
        self.conf = conf
        self.__setup_metadata()

    def __setup_metadata(self):
        # General Attributes
        self.Source_Software = f'lisflood-utilities.gridding v{version} (interpolation: {self.conf.interpolation_mode})'
        self.reference = self.conf.get_config_field('GENERIC','NETCDF_REFERENCE')
        self.title = self.conf.get_config_field('GENERIC','NETCDF_TITLE')
        self.keywords = self.conf.get_config_field('GENERIC','NETCDF_KEYWORDS')
        self.source = self.conf.get_config_field('GENERIC','NETCDF_SOURCE')
        self.institution = self.conf.get_config_field('GENERIC','NETCDF_INSTITUTION')
        self.comment = self.conf.get_config_field('GENERIC','NETCDF_COMMENT')
        self.var_code = self.conf.get_config_field('PROPERTIES', 'VAR_CODE')
        self.var_standard_name = self.conf.get_config_field('PROPERTIES', 'STANDARD_NAME')
        self.var_long_name = self.conf.get_config_field('PROPERTIES', 'LONG_NAME')
        self.cell_methods = self.conf.get_config_field('PROPERTIES', 'CELL_METHODS')
        self.units = self.conf.get_config_field('PROPERTIES', 'UNIT')

    def setNaN(self, value, defaultNaN=np.nan):
        try:
            value[value==1e31] = defaultNaN
        except Exception as e:
            self.print_msg(f"value==1e31 : {str(e)}")
        try:
            value[value==self.conf.VALUE_NAN] = defaultNaN
        except Exception as e:
            self.print_msg(f"value=={self.conf.VALUE_NAN} : {str(e)}")
        try:
            value[value==-32768.0] = defaultNaN
        except Exception as e:
            self.print_msg(f"value==-32768.0 : {str(e)}")
        try:
            value[value==31082] = defaultNaN
        except Exception as e:
            self.print_msg(f"value==31082 : {str(e)}")
        return value

    def open(self, out_filename: Path):
        self.filepath = out_filename
        self.time_created = timex.ctime(timex.time())

    def print_grid_statistics(self, grid: np.ndarray):
        grid_min = np.nanmin(grid)
        grid_max = np.nanmax(grid)
        grid_mean = np.nanmean(grid)
        grid_percentile_10 = np.nanpercentile(grid, 10)
        grid_percentile_90 = np.nanpercentile(grid, 90)
        stats_string = (
            f'#APP_STATS: {{"TIMESTAMP": "{self.current_timestamp}", "VAR_CODE": "{self.conf.var_code}", '
            f'"MINIMUM_VALUE": {grid_min:.2f}, "MAXIMUM_VALUE": {grid_max:.2f}, '
            f'"MEAN_VALUE": {grid_mean:.2f}, "PERCENTILE_10": {grid_percentile_10:.2f}, '
            f'"PERCENTILE_90": {grid_percentile_90:.2f}}}'
        )
        self.print_msg(stats_string)
    
    def setup_grid(self, grid: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def write(self, grid: np.ndarray, timestamp: datetime = None, print_stats: bool = True):
        raise NotImplementedError

    def write_timestep(self, grid: np.ndarray, timestep: int = -1):
        raise NotImplementedError

    def opened(self) -> bool:
        return not self.filepath is None

    def close(self):
        self.filepath = None


class NetCDFWriter(OutputWriter):
    NETCDF_DATASET_FORMAT = 'NETCDF4_CLASSIC'
    NETCDF_CONVENTIONS = 'CF-1.6'
    NETCDF_VAR_DATA_TYPE = 'f8' # np.float64
    NETCDF_COORDINATES_DATA_TYPE = 'i4' # np.int32
    NETCDF_VAR_TIME_CALENDAR_TYPE = 'proleptic_gregorian'
    COMPRESSION_LEVEL = 4

    def __init__(self, conf: Config, overwrite_file: bool = False, quiet_mode: bool = False):
        super().__init__(conf, overwrite_file, quiet_mode)
        self.nf = None
        self.netcdf_var_time_unit_pattern = self.conf.get_config_field('VAR_TIME','UNIT_PATTERN')
        self.netcdf_var_time = self.conf.get_config_field('VAR_TIME','NAME')
        self.time_frequency = int(self.conf.get_config_field('VAR_TIME','FREQUENCY'))
        self.calendar_type = self.NETCDF_VAR_TIME_CALENDAR_TYPE
        self.calendar_time_unit = self.conf.start_date.strftime(self.netcdf_var_time_unit_pattern)
        # Write index not set yet
        self.write_idx = -1
        self.current_timestamp = None
        self.is_new_file = False

    def open(self, out_filename: Path):
        super().open(out_filename)
        if not os.path.isfile(self.filepath):
            self.nf = Dataset(self.filepath, 'w', format=self.NETCDF_DATASET_FORMAT)
            self.__setup_netcdf_metadata()
            self.is_new_file = True
        elif self.overwrite_file:
            self.nf = Dataset(self.filepath, 'r+', clobber=True, format=self.NETCDF_DATASET_FORMAT)
        else:
            raise ArgumentTypeError(f'File {self.filepath} already exists. Use --force flag to append.')

    def setup_grid(self, grid: np.ndarray, print_stats: bool = True) -> np.ndarray:
        values = self.setNaN(copy.deepcopy(grid))
        values[values < self.conf.value_min_packed] = np.nan
        values[values > self.conf.value_max_packed] = np.nan
        values[values != self.conf.VALUE_NAN] *= self.conf.scale_factor
        values[values != self.conf.VALUE_NAN] += self.conf.add_offset
        if print_stats:
            self.print_grid_statistics(values)
        values[np.isnan(values)] = self.conf.VALUE_NAN * self.conf.scale_factor + self.conf.add_offset
        return values   

    def write(self, grid: np.ndarray, timestamp: datetime = None, print_stats: bool = True):
        timestep = -1
        if timestamp is not None:
            self.current_timestamp = timestamp.strftime(FileUtils.DATE_PATTERN_SEPARATED)
            timestep = date2num(timestamp, self.calendar_time_unit, self.calendar_type)
        else:
            self.current_timestamp = None
        cur_grid = self.setup_grid(grid, print_stats)
        self.write_timestep(cur_grid, timestep)

    def write_timestep(self, grid: np.ndarray, timestep: int = -1):
        if timestep >= 0:
            if not self.opened():
                raise Exception("netCDF Dataset was not initialized. If file already exists, use --force flag to append.")
            self.__set_write_index(timestep)
            self.nf.variables[self.netcdf_var_time][self.write_idx] = timestep
            self.nf.variables[self.var_code][self.write_idx, :, :] = grid

    def __set_write_index(self, timestep: int):
        if not self.is_new_file:
            # Writing into an existing file need to calculate write index
            time_array = self.nf.variables[self.netcdf_var_time][:].tolist()
            self.write_idx = time_array.index(timestep)
        elif self.write_idx >= 0:
            self.write_idx += 1
        else:
            self.write_idx = 0

    def opened(self) -> bool:
        return not self.nf is None

    def close(self):
        super().close()
        if self.nf is not None:
            self.nf.close()
        self.nf = None

    def __set_property(self, netcdf_var: NetCDF4Variable, property: str, section: str, option: str, conversion_function=str):
        value = self.conf.get_config_field(section, option).strip()
        if len(value) > 0:
            setattr(netcdf_var, property, conversion_function(value))

    def __setup_netcdf_metadata(self, start_date: datetime = None):
        # General Attributes
        self.nf.history = f'Created {self.time_created}'
        self.nf.Conventions = self.NETCDF_CONVENTIONS
        self.nf.Source_Software = self.Source_Software
        self.nf.reference = self.reference
        self.nf.title = self.title
        self.nf.keywords = self.keywords
        self.nf.source = self.source
        self.nf.institution = self.institution
        self.nf.comment = self.comment

        netcdf_var_x = self.conf.get_config_field('VAR_X','NAME')
        netcdf_var_y = self.conf.get_config_field('VAR_Y','NAME')
        # netcdf_var_time = self.conf.get_config_field('VAR_TIME','NAME')

#         ncols = int(self.conf.get_config_field('DIMENSION','COLUMNS'))
#         nrows = int(self.conf.get_config_field('DIMENSION','ROWS'))

        #Dimensions
        self.nf.createDimension(netcdf_var_x, self.conf.dem_ncols)
        self.nf.createDimension(netcdf_var_y, self.conf.dem_nrows)
        self.nf.createDimension(self.netcdf_var_time, None)

        #Variables
        longitude = self.nf.createVariable(netcdf_var_x, self.NETCDF_VAR_DATA_TYPE, (netcdf_var_x, ))
        longitude.standard_name= self.conf.get_config_field('VAR_X','STANDARD_NAME')
        longitude.long_name= self.conf.get_config_field('VAR_X','LONG_NAME')
        longitude.units = self.conf.get_config_field('VAR_X','UNIT')

        latitude = self.nf.createVariable(netcdf_var_y, self.NETCDF_VAR_DATA_TYPE, (netcdf_var_y, ))
        latitude.standard_name= self.conf.get_config_field('VAR_Y','STANDARD_NAME')
        latitude.long_name= self.conf.get_config_field('VAR_Y','LONG_NAME')
        latitude.units = self.conf.get_config_field('VAR_Y','UNIT')

        time = self.nf.createVariable(self.netcdf_var_time, self.NETCDF_VAR_DATA_TYPE, (self.netcdf_var_time, ))
        time.standard_name = self.conf.get_config_field('VAR_TIME','STANDARD_NAME')
        time.frequency = self.time_frequency
        time.calendar = self.calendar_type
        time.units = self.calendar_time_unit

        proj = self.nf.createVariable(self.conf.get_config_field('PROJECTION','GRID_MAPPING'), self.NETCDF_COORDINATES_DATA_TYPE)
        self.__set_property(proj, 'grid_mapping_name', 'PROJECTION', 'GRID_MAPPING')
        self.__set_property(proj, 'false_easting', 'PROJECTION', 'FALSE_EASTING', float)
        self.__set_property(proj, 'false_northing', 'PROJECTION', 'FALSE_NORTHING', float)
        self.__set_property(proj, 'longitude_of_projection_origin', 'PROJECTION', 'ORIGIN_LONGITUDE', float)
        self.__set_property(proj, 'latitude_of_projection_origin', 'PROJECTION', 'ORIGIN_LATITUDE', float)
        self.__set_property(proj, 'semi_major_axis', 'PROJECTION', 'SEMI_MAJOR_AXIS', float)
        self.__set_property(proj, 'inverse_flattening', 'PROJECTION', 'INVERSE_FLATTENING', float)
        self.__set_property(proj, 'proj4_params', 'PROJECTION', 'PARAMS')
        self.__set_property(proj, 'EPSG_code', 'PROJECTION', 'EPSG_CODE')
        self.__set_property(proj, 'spatial_ref', 'PROJECTION', 'STRING')
        # self.__set_property(proj, 'longitude_of_prime_meridian', 'PROJECTION', 'LONGITUDE_PRIME_MERIDIAN', float)
        # self.__set_property(proj, 'GeoTransform', 'PROJECTION', 'GEO_TRANSFORM', self.__get_tuple)

        var_data_type_packed = self.conf.get_config_field('PROPERTIES', 'DATA_TYPE_PACKED')

        value = self.nf.createVariable(self.var_code, var_data_type_packed, (self.netcdf_var_time, netcdf_var_y, netcdf_var_x),
                                       zlib=True, complevel=self.COMPRESSION_LEVEL, fill_value=self.conf.VALUE_NAN) # , least_significant_digit=2)
        value.missing_value=self.conf.VALUE_NAN
        value.standard_name = self.var_standard_name
        value.long_name = self.var_long_name
        value.cell_methods = self.cell_methods
        value.units = self.units
        value.valid_min = self.conf.value_min_packed
        if self.conf.value_max != self.conf.VALUE_NAN:
            value.valid_max = self.conf.value_max_packed
        value.scale_factor = self.conf.scale_factor
        value.add_offset = self.conf.add_offset
        value.set_auto_maskandscale(True)
        self.__set_property(value, 'grid_mapping', 'PROJECTION', 'GRID_MAPPING')
        self.__set_property(value, 'esri_pe_string', 'PROJECTION', 'STRING')
        # self.__set_property(value, 'longitude_of_prime_meridian', 'PROJECTION', 'LONGITUDE_PRIME_MERIDIAN')
        # self.__set_property(value, 'GeoTransform', 'PROJECTION', 'GEO_TRANSFORM')

        latitude[:] = self.conf.dem.lat_values
        longitude[:] = self.conf.dem.lon_values

        self.print_msg(f'var-name: {self.var_standard_name}')
        self.print_msg(f'data-type-packed: {var_data_type_packed}')
        self.print_msg(f'scale: {value.scale}')
        self.print_msg(f'scale-factor: {value.scale_factor}')
        self.print_msg(f'offset: {value.add_offset}')
        self.print_msg(f'valid-min: {value.valid_min}')
        if self.conf.value_max != self.conf.VALUE_NAN:
            self.print_msg(f'valid-max: {value.valid_max}')
        else:
            self.print_msg('valid-max: N/A')


class GDALWriter(OutputWriter):

    def __init__(self, conf: Config, overwrite_file: bool = False, quiet_mode: bool = False):
        super().__init__(conf, overwrite_file, quiet_mode)
        self.current_timestamp = None
        self.driver_gtiff = gdal.GetDriverByName("GTiff")

    def setup_dataset_metadata(self, ds: gdal.Dataset) -> gdal.Dataset:
        ds.SetMetadataItem('SCALE', f'{self.conf.scale_factor}')
        ds.SetMetadataItem('OFFSET', f'{self.conf.add_offset}')
        ds.SetMetadataItem('NODATA', f'{self.conf.VALUE_NAN}')
        ds.SetMetadataItem('Source_Software', self.Source_Software)
        ds.SetMetadataItem('Created', f'{self.time_created}')
        ds.SetMetadataItem('reference', f'{self.reference}')
        ds.SetMetadataItem('title', f'{self.title}')
        ds.SetMetadataItem('keywords', f'{self.keywords}')
        ds.SetMetadataItem('source', f'{self.source}')
        ds.SetMetadataItem('institution', f'{self.institution}')
        ds.SetMetadataItem('comment', f'{self.comment}')
        ds.SetMetadataItem('var_code', f'{self.var_code}')
        ds.SetMetadataItem('var_standard_name', f'{self.var_standard_name}')
        ds.SetMetadataItem('var_long_name', f'{self.var_long_name}')
        ds.SetMetadataItem('cell_methods', f'{self.cell_methods}')
        ds.SetMetadataItem('units', f'{self.units}')
        if self.current_timestamp is not None:
            ds.SetMetadataItem('Timestamp', f'{self.current_timestamp}')
        return ds
    
    def setup_grid(self, grid: np.ndarray, print_stats: bool = True) -> np.ndarray:
        if print_stats:
            values = self.setNaN(copy.deepcopy(grid))
            self.print_grid_statistics(values)
        return grid

    def write(self, grid: np.ndarray, timestamp: datetime = None, print_stats: bool = True):
        if timestamp is not None:
            self.current_timestamp = timestamp.strftime(FileUtils.DATE_PATTERN_SEPARATED)
        else:
            self.current_timestamp = None
        cur_grid = self.setup_grid(grid, print_stats)
        self.write_timestep(cur_grid)
        self.current_timestamp = None

    def write_timestep(self, grid: np.ndarray, timestep: int = -1):
        if not self.opened():
            raise Exception("Dataset was not initialized.")
        if self.overwrite_file or not self.filepath.is_file():
            self.print_msg(f'Generating file: {self.filepath}')
            size_lats, size_lons = grid.shape
            ds = self.driver_gtiff.Create(str(self.filepath), size_lons, size_lats, 1, gdal.GDT_Int16, options=['COMPRESS=LZW'])
            # Upper Left x, East-West px resolution, rotation, Upper Left y, rotation, North-South px resolution
            ds.SetGeoTransform( [self.conf.COORDINATE_MIN_X, self.conf.CELL_SIZE_X, 0, self.conf.COORDINATE_MAX_Y, 0, self.conf.CELL_SIZE_Y ] )
            # Set CRS
            srs = osr.SpatialReference()
            srs.SetWellKnownGeogCS("WGS84")
            ds.SetProjection( srs.ExportToWkt() )
            # Write the band
            ds.GetRasterBand(1).SetNoDataValue(self.conf.VALUE_NAN)
            ds.GetRasterBand(1).SetScale(self.conf.scale_factor)
            ds.GetRasterBand(1).SetOffset(self.conf.add_offset)
            ds = self.setup_dataset_metadata(ds)
            ds.GetRasterBand(1).WriteArray(grid.astype(np.int16))
            ds.FlushCache()
            ds = None
            self.print_msg(f'Wrote file: {self.filepath}')
        else:
            self.print_msg(f'File already exists: {self.filepath}')
        self.close()

