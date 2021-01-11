"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

Module containing the code for the NetCDFWriter class
"""

import datetime
import time

from nine import IS_PYTHON2
if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from ..readers import PCRasterMap
from .. import logger

wgs84_lons = np.arange(-180, 180, 0.1)
wgs84_lats = np.arange(-60, 90, 0.1)


class NetCDFWriter:
    """
    This class manages all aspects concerning definition and writing of a NetCDF4 file.
    """

    WKT_STRINGS = {
        'ETRS89': 'PROJCS["JRC_LAEA_ETRS-DEF",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]',
        'WGS84': 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]]',
        'GISCO': 'PROJCS["PCS_Lambert_Azimuthal_Equal_Area",GEOGCS["GCS_User_Defined",DATUM["D_User_Defined",SPHEROID["User_Defined_Spheroid",6378388.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",9.0],PARAMETER["Latitude_Of_Origin",48.0],UNIT["Meter",1.0]]',
    }

    def __init__(self, filename, is_mapstack=True, **metadata):
        """
        :param filename: Path to netCDF file to create. If file exists, open it in append mode
        :param is_mapstack: True if output file is a mapstack
        :param metadata: metadata dict
            format: NETCDF4 or NETCDF_CLASSIC
            time:
                calendar: netcdf calendar (e.g. proleptic_gregorian)
                units: time unit (days since 1999-01-01)
            variable: dictionary
              shortname: name of variable
              description: description
              longname: long name
              units: variable units
              least_significant_digit: (to use in conjuction with compression:
                    power of ten of the smallest decimal place in the data that is a reliable value.
                    For example if the data has a precision of 0.1, then setting least_significant_digit=1
                    will cause data the data to be quantized using numpy.around(scale*data)/scale,
                    where scale = 2**bits, and bits is determined so that a precision of 0.1 is retained
                    (in this case bits=4))
              compression: compression level [1-9]
            source: Organisation producing the file (ECMWF)
            reference: Organisation reference (JRC E1)
            geographical: dictionary
                datum: [WGS84|ETRS89|GISCO]
            dtype: internal numpy type (e.g. float64, int8 etc.)
            rows: number of rows in 2D values array
            cols: number of cols in 2D values array
            lats: 1D latitude array
            lons: 1D longitude array
        """
        self.name = '{}.nc'.format(filename) if not filename.endswith('.nc') else filename
        self.path = Path(self.name)
        self.append = self.path.exists()
        self.metadata = metadata
        self.is_mapstack = is_mapstack
        self.hour = float(self.metadata.get('time', {}).get('hour') or 0)

        # you can pass the MV to set in netcdf files directly in yaml configuration, otherwuse np.nan is used
        self.mv = self.metadata['variable'].get('mv')
        if self.mv is not None:
            self.mv = int(self.mv) if np.issubdtype(self.metadata['dtype'], np.integer) else float(self.mv)
        else:
            self.mv = np.nan
        self.current_idx1 = 0
        self.current_idx2 = 0
        self.time, self.variable = self._init_dataset()

        self.hour_timestep = self.hour / 24.
        self.values = []
        self.timesteps = []
        self.current_count = 0

    def _init_dataset(self):
        """
        Create dimensions and variables

        :return: (time, values) netCDF4 variables in a tuple
        """
        self.frmt = self.metadata.get('format', 'NETCDF4')
        time_nc = None
        if not self.append:
            self.nf = Dataset(self.name, 'w', format=self.frmt)
            self.nf.history = 'Created {}'.format(time.ctime(time.time()))
            self.nf.Conventions = 'CF-1.7'
            self.nf.Source_Software = 'JRC.E1 lisfloodutilities - pcr2nc'
            self.nf.source = self.metadata.get('source', '')
            self.nf.reference = self.metadata.get('reference', '')

            # Dimensions
            if self.is_mapstack:
                self.nf.createDimension('time', None)
            self.nf.createDimension('lat', self.metadata.get('rows') or self.metadata['lats'].size)
            self.nf.createDimension('lon', self.metadata.get('cols') or self.metadata['lons'].size)

            # define coordinates variables by calling one of the define_* functions
            datum = self.metadata.get('geographical', {}).get('datum', '').lower()
            if not datum:
                datum = self._guess_datum(self.metadata['lats'], self.metadata['lons']) if 'lats' in self.metadata else 'wgs84'
            datum_function = 'define_{}'.format(datum)
            post_datum_function = '{}_post'.format(datum_function)
            getattr(self, datum_function)()  # call the function

            vardimensions = ('lat', 'lon')
            if self.is_mapstack:
                # time variable
                time_units = self.metadata['time'].get('units', '')
                time_nc = self.nf.createVariable('time', 'f8', ('time',))
                time_nc.standard_name = 'time'
                if time_units:
                    if str(self.hour) != '24':
                        time_nc.units = '{} {}:00'.format(time_units, str(self.hour).zfill(2))
                    else:
                        # observation is at 24h...need to rotate one day more
                        start_date = datetime.datetime.strptime(time_units[-10:], '%Y-%m-%d')  # 'days since 1996-01-01'
                        start_date = start_date + datetime.timedelta(days=1)
                        time_nc.units = 'days since {} 00:00'.format(start_date.strftime('%Y-%m-%d'))
                    time_nc.calendar = self.metadata['time'].get('calendar', 'proleptic_gregorian') if 'time' in self.metadata else 'proleptic_gregorian'
                vardimensions = ('time', 'lat', 'lon')

            # data variable
            complevel = self.metadata['variable'].get('compression')
            additional_args = {'zlib': bool(complevel)}
            if complevel:
                logger.info('Applying compression level %s', str(complevel))
                additional_args['complevel'] = complevel
                if np.issubdtype(self.metadata['dtype'], np.floating):
                    additional_args['least_significant_digit'] = self.metadata.get('least_significant_digit')

            values_nc = self.nf.createVariable(self.metadata['variable'].get('shortname', ''),
                                               self.metadata['dtype'], vardimensions,
                                               fill_value=self.mv, **additional_args)
            getattr(self, post_datum_function)(values_nc)

            values_nc.standard_name = self.metadata['variable'].get('shortname', '')
            values_nc.long_name = self.metadata['variable'].get('longname', '')
            values_nc.units = self.metadata['variable'].get('units', '')
        else:
            # open in append mode
            self.nf = Dataset(self.name, 'a', format=self.frmt)
            values_nc = self.nf.variables[self.metadata['variable'].get('shortname', '')]
            if self.is_mapstack:
                time_nc = self.nf.variables['time']
                self.current_idx1 = self.current_idx2 = time_nc.size
        return time_nc, values_nc

    def add_to_stack(self, amap, time_step=None):
        """
        Add a PCRaster map or numpy 2D array to the NetCDF4 file.
        :param time_step: int, it's basically the extension of pcraster map file
            For single files (ie not time series) time_step is None
        :param amap: numpy.ndarray or PCRasterMap object
        """
        if isinstance(amap, PCRasterMap):
            # PCRasterMap
            logger.info('Adding %s - timestep %s - hour %s', amap.filename, str(time_step), self.hour_timestep)
            values = amap.data
            values[values == amap.mv] = self.mv
        else:
            # amap is a simple numpy array
            logger.info('Adding array - timestep %s - hour %s', str(time_step), self.hour_timestep)
            values = amap
            values[values == np.nan] = self.mv
        self.values.append(values)

        if time_step is not None:
            self.timesteps.append(float(time_step))

        self.current_count += 1

        if self.current_count == 20:
            self.current_idx2 += self.current_count
            logger.info('Writing a chunk into output file...')
            dtype = self.values[0].dtype
            if self.is_mapstack:
                self.variable[self.current_idx1:self.current_idx2, :, :] = np.array(self.values, dtype=dtype)
            else:
                self.variable[:, :] = np.array(self.values, dtype=dtype)
            # update slicing indexes
            self.current_idx1 = self.current_idx2
            # reset
            self.values = []
            self.current_count = 0

    def finalize(self, timesteps=None):
        """
        Write last maps to the stack and close the NetCDF4 dataset.
        """
        logger.info('Writing %s', self.name)
        if self.is_mapstack:
            if timesteps is None:
                timesteps = self.timesteps
            self.time[:] = np.array(timesteps, dtype=np.float64) if not isinstance(timesteps, np.ndarray) else timesteps

        if self.values:
            dtype = self.values[0].dtype
            if self.is_mapstack:
                self.current_idx2 += self.current_count
                self.variable[self.current_idx1:self.current_idx2, :, :] = np.array(self.values, dtype=dtype)
            else:
                self.variable[:, :] = np.array(self.values, dtype=dtype)
        self.nf.close()

    def define_wgs84(self):
        """
        Define WGS84 reference system
        """
        # coordinates variables
        logger.info('Defining WGS84 coordinates variables')
        longitude = self.nf.createVariable('lon', 'f8', ('lon',))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'longitude coordinate'
        longitude.units = 'degrees_east'

        latitude = self.nf.createVariable('lat', 'f8', ('lat',))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'latitude coordinate'
        latitude.units = 'degrees_north'
        if 'lons' in self.metadata:
            longitude[:] = self.metadata['lons']
            latitude[:] = self.metadata['lats']
        else:
            latitude[:] = wgs84_lons
            latitude[:] = wgs84_lats

    def define_wgs84_post(self, values_var):
        values_var.coordinates = 'lon lat'
        values_var.esri_pe_string = self.WKT_STRINGS.get(self.metadata.get('geographical', {}).get('datum', 'WGS84').upper(), '')

    def define_etrs89(self):
        """
        Define a ETRS89 reference system
        """
        logger.info('Defining ETRS89 coordinates variables')
        x = self.nf.createVariable('x', 'f8', ('lon',))
        y = self.nf.createVariable('y', 'f8', ('lat',))
        x.standard_name = 'projection_x_coordinate'
        x.long_name = 'x coordinate of projection'
        x.units = 'Meter'

        y.standard_name = 'projection_y_coordinate'
        y.long_name = 'y coordinate of projection'
        y.units = 'Meter'
        x[:] = self.metadata['lons']
        y[:] = self.metadata['lats']

        proj = self.nf.createVariable('laea', 'i4')
        proj.grid_mapping_name = 'lambert_azimuthal_equal_area'
        proj.false_easting = 4321000.0
        proj.false_northing = 3210000.0
        proj.longitude_of_projection_origin = 10.0
        proj.latitude_of_projection_origin = 52.0
        proj.semi_major_axis = 6378137.0
        proj.inverse_flattening = 298.257223563
        proj.proj4_params = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
        proj.EPSG_code = "EPSG:3035"

    def define_etrs89_post(self, values_var):
        values_var.coordinates = 'x y'
        values_var.grid_mapping = 'lambert_azimuthal_equal_area'
        values_var.esri_pe_string = self.WKT_STRINGS.get(self.metadata.get('geographical', {}).get('datum', 'ETRS89').upper(), '')

    def define_gisco(self):
        """
        It defines a custom LAEA ETRS89 GISCO reference system
        """
        logger.info('Defining GISCO coordinates variables')
        x = self.nf.createVariable('x', 'f8', ('lon',))
        y = self.nf.createVariable('y', 'f8', ('lat',))
        x.standard_name = 'projection_x_coordinate'
        x.long_name = 'x coordinate of projection'
        x.units = 'Meter'

        y.standard_name = 'projection_y_coordinate'
        y.long_name = 'y coordinate of projection'
        y.units = 'Meter'
        x[:] = self.metadata['lons']
        y[:] = self.metadata['lats']

        proj = self.nf.createVariable('laea', 'i4')
        proj.grid_mapping_name = 'lambert_azimuthal_equal_area'
        proj.false_easting = 0.0
        proj.false_northing = 0.0
        proj.longitude_of_projection_origin = 9.0
        proj.latitude_of_projection_origin = 48.0
        proj.semi_major_axis = 6378388.0
        proj.inverse_flattening = 0.0
        proj.proj4_params = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
        proj.EPSG_code = "EPSG:3035"

    def define_gisco_post(self, values_var):
        values_var.coordinates = 'x y'
        values_var.grid_mapping = 'lambert_azimuthal_equal_area'
        values_var.esri_pe_string = self.WKT_STRINGS.get(self.metadata['geographical'].get('datum', 'WGS84').upper(), '')

    def _guess_datum(self, lats, lons):
        # very naive version that works for etrs89 vs wgs84
        int_coords = abs(int(lats[1] - lats[0])) == abs(float(lats[1] - lats[0]))
        return 'etrs89' if int_coords else 'wgs84'
