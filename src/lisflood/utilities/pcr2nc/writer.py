"""
Copyright 2019 European Union

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

import numpy as np
from netCDF4 import Dataset

from lisflood.utilities.pcr2nc.reader import PCRasterReader


class NetCDFWriter:
    """
    This class manages all aspects concerning definition and writing of a NetCDF4 file.
    """

    FORMATS = {'classic': 'NETCDF4_CLASSIC', 'normal': 'NETCDF4'}
    DATUM = {
        'ETRS89': 'PROJCS["JRC_LAEA_ETRS-DEF",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],UNIT["Meter",1.0]]',
        'WGS84': 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]]',
        'GISCO': 'PROJCS["PCS_Lambert_Azimuthal_Equal_Area",GEOGCS["GCS_User_Defined",DATUM["D_User_Defined",SPHEROID["User_Defined_Spheroid",6378388.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",9.0],PARAMETER["Latitude_Of_Origin",48.0],UNIT["Meter",1.0]]',
    }

    def __init__(self, filename, nc_metadata, pcr_metadata, mapstack=True):
        self.name = '{}.nc'.format(filename) if not filename.endswith('.nc') else filename
        self.nc_metadata = nc_metadata
        self.hour = self.nc_metadata.get('time', {}).get('hour', '00')

        self.pcr_metadata = pcr_metadata
        self.is_mapstack = mapstack
        self.time, self.variable = self._init_dataset()

        self.hour_timestep = float(self.hour) / 24
        self.values = []
        self.timesteps = []
        self.current_count = 0
        self.current_idx1 = 0
        self.current_idx2 = 0

    def _init_dataset(self):
        # frmt = 'normal' if np.issubdtype(self.pcr_metadata['dtype'], np.unsignedinteger) else 'classic'
        self.frmt = self.nc_metadata.get('format', 'NETCDF4')
        self.nf = Dataset(self.name, 'w', format=self.frmt)
        self.nf.history = 'Created {}'.format(time.ctime(time.time()))
        self.nf.Conventions = 'CF-1.7'
        self.nf.Source_Software = 'JRC.E1 lisflood nexus wefe - pcr2nc'
        self.nf.source = self.nc_metadata.get('source')
        self.nf.reference = self.nc_metadata.get('reference')

        # Dimensions
        if self.is_mapstack:
            self.nf.createDimension('time', None)
        self.nf.createDimension('yc', self.pcr_metadata['rows'])
        self.nf.createDimension('xc', self.pcr_metadata['cols'])

        # define coordinates variables by calling one of the define_* functions
        datum_function = 'define_{}'.format(self.nc_metadata['geographical']['datum'].lower())
        post_datum_function = '{}_post'.format(datum_function)
        getattr(self, datum_function)()

        time_nc = None
        vardimensions = ('yc', 'xc')
        if self.is_mapstack:
            # time variable
            time_units = self.nc_metadata['time'].get('units', '')
            time_nc = self.nf.createVariable('time', 'f8', ('time',))
            time_nc.standard_name = 'time'
            if str(self.hour) != '24':
                time_nc.units = '{} {}:00'.format(time_units, str(self.hour).zfill(2))
            else:
                # observation is at 24h...need to rotate one day more
                start_date = datetime.datetime.strptime(time_units[-10:], '%Y-%m-%d')  # 'days since 1996-01-01'
                start_date = start_date + datetime.timedelta(days=1)
                time_nc.units = 'days since {} 00:00'.format(start_date.strftime('%Y-%m-%d'))
            time_nc.calendar = self.nc_metadata['time'].get('calendar', 'proleptic_gregorian')
            vardimensions = ('time', 'yc', 'xc')

        # data variable
        complevel = self.nc_metadata['variable'].get('compression')
        additional_args = {'zlib': bool(complevel)}
        if complevel:
            print('Applying compression level', str(complevel))
            additional_args['complevel'] = complevel
        if np.issubdtype(self.pcr_metadata['dtype'], np.floating):
            additional_args['least_significant_digit'] = self.nc_metadata.get('least_significant_digit', None)

        values_nc = self.nf.createVariable(self.nc_metadata['variable'].get('shortname', ''),
                                           self.pcr_metadata['dtype'],
                                           vardimensions, **additional_args)
        getattr(self, post_datum_function)(values_nc)

        values_nc.standard_name = self.nc_metadata['variable'].get('shortname', '')
        values_nc.long_name = self.nc_metadata['variable'].get('longname', '')
        values_nc.units = self.nc_metadata['variable'].get('units', '')
        return time_nc, values_nc

    def add_to_stack(self, pcr_map, time_step=None):
        """
        Add a PCRaster map to the NetCDF4 file.
        :param time_step: int, it's basically the extension of pcraster map file
            For single files (ie not time series) time_step is None
        :param pcr_map: PCRasterMap object
        """
        print('Adding', pcr_map.filename, 'timestep', str(time_step), 'hour', self.hour_timestep)
        values = pcr_map.data
        if not np.issubdtype(values.dtype, np.integer):
            values[values == pcr_map.mv] = np.nan
        self.values.append(values)
        if time_step:
            self.timesteps.append(float(time_step))
        self.current_count += 1

        if self.current_count == 20:
            self.current_idx2 += self.current_count
            print('Writing a chunk...')
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

    def finalize(self):
        """
        Write last maps to the stack and close the NetCDF4 dataset.
        """
        print('Writing...', self.name)
        if self.is_mapstack:
            print('Writing time dimension...', self.timesteps[:4], '...', self.timesteps[-4:])
            self.time[:] = np.array(self.timesteps, dtype=np.float64)

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
        print('Defining WGS84 coordinates variables')
        longitude = self.nf.createVariable('lon', 'f8', ('xc',))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'longitude coordinate'
        longitude.units = 'degrees_east'

        latitude = self.nf.createVariable('lat', 'f8', ('yc',))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'latitude coordinate'
        latitude.units = 'degrees_north'
        longitude[:] = self.pcr_metadata['lons']
        latitude[:] = self.pcr_metadata['lats']

    def define_wgs84_post(self, values_var):
        values_var.coordinates = 'lon lat'
        values_var.esri_pe_string = self.DATUM.get(self.nc_metadata['geographical'].get('datum', 'WGS84').upper(), '')

    def define_etrs89(self):
        """
        Define a ETRS89 reference system
        """
        print('Defining ETRS89 coordinates variables')
        # Variables
        x = self.nf.createVariable('x', 'f8', ('xc',))
        y = self.nf.createVariable('y', 'f8', ('yc',))
        x.standard_name = 'projection_x_coordinate'
        x.long_name = 'x coordinate of projection'
        x.units = 'Meter'

        y.standard_name = 'projection_y_coordinate'
        y.long_name = 'y coordinate of projection'
        y.units = 'Meter'
        x[:] = self.pcr_metadata['lons']
        y[:] = self.pcr_metadata['lats']

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
        values_var.esri_pe_string = self.DATUM.get(self.nc_metadata['geographical'].get('datum', 'WGS84').upper(), '')

    def define_gisco(self):
        """
        It defines a custom LAEA ETRS89 GISCO reference system
        """
        print('Defining GISCO coordinates variables')
        # Variables
        x = self.nf.createVariable('x', 'f8', ('xc',))
        y = self.nf.createVariable('y', 'f8', ('yc',))
        x.standard_name = 'projection_x_coordinate'
        x.long_name = 'x coordinate of projection'
        x.units = 'Meter'

        y.standard_name = 'projection_y_coordinate'
        y.long_name = 'y coordinate of projection'
        y.units = 'Meter'
        x[:] = self.pcr_metadata['lons']
        y[:] = self.pcr_metadata['lats']

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
        values_var.esri_pe_string = self.DATUM.get(self.nc_metadata['geographical'].get('datum', 'WGS84').upper(), '')


def convert(config):
    input_set = config['input_set']
    reader = PCRasterReader(input_set)
    pcr_metadata = reader.get_metadata_from_set()
    writer = NetCDFWriter(config.get('output_filename') or config.get('variable'),
                          config['metadata'],
                          pcr_metadata,
                          mapstack=not reader.input_is_single_file())
    for pcr_map, time_step in reader.fileset:
        writer.add_to_stack(pcr_map, time_step)
        pcr_map.close()
    writer.finalize()
