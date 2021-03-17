"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import os

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities.readers import PCRasterMap
from lisfloodutilities.pcr2nc import convert

from . import TestWithCleaner


class TestPcr2nc(TestWithCleaner):
    def test_convert(self):
        dataset = 'tests/data/folder_d'
        out = 'tests/data/pcr2nc_test.nc'
        self.cleanups.append((os.unlink, (out,)))

        metadata = {
            'variable': {
                'shortname': 'x',
                'units': 'm',
                'least_significant_digit': 2,
                'compression': 9
            },
            'geographical': {
                'datum': 'WGS84'
            },
            'source': 'JRC E1',
            'reference': 'JRC E1',
            'time': {
                'calendar': 'proleptic_gregorian',
                'units': 'days since 1999-01-01'
            }
        }
        map_0 = PCRasterMap('tests/data/folder_d/map.001')
        map_1 = PCRasterMap('tests/data/folder_d/map.002')
        map_2 = PCRasterMap('tests/data/folder_d/map.003')
        map_3 = PCRasterMap('tests/data/folder_d/map.004')
        convert(dataset, out, metadata)
        with Dataset(out) as nc:
            time_arr = nc.variables['time'][:]
            lat_arr = nc.variables['lat'][:]
            lon_arr = nc.variables['lon'][:]
            assert time_arr.size == 4
            assert (lat_arr.size, lon_arr.size) == (35, 35)
            var_0 = nc.variables['x'][0, :, :]
            var_1 = nc.variables['x'][1, :, :]
            var_2 = nc.variables['x'][2, :, :]
            var_3 = nc.variables['x'][3, :, :]
            assert np.allclose(map_0.data, var_0)
            assert np.allclose(map_1.data, var_1)
            assert np.allclose(map_2.data, var_2)
            assert np.allclose(map_3.data, var_3)
        map_0.close()
        map_1.close()
        map_2.close()
        map_3.close()
