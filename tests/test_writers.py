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
from numpy import ma

from netCDF4 import Dataset

from lisfloodutilities.writers import NetCDFWriter, PCRasterWriter
from lisfloodutilities.readers import PCRasterMap


class TestPCRasterWriter:
    filename = 'tests/data/test.map'
    clonemap_filename = 'tests/data/asia.map'

    def test_simple(self):
        writer = PCRasterWriter(clonemap=self.clonemap_filename)
        clonemap = PCRasterMap(self.clonemap_filename)
        values = np.random.random(clonemap.data.shape)
        writer.write(self.filename, values)
        writer.close()

        outmap = PCRasterMap(self.filename)
        values_pcr = ma.masked_where(writer._mask, outmap.data, copy=False)
        values_pcr = ma.filled(values_pcr, outmap.mv)

        values = ma.masked_where(writer._mask, values, copy=False)
        values = ma.filled(values, writer.mv)

        assert outmap.mv == writer.mv
        assert values_pcr.shape == values.shape
        assert values_pcr.shape == (writer.rows, writer.cols)
        assert np.ma.allclose(values_pcr, values, masked_equal=True)
        os.unlink(self.filename)


class TestNetcdfWriter:
    filename = 'tests/data/test.nc'
    step = 0.1
    start_x = -180
    end_x = 180
    start_y = -60
    end_y = 90
    lons = np.arange(start_x, end_x, step)
    lats = np.arange(start_y, end_y, step)
    cols, rows = lats.size, lons.size

    def test_simple(self):

        metadata = {
            'variable': {
                'shortname': 'x',
                'units': 'm'
            },
            'geographical': {
                'datum': 'WGS84'
            },
            'source': 'JRC E1',
            'reference': 'JRC E1',
            'dtype': np.float,
            'lons': self.lons,
            'lats': self.lats
        }

        w = NetCDFWriter(self.filename, is_mapstack=False, **metadata)
        a = np.random.uniform(size=(self.cols, self.rows))
        w.add_to_stack(a, None)
        w.finalize()
        with Dataset(self.filename) as f:
            assert np.allclose(f.variables['x'][:, :], a)
            assert f.variables['x'].units == metadata['variable']['units']
        os.unlink(self.filename)

    def test_mapstack(self):
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
            'dtype': np.float,
            'lons': self.lons,
            'lats': self.lats,
            'time': {
                'calendar': 'proleptic_gregorian',
                'units': 'days since 1999-01-01'
            }
        }
        w = NetCDFWriter(self.filename, is_mapstack=True, **metadata)
        values_to_compare = []
        for i in range(1, 10):
            a = np.random.uniform(size=(self.cols, self.rows))
            w.add_to_stack(a, i)
            values_to_compare.append(a)
        w.finalize()

        with Dataset(self.filename) as f:
            for i, values in enumerate(values_to_compare):
                nc_values = f.variables['x'][i, :, :]
                assert nc_values.shape == values.shape
                assert np.allclose(nc_values, values)
        os.unlink(self.filename)
