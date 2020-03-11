import os

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities.writers import NetCDFWriter


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
