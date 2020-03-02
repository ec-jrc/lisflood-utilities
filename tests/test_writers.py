import numpy as np
from netCDF4 import Dataset

from lisfloodutilities.writers import NetCDFWriter


class TestNetcdfWriter:

    def test_simple(self):
        filename = 'test.nc'
        step = 0.1
        start_x = -180
        end_x = 180
        start_y = -60
        end_y = 90
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
            # 'rows': 100,
            # 'cols': 100,
            'lons': np.arange(start_x, end_x, step),
            'lats': np.arange(start_y, end_y, step),
        }
        cols, rows = metadata['lats'].size, metadata['lons'].size
        w = NetCDFWriter(filename, is_mapstack=False, **metadata)
        a = np.random.uniform(size=(cols, rows))
        w.add_to_stack(a, None)
        w.finalize()
        with Dataset(filename) as f:
            assert np.allclose(f.variables['x'][:, :], a)
            assert f.variables['x'].units == metadata['variable']['units']
