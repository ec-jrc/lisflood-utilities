import numpy as np
from osgeo.gdal import Band

from lisfloodutilities.readers import PCRasterMap, PCRasterReader


class TestPCRasterMap:

    def test_init(self):
        map = PCRasterMap('tests/data/folder_b/1.map')
        assert map.filename == 'tests/data/folder_b/1.map'
        assert map.cols == 35
        assert map.rows == 35
        assert map.min < 0.33
        assert map.max > 0.378

    def test_eq(self):
        map_a = PCRasterMap('tests/data/folder_a/1.map')
        map_b = PCRasterMap('tests/data/folder_a/1.map')
        map_c = PCRasterMap('tests/data/folder_a/2.map')
        assert map_a == map_b
        assert map_a != map_c

    def test_data(self):
        map_a = PCRasterMap('tests/data/folder_a/1.map')
        map_b = PCRasterMap('tests/data/folder_c/1.map')
        assert map_a.data.shape == (35, 35)
        assert map_a.mv == -3.4028234663852886e+38
        assert np.allclose(map_a.data, map_b.data)
        assert map_a.data[15, 25] == -3.4028234663852886e+38
        assert map_a.geo_transform == (5450000.0, 5000.0, 0.0, 2615000.0, 0.0, -5000.0)
        assert isinstance(map_a.band, Band)
        assert map_a.lats.shape == (35, )
        assert map_a.lons.shape == (35,)
        assert map_a.lats[0] == 2612500.
        assert map_a.lats[-1] == 2442500.
        map_a.close()
        assert map_a.band is None
        assert map_b.band is not None and isinstance(map_b.band, Band)


class TestPCRasterReader:

    def test_input_singlefile(self):
        pass
    def test_input_directory(self):
        pass
    def test_input_wildcard(self):
        pass
