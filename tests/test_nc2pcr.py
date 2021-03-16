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

from lisfloodutilities.readers import NetCDFMap, PCRasterMap
from lisfloodutilities.nc2pcr import convert
from lisfloodutilities.compare.pcr import PCRComparator


class TestNc2Pcr:
    def test_convert(self):
        dataset = 'tests/data/cutmaps/ldd_eu.nc'
        out = 'tests/data/nc2pcr_test.map'
        clonemap = 'tests/data/cutmaps/area_eu.map'

        convert(dataset, out, clonemap=clonemap, is_ldd=True)
        assert os.path.exists(out)

        mappcr_0 = PCRasterMap(out)
        map_0 = NetCDFMap(dataset)

        variable = {n: v for n, v in map_0.data}['ec_ldd_repaired']
        # test non masked values
        values = np.ma.masked_where(variable.values == map_0.mv, variable.values, copy=False)
        mask = np.ma.getmask(values)
        values_pcr = np.ma.masked_where(mask, mappcr_0.data, copy=False)
        assert values.shape == values_pcr.shape
        assert np.ma.allclose(values, values_pcr, masked_equal=True)
        os.unlink(out)
        mappcr_0.close()
        map_0.close()

    def test_convert_noclonemap(self):
        dataset = 'tests/data/cutmaps/ldd_eu.nc'
        out = 'tests/data/nc2pcr_test.map'

        convert(dataset, out, is_ldd=True)
        assert os.path.exists(out)

        mappcr_0 = PCRasterMap(out)
        map_0 = NetCDFMap(dataset)

        variable = {n: v for n, v in map_0.data}['ec_ldd_repaired']
        # test non masked values
        values = np.ma.masked_where(variable.values == map_0.mv, variable.values, copy=False)
        mask = np.ma.getmask(values)
        values_pcr = np.ma.masked_where(mask, mappcr_0.data, copy=False)
        assert values.shape == values_pcr.shape
        assert np.ma.allclose(values, values_pcr, masked_equal=True)
        os.unlink(out)
        mappcr_0.close()
        map_0.close()

    def test_convert_noclone_same_values(self):

        dataset = 'tests/data/cutmaps/ldd_eu.nc'
        out1 = 'tests/data/nc2pcr_test1.map'  # produced without clonemap
        out2 = 'tests/data/nc2pcr_test2.map'  # produced with clonemap
        clonemap = 'tests/data/cutmaps/area_eu.map'

        convert(dataset, out1, is_ldd=True)
        convert(dataset, out2, clonemap=clonemap, is_ldd=True)

        assert os.path.exists(out1)
        assert os.path.exists(out2)
        # os.unlink(out1)
        # os.unlink(out2)

        compare = PCRComparator(array_equal=True, for_testing=True)
        compare.compare_files(out1, out2)
