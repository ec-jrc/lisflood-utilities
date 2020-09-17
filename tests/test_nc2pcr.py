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


class TestNc2Pcr:
    def test_convert(self):
        dataset = 'tests/data/cutmaps/ldd.nc'
        out = 'tests/data/nc2pcr_test.map'
        outfilename = 'tests/data/nc2pcr_test_ec_ldd_repaired.map'
        clonemap = 'tests/data/cutmaps/areawithdeadsea.map'

        convert(dataset, clonemap, out)
        assert os.path.exists(outfilename)

        mappcr_0 = PCRasterMap(outfilename)
        map_0 = NetCDFMap(dataset)

        variable = {n: v for n, v in map_0.data}['ec_ldd_repaired']
        values = variable.values
        values = np.where(values == map_0.mv, mappcr_0.mv, values)
        assert values.shape == mappcr_0.data.shape
        assert np.allclose(values, mappcr_0.data, equal_nan=True)
        # assert map_0.mv == mappcr_0.mv
        os.unlink(outfilename)
        mappcr_0.close()
        map_0.close()
