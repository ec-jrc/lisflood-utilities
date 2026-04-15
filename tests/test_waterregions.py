import os

import pytest

from lisfloodutilities.readers import PCRasterMap

from lisfloodutilities.waterregions.define_waterregions import define_waterregions
from lisfloodutilities.waterregions.verify_waterregions import verify_waterregions

import pcraster as pcr

class TestWaterRegions:

    @pytest.mark.defwaterregions
    def test_define_waterregions(self):
        calib_points = 'tests/data/waterregions/calib_points_test.txt'
        countries_id = 'tests/data/waterregions/countries_id_test.map'
        ldd = 'tests/data/waterregions/ldd_test.map'
        waterregion_init = 'tests/data/waterregions/waterregions_initial_test.map'
        output_wr = 'tests/data/waterregions/waterregions_output_test.map'
        reference_output = 'tests/data/waterregions/waterregions_reference_output.map'

        define_waterregions(calib_points, countries_id, ldd, waterregion_init, output_wr, tmpdir='tests/data/waterregions/')

        out = PCRasterMap(output_wr)
        reference = PCRasterMap(reference_output)

        assert out == reference

        out.close()
        reference.close()
        os.unlink(output_wr)
        
    def test_verify_waterregions(self):
        calib_catchments = 'tests/data/waterregions/calib_catchments_test.nc'
        waterregions = 'tests/data/waterregions/waterregions_test.nc'
        output_message = verify_waterregions(calib_catchments, waterregions)
        assert output_message == 'PASSED: Each calibration catchment contains only a finite number of water regions.'

    def test_verify_waterregions_fail(self):
        calib_catchments = 'tests/data/waterregions/calib_catchments_test.nc'
        waterregions = 'tests/data/waterregions/waterregions_test_NOTc.nc'
        output_message = verify_waterregions(calib_catchments, waterregions)
        assert output_message != 'PASSED: Each calibration catchment contains only a finite number of water regions.'
