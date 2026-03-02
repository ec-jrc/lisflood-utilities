import os

import pytest


from lisfloodutilities.waterregions.define_waterregions import define_waterregions, parse_metadata
from lisfloodutilities.waterregions.verify_waterregions import verify_waterregions

from pathlib import Path
from netCDF4 import Dataset
import numpy as np

cur_folder = Path(os.path.dirname(os.path.realpath(__file__)), 'data/waterregions')

class TestWaterRegions:

    # @pytest.mark.defwaterregions
    def test_define_waterregions(self):
        calib_points = Path(cur_folder, 'calib_points_test.txt')
        countries_id = Path(cur_folder, 'countries_id_test.nc')
        ldd = Path(cur_folder, 'ldd_test.nc')
        waterregion_init = Path(cur_folder, 'waterregions_initial_test.nc')
        output_wr = Path(cur_folder, 'waterregions_output_test.nc')
        reference_output = Path(cur_folder, 'waterregions_reference_output.nc')
        temp_dir = cur_folder
        metadata_path = Path(cur_folder, 'metadata_test.yaml')
        metadata = parse_metadata(metadata_path)

        # Clean up any existing output file before running the test
        if output_wr.exists():
            os.unlink(output_wr.as_posix())

        define_waterregions(calib_points, countries_id, ldd, waterregion_init,
                            output_wr, tmpdir=temp_dir, metadata_parsed=metadata)

        # Use context managers to ensure files are closed properly
        with Dataset(output_wr) as out_ds, Dataset(reference_output) as ref_ds:
            out_waterregions = out_ds.variables['wregion'][:]
            ref_waterregions = ref_ds.variables['wregion'][:]

            # Verify shapes match
            assert out_waterregions.shape == ref_waterregions.shape, (
                f"Shape mismatch: {out_waterregions.shape} vs {ref_waterregions.shape}"
            )

            # Verify data equality with detailed diff on failure
            np.testing.assert_array_equal(
                out_waterregions,
                ref_waterregions,
                err_msg="Water‑region data differ between output and reference"
            )
        os.unlink(output_wr.as_posix())

    def test_verify_waterregions(self):
        calib_catchments = Path(cur_folder, 'calib_catchments_test.nc')
        waterregions = Path(cur_folder, 'waterregions_test.nc')
        output_message = verify_waterregions(calib_catchments, waterregions)
        assert output_message == 'PASSED: Each calibration catchment contains only a finite number of water regions.'

    def test_verify_waterregions_fail(self):
        calib_catchments = Path(cur_folder, 'calib_catchments_test.nc')
        waterregions = Path(cur_folder, 'waterregions_test_NOTc.nc')
        output_message = verify_waterregions(calib_catchments, waterregions)
        assert output_message != 'PASSED: Each calibration catchment contains only a finite number of water regions.'

test_define_waterregions = TestWaterRegions().test_define_waterregions()
