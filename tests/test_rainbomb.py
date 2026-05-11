import os

import pytest

from lisfloodutilities.rainbomb.rainbomb_correction import main

from pathlib import Path
import xarray as xr
import numpy as np

program_folder = os.path.dirname(os.path.realpath(__file__))
cur_folder = Path(program_folder, 'data/rainbomb')

class TestRainBomb:

    def test_define_rainbomb(self):
        # input_file = Path(cur_folder, 'input.grb')
        output_file = Path(cur_folder, 'output.grb')
        reference_output = Path(cur_folder, 'output_reference.grb')

        argv = [f'{program_folder}/../src/lisfloodutilities/rainbomb/rainbomb_correction.py',
                '--parent_dir', f'{cur_folder}/config',
                '--input_file', f'{cur_folder}/input.grb',
                '--output_file', f'{cur_folder}/output.grb',
                '--set-grib-date']

        # Clean up any existing output file before running the test
        if output_file.exists():
            os.unlink(output_file.as_posix())

        # Run the rainbomb correction script
        main(argv)

        # Use xarray with cfgrib engine to read GRIB files
        out_ds = xr.open_dataset(output_file, engine='cfgrib')
        ref_ds = xr.open_dataset(reference_output, engine='cfgrib')

        # Get the data variable (GRIB files may have different variable names)
        out_data = out_ds[list(out_ds.data_vars)[0]].values
        ref_data = ref_ds[list(ref_ds.data_vars)[0]].values

        # Close datasets
        out_ds.close()
        ref_ds.close()

        # Verify shapes match
        assert out_data.shape == ref_data.shape, (
            f"Shape mismatch: {out_data.shape} vs {ref_data.shape}"
        )

        # Verify data equality with detailed diff on failure
        np.testing.assert_array_equal(
            out_data,
            ref_data,
            err_msg="Corrected rainbomb data differ between output and reference"
        )

        os.unlink(output_file.as_posix())
        # unlink all *.idx files created by cfgrib
        for idx_file in cur_folder.glob('*.idx'):
            os.unlink(idx_file.as_posix())
