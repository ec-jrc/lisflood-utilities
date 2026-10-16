import os
import unittest

import numpy as np
import xarray as xr

from lisfloodutilities.adjust_trends.adjust_trends import main


class TestAdjustTrends(unittest.TestCase):
    def test_main(self):
        # Get the absolute paths for input, output, and reference files
        base_path = os.path.dirname(__file__)
        input_file = os.path.join(base_path, 'data/adjust_trends/input/dis_yearmax.nc')
        output_file = os.path.join(base_path, 'data/adjust_trends/output/dis_yearmax_adjusted.nc')
        reference_file = os.path.join(base_path, 'data/adjust_trends/reference/dis_yearmax_adjusted.nc')

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Call the main function with the argv array of parameters
        argv = [
            'adjust_trends',
            input_file,
            output_file,
            '--variable', 'var7',
            '--reference-year', '2010',
        ]
        ret = main(argv)
        self.assertEqual(ret, 0)

        # Read output and reference netcdf files
        output_ds = xr.open_dataset(output_file)
        reference_ds = xr.open_dataset(reference_file)

        # Compare corresponding variables of output and reference
        for var in output_ds.data_vars:
            self.assertIn(var, reference_ds.data_vars, f"Variable {var} not found in reference")
            np.testing.assert_allclose(
                output_ds[var].values,
                reference_ds[var].values,
                rtol=1e-5,
                err_msg=f"Variable {var} does not match reference"
            )

        # Close datasets
        output_ds.close()
        reference_ds.close()

        # Clean up the generated output
        os.remove(output_file)

    def test_main_diagnostics(self):
        # Get the absolute paths for input, output, and reference files
        base_path = os.path.dirname(__file__)
        input_file = os.path.join(base_path, 'data/adjust_trends/input/dis_yearmax.nc')
        output_file = os.path.join(base_path, 'data/adjust_trends/output/dis_yearmax_adjusted_diag.nc')
        reference_file = os.path.join(base_path, 'data/adjust_trends/reference/dis_yearmax_adjusted_diag.nc')

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Call the main function with the argv array of parameters, including --diagnostics
        argv = [
            'adjust_trends',
            input_file,
            output_file,
            '--variable', 'var7',
            '--reference-year', '2010',
            '--diagnostics',
        ]
        ret = main(argv)
        self.assertEqual(ret, 0)

        # Read output and reference netcdf files
        output_ds = xr.open_dataset(output_file)
        reference_ds = xr.open_dataset(reference_file)

        # Compare corresponding variables of output and reference
        for var in output_ds.data_vars:
            self.assertIn(var, reference_ds.data_vars, f"Variable {var} not found in reference")
            np.testing.assert_allclose(
                output_ds[var].values,
                reference_ds[var].values,
                rtol=1e-5,
                err_msg=f"Variable {var} does not match reference"
            )

        # Close datasets
        output_ds.close()
        reference_ds.close()

        # Clean up the generated output
        os.remove(output_file)


if __name__ == '__main__':
    unittest.main()
