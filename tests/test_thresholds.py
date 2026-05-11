import os
import numpy as np
import unittest
import xarray as xr

from lisfloodutilities.thresholds import lmoments, gumbel_parameters_moments, \
    gumbel_parameters_lmoments, gumbel_function
from lisfloodutilities.thresholds.thresholds import main


class TestLMoments(unittest.TestCase):
    def test_lmoments_shapes(self):
        # Test with input array of shape (100, 5)
        random = np.random.randn(100, 5)
        lambda_coef = lmoments(random)
        self.assertEqual(lambda_coef.shape, (2, 5), "Expected shape (2, 5) for input array of shape (100, 5)")

        # Test with input array of shape (200, 2)
        random = np.random.randn(200, 2)
        lambda_coef = lmoments(random)
        self.assertEqual(lambda_coef.shape, (2, 2), "Expected shape (2, 2) for input array of shape (200, 2)")

        # Test with input array of shape (3,)
        random = np.random.randn(3)
        with self.assertRaises(ValueError):
            lambda_coef = lmoments(random)

    def test_lmoments(self):
        values = np.arange(1, 100)
        lambda_coef = lmoments(values)
        print(lambda_coef)
        self.assertEqual(lambda_coef[0], 50., "Expected lambda_1 to be 50 for input values from 1 to 99")
        self.assertEqual(lambda_coef[1], 16.66666666666667, "Expected lambda_2 to be 16.66666666666667 for input values from 1 to 99")


class TestGumbelParameters(unittest.TestCase):
    def setUp(self):
        self.DECIMAL_CASES = 2 # Due to randomness of the input data, we round the output to 2 decimal cases for testing
        self.N_ELEMS = 100
        random = np.random.randn(self.N_ELEMS, 5)
        self.random = xr.DataArray(
            random,
            dims=['time', 'feature'],
            coords={
                'time': np.arange(self.N_ELEMS),
                'feature': ['a', 'b', 'c', 'd', 'e']
                }
            )
        values = np.arange(1, self.N_ELEMS)
        self.values = xr.DataArray(
            values,
            dims=['time'],
            coords={'time': np.arange(1, self.N_ELEMS)}
        )

    def test_moments_random(self):
        # Test output shapes
        parameters = gumbel_parameters_moments(self.random)
        sigma, mu = parameters['sigma'], parameters['mu']
        self.assertEqual(sigma.shape, (5,), f"Expected shape (5,) for sigma with input array of shape ({self.N_ELEMS}, 5)")
        self.assertEqual(mu.shape, (5,), f"Expected shape (5,) for mu with input array of shape ({self.N_ELEMS}, 5)")

    def test_lmoments_random(self):
        # Test output shapes
        parameters = gumbel_parameters_lmoments(self.random)
        sigma, mu = parameters['sigma'], parameters['mu']
        self.assertEqual(sigma.shape, (5,), f"Expected shape (5,) for sigma with input array of shape ({self.N_ELEMS}, 5)")
        self.assertEqual(mu.shape, (5,), f"Expected shape (5,) for mu with input array of shape ({self.N_ELEMS}, 5)")

    def test_moments(self):
        parameters = gumbel_parameters_moments(self.values)
        sigma = np.round(parameters['sigma'].values, self.DECIMAL_CASES)
        mu = np.round(parameters['mu'].values, self.DECIMAL_CASES)
        print('test_moments', sigma, mu)
        # Test output values
        test_sigma = np.round(22.395085599960808, self.DECIMAL_CASES)
        test_mu = np.round(37.07355659170262, self.DECIMAL_CASES)
        self.assertEqual(sigma, test_sigma, "Expected sigma value does not match")
        self.assertEqual(mu, test_mu, "Expected mu value does not match")

    def test_lmoments(self):
        parameters = gumbel_parameters_lmoments(self.values)
        sigma = np.round(parameters['sigma'].values, self.DECIMAL_CASES)
        mu = np.round(parameters['mu'].values, self.DECIMAL_CASES)
        print('test_lmoments', sigma, mu)
        # Test output values
        test_sigma = np.round(24.044917348149397, self.DECIMAL_CASES)
        test_mu = np.round(36.12127370664817, self.DECIMAL_CASES)
        self.assertEqual(sigma, test_sigma, "Expected sigma value does not match")
        self.assertEqual(mu, test_mu, "Expected mu value does not match")


class TestGumbelFunction(unittest.TestCase):
    def test_gumbel_function(self):
        y = np.array([2, 3, 4, 5])
        sigma = 1.0
        mu = 0.0
        result = gumbel_function(y, sigma, mu)
        if isinstance(result, np.ndarray):
            self.assertEqual(result.shape, (4,), "Expected shape (4,) for gumbel function output")
        print(result)
        truth = np.array([0.36651292, 0.90272046, 1.24589932, 1.49993999])
        np.testing.assert_allclose(result, truth, rtol=1e-6, err_msg="Expected gumbel function output does not match")

class TestFullTool(unittest.TestCase):
    def test_without_split(self):
        # Test the full tool without splitting the data
        # Get the absolute paths for input, output, and reference files
        base_path = os.path.dirname(__file__)
        input_file = os.path.join(base_path, 'data/thresholds/input/dis_yearmax.nc')
        output_file = os.path.join(base_path, 'data/thresholds/output/dis_THRESHOLDS.nc')
        reference_file = os.path.join(base_path, 'data/thresholds/reference/dis_THRESHOLDS.nc')

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Call the main function with the argv array of parameters
        argv = [
            'thresholds',
            '--discharge', input_file,
            '--output', output_file
        ]
        main(argv)

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

    def test_with_split(self):
        # Test the full tool with splitting the data
        # Get the absolute paths for input, output, and reference files
        base_path = os.path.dirname(__file__)
        input_file = os.path.join(base_path, 'data/thresholds/input/dis_yearmax.nc')
        output_file = os.path.join(base_path, 'data/thresholds/output/dis_THRESHOLDS.nc')
        reference_dir = os.path.join(base_path, 'data/thresholds/reference')

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Call the main function with the argv array of parameters including --split
        argv = [
            'thresholds',
            '--discharge', input_file,
            '--output', output_file,
            '--split'
        ]
        main(argv)

        # Get the base path without extension
        base_path = output_file.rsplit('.nc', 1)[0]

        # List of expected split files based on the reference directory
        expected_vars = ['mu', 'sigma', 'rp_1.5', 'rp_2', 'rp_5', 'rp_10', 'rp_20', 'rp_50', 'rp_100', 'rp_200', 'rp_500']

        # Compare each split output file with the corresponding reference file
        for var in expected_vars:
            split_output_file = f"{base_path}_{var}.nc"
            reference_file = os.path.join(reference_dir, f"dis_THRESHOLDS_{var}.nc")

            # Read output and reference netcdf files
            output_ds = xr.open_dataset(split_output_file)
            reference_ds = xr.open_dataset(reference_file)

            # Compare the data variable
            output_var = list(output_ds.data_vars)[0]
            self.assertIn(output_var, reference_ds.data_vars, f"Variable {output_var} not found in reference")

            np.testing.assert_allclose(
                output_ds[output_var].values,
                reference_ds[output_var].values,
                rtol=1e-5,
                err_msg=f"Variable {output_var} does not match reference"
            )

            # Close datasets
            output_ds.close()
            reference_ds.close()
