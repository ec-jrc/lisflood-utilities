import numpy as np
import unittest

from lisfloodutilities.thresholds import lmoments, gumbel_parameters_moments, \
    gumbel_parameters_lmoments, gumbel_function


class TestLMoments(unittest.TestCase):
    def test_lmoments_shapes(self):
        # Test with input array of shape (10, 5)
        random = np.random.randn(100, 5)
        lambda_coef = lmoments(random)
        self.assertEqual(lambda_coef.shape, (2, 5))

        # Test with input array of shape (2, 2)
        random = np.random.randn(200, 2)
        lambda_coef = lmoments(random)
        self.assertEqual(lambda_coef.shape, (2, 2))

        # Test with input array of shape (3,)
        random = np.random.randn(3)
        with self.assertRaises(ValueError):
            lambda_coef = lmoments(random)

    def test_lmoments(self):
        values = np.arange(1, 100)
        lambda_coef = lmoments(values)
        print(lambda_coef)
        self.assertEqual(lambda_coef[0], 50.)
        self.assertEqual(lambda_coef[1], 16.66666666666667)


class TestGumbelParameters(unittest.TestCase):
    def setUp(self):
        self.random = np.random.randn(100, 5)
        self.values = np.arange(1, 100)

    def test_moments_random(self):
        # Test output shapes
        sigma, mu = gumbel_parameters_moments(self.random)
        self.assertEqual(sigma.shape, (5,))
        self.assertEqual(mu.shape, (5,))

    def test_lmoments_random(self):
        # Test output shapes
        sigma, mu = gumbel_parameters_lmoments(self.random)
        self.assertEqual(sigma.shape, (5,))
        self.assertEqual(mu.shape, (5,))

    def test_moments(self):
        sigma, mu = gumbel_parameters_moments(self.values)
        # Test output values
        self.assertEqual(sigma, 22.395085599960808)
        self.assertEqual(mu, 37.07355659170262)

    def test_lmoments(self):
        sigma, mu = gumbel_parameters_lmoments(self.values)
        # Test output values
        self.assertEqual(sigma, 24.044917348149397)
        self.assertEqual(mu, 36.12127370664817)


class TestGumbelFunction(unittest.TestCase):
    def test_gumbel_function(self):
        y = np.array([2, 3, 4, 5])
        sigma = 1.0
        mu = 0.0
        result = gumbel_function(y, sigma, mu)
        self.assertEqual(result.shape, (4,))
        print(result)
        truth = np.array([0.36651292, 0.90272046, 1.24589932, 1.49993999])
        np.testing.assert_allclose(result, truth)
