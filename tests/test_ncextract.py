import unittest
# from lisfloodutilities.compare.nc import NetCDFComparator
from lisfloodutilities.ncextract import read_points, read_inputmaps, extract_timeseries
import numpy as np
import xarray as xr
from datetime import datetime


class TestExtract(unittest.TestCase):

    def compare_datasets(
        self,
        dataset1: xr.Dataset,
        dataset2: xr.Dataset
    ) -> bool:
        """
        Compare the content of two xarray.Datasets and return True if they are identical, False otherwise.
        """
        # Check if both datasets have the same dimensions
        if set(dataset1.dims) != set(dataset2.dims):
            return False

        # Check if both datasets have the same coordinates
        if not all((dataset1.coords[dim] == dataset2.coords[dim]).all() for dim in dataset1.coords):
            return False

        # Check if both datasets have the same variables
        if set(dataset1.data_vars) != set(dataset2.data_vars):
            return False

        # Check if the variables' values are the same
        for var in dataset1.data_vars:
            if not np.allclose(dataset1[var].values, dataset2[var].values):
                return False
            
        # all tests passed
        return True 

    def test_ncextract(self):
    
        # config
        inputcsv = './data/ncextract/stations.csv'
        data_dir = './data/ncextract/datasets'
        expected_file = './data/ncextract/expected.nc'
        start = datetime(2018, 10, 1)
        end = datetime(2019, 9, 30)
        
        # read expected results
        expected = xr.open_dataset(expected_file)
        
        # read points of interest
        poi = read_points(inputcsv)
        
        # read maps
        maps = read_inputmaps(data_dir, start=start, end=end)
        
        # extract timeseries
        output = extract_timeseries(maps, poi)
        
        self.assertTrue(self.compare_datasets(output, expected))
