import unittest
# from lisfloodutilities.compare.nc import NetCDFComparator
from lisfloodutilities.ncextract import read_points, read_inputmaps, read_ldd, extract_timeseries, rename_geographic_coords
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
        inputcsv = 'tests/data/ncextract/reservoirs.csv'
        data_dir = 'tests/data/ncextract/datasets'
        ldd_file = 'tests/data/ncextract/ldd.nc'
        expected_file = 'tests/data/ncextract/expected.nc'
        start = datetime(2018, 10, 2)
        end = datetime(2019, 10, 1)

        # read expected results
        expected = xr.open_dataset(expected_file)

        # read maps
        maps = read_inputmaps(data_dir, start=start, end=end)['dis24']

        # read points of interest
        poi = read_points(inputcsv)
        poi = rename_geographic_coords(poi, maps)

        # read LDD
        ldd = read_ldd(ldd_file)
        ldd = rename_geographic_coords(ldd, maps)

        # extract outflow timeseries
        print('Extracting reservoir outflow...')
        outflow = extract_timeseries(maps, poi, inflow=False)
        outflow.name = 'outflow'
        
        # extract inflow timeseries
        print('Extracting reservoir inflow...')
        inflow = extract_timeseries(maps, poi, inflow=True, ldd=ldd)
        inflow.name = 'inflow'
        
        # merge both extractions
        output = xr.merge((outflow, inflow))

        self.assertTrue(self.compare_datasets(output, expected))