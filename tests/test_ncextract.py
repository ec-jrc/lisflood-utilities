import unittest

from lisfloodutilities.ncextract import extract

class TestExtract(unittest.TestCase):
    def test_extract_csv(self):
        inputcsv = 'tests/data/ncextract/stations.csv'
        datasets = 'tests/data/ncextract'
        outputfile = 'tests/data/output.csv'
        expected = 'tests/data/expected.csv'
        extract(inputcsv, datasets, outputfile)
        assert outputfile == expected
        
    def test_extract_nc(self):
        inputcsv = 'tests/data/ncextract/stations.csv'
        datasets = 'tests/data/ncextract'
        outputfile = 'tests/data/output.nc'
        expected = 'tests/data/expected.csv'
        extract(inputcsv, datasets, outputfile, nc=True)
        assert outputfile == expected