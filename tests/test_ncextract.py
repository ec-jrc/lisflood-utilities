import unittest
from lisfloodutilities.compare.nc import NetCDFComparator

from lisfloodutilities.ncextract import extract
import csv

class TestExtract(unittest.TestCase):
    def compare_csv_files(self, file1, file2):
        """
        Compare the content of two CSV files and return True if they are identical, False otherwise.
        """
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            reader1 = csv.reader(f1)
            reader2 = csv.reader(f2)
            
            for row1, row2 in zip(reader1, reader2):
                if row1 != row2:
                    return False
            
            # Check if both files have the same number of rows
            if len(list(reader1)) != len(list(reader2)):
                return False
            
        return True

    def test_extract_csv(self):
        inputcsv = 'tests/data/ncextract/stations.csv'
        datasets = 'tests/data/ncextract/datasets'
        outputfile = 'tests/data/output.csv'
        expected = 'tests/data/ncextract/expected.csv'
        extract(inputcsv, datasets, outputfile, nc=False)
        assert self.compare_csv_files(outputfile, expected)
        
    # def test_extract_nc(self):
    #     inputcsv = 'tests/data/ncextract/stations.csv'
    #     datasets = 'tests/data/ncextract/datasets'
    #     outputfile = 'tests/data/output.nc'
    #     expected = 'tests/data/ncextract/expected.nc'
    #     extract(inputcsv, datasets, outputfile, nc=True)
    #     comp = NetCDFComparator(None, for_testing=True)
    #     comp.compare_files(outputfile, expected)
    #     assert comp.errors == None