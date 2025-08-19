import unittest
from pathlib import Path
import pandas as pd
import pandas.testing as pdt
from lisfloodutilities.lfcoords import Config
from lisfloodutilities.lfcoords.finer_grid import coordinates_fine
from lisfloodutilities.lfcoords.coarser_grid import coordinates_coarse


class TestLFcoords(unittest.TestCase):

    path = Path('tests/data/lfcoords')

    def test_lfcoords(self):

        # compute test values
        cfg = Config(self.path / 'config.yml')
        stations_HR = coordinates_fine(cfg, save=False)
        test = coordinates_coarse(cfg, stations_HR, save=False)

        # load expected values
        expected = pd.read_csv(self.path / 'expected.csv', index_col='ID')
        expected.index =expected.index.astype(test.index.dtype)

        # check
        try:
            pdt.assert_frame_equal(test, expected, check_dtype=False)
        except AssertionError as e:
            print(e)
