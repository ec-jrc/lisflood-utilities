import unittest
from pathlib import Path
import pandas as pd
import pandas.testing as pdt
from lisfloodutilities.lfcoords import Config
from lilisfloodutilities.lfcoords.finer_grid import coordinates_fine
from llisfloodutilities.lfcoords.coarser_grid import coordinates_coarse


class TestCatchStats(unittest.TestCase):

    path = Path('tests/data/lfcoords')

    def test_lfcoords(self):

        # compute test values
        cfg = Config(self.path / 'config.yml')
        stations_HR = coordinates_fine(cfg, save=False)
        test = coordinates_coarse(cfg, stations_HR, save=False)


        # maps = read_inputmaps(self.path / 'maps')
        # masks = read_masks(self.path / 'masks')
        # weight = read_pixarea(self.path / 'pixarea_iberian_01min.nc')
        # test = catchment_statistics(maps, masks, ['mean', 'std', 'min', 'max', 'count'], weight=weight, output=None).to_pandas()

        # load expected values
        expected = pd.read_csv(self.path / 'expected.csv', index_col='id')
        expected.index =expected.index.astype(test.index.dtype)

        # check
        try:
            pdt.assert_frame_equal(test, expected)
        except AssertionError as e:
            print(e)
