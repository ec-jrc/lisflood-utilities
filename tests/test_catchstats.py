import unittest
from pathlib import Path
import pandas as pd
import pandas.testing as pdt
import sys
sys.path.append('../src/lisfloodutilities/catchstats')
from catchstats import _read_inputmaps, _read_masks, _read_pixarea, catchment_statistics

class TestCatchStats(unittest.TestCase):

    path = Path('data/catchstats')

    def test_catchstats(self):

        # compute test values
        maps = _read_inputmaps(path / 'maps')
        masks = _read_masks(path / 'masks')
        weight = _read_pixarea(path / 'pixarea_iberian_01min.nc')
        test = catchment_statistics(maps, masks, ['mean', 'std', 'min', 'max', 'count'], weight=weight, output=None).to_pandas()

        # load expected values
        expected = pd.read_csv(path / 'expected.csv', index_col='id')
        expected.index =expected.index.astype(test.index.dtype)

        # check
        try:
            pdt.assert_frame_equal(test, expected)
        except AssertionError as e:
            print(e)

