import os
import sys
import subprocess
import numpy as np
import xarray as xr

from lisfloodutilities.twsmaps.twsmaps import delete_double


class TestTwsMaps:
    """
    Test the twsmaps lakes/reservoirs extent map generation by comparing the
    output of the command against a reference map stored in tests/data/twsmaps.
    See .kiro/steering/twsmaps-testing.md for how the reference maps are generated.

    The command is run as a subprocess (as it is meant to be used from the command
    line). This also isolates the test from a harmless segmentation fault that the
    GDAL/netCDF libraries can raise at interpreter shutdown.
    """

    data_dir = os.path.join(os.path.dirname(__file__), 'data', 'twsmaps')
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def _run_and_compare(self, case_dir):
        """Run twsmaps on the inputs in case_dir and compare with case_dir/reference.nc."""
        out_file = os.path.join(case_dir, 'out.nc')
        ref_file = os.path.join(case_dir, 'reference.nc')

        # generate the extent map with the command (run as a subprocess)
        cmd = [
            sys.executable, '-m', 'lisfloodutilities.twsmaps.twsmaps',
            '-d', 'GloFAS', '-t', 'reservoir', '-g', '-o', out_file,
            '--file-shp1', os.path.join(case_dir, 'grand.shp'),
            '--file-shp2', os.path.join(case_dir, 'hylak.shp'),
            '--file-shp3', os.path.join(case_dir, 'glwd.shp'),
            '--file-tab', os.path.join(case_dir, 'table.xlsx'),
            '--file-loc', os.path.join(case_dir, 'locations.nc'),
        ]
        subprocess.run(cmd, check=False, cwd=self.root_dir)

        assert os.path.exists(out_file), f'twsmaps did not produce the output map "{out_file}".'

        # compare the generated map with the reference one
        reference = xr.open_dataset(ref_file)
        generated = xr.open_dataset(out_file)
        all_equal = reference.equals(generated)
        reference.close()
        generated.close()

        # clean up the generated output
        try:
            os.remove(out_file)
        except FileNotFoundError:
            pass

        fail_message = (f'twsmaps extent map generation failed. Please check the differences '
                        f'between the generated map "{out_file}" and the reference map "{ref_file}".')
        assert all_equal, fail_message

    def test_twsmaps(self):
        """Standard case: reservoir with a polygon in the shapefile."""
        self._run_and_compare(self.data_dir)

    def test_twsmaps_polygon_not_found(self):
        """
        A reservoir that is present in the ID table but has no polygon in any shapefile
        must still be assigned to its outlet location pixel. This checks the fallback
        added to twsmaps ("POLYGON NOT FOUND, lake/res assigned to outlet loc pixel").
        """
        self._run_and_compare(os.path.join(self.data_dir, 'polygon_not_found'))


class TestDeleteDouble:
    """
    Unit tests for the delete_double workaround that removes duplicated source IDs
    (e.g. GRAND_ID) from the locations map, keeping only the LISFLOOD ID with the
    largest area in the ID match table.
    """

    @staticmethod
    def _loc_map():
        loc = np.zeros((3, 3), dtype=np.int_)
        loc[1, 1] = 400
        loc[2, 2] = 401
        return loc

    def test_no_double_ids(self):
        """A table without duplicated source IDs must not modify the locations map."""
        loc = self._loc_map()
        id_src = np.array(['GRAND', 'HYLAK'])
        result = delete_double(np.array([400, 401]), np.array([1, 2]), loc, 'GRAND', id_src, np.array([10., 20.]))
        assert np.array_equal(result, loc)

    def test_double_id_keeps_largest_area(self):
        """Two LISFLOOD IDs sharing a source ID: the largest-area one is kept."""
        loc = self._loc_map()
        id_src = np.array(['GRAND', 'GRAND'])
        result = delete_double(np.array([400, 401]), np.array([7, 7]), loc, 'GRAND', id_src, np.array([10., 50.]))
        assert result[1, 1] == 0
        assert result[2, 2] == 401

    def test_double_id_other_source(self):
        """Duplicated IDs whose rows do not use the given source must be left untouched."""
        loc = self._loc_map()
        id_src = np.array(['HYLAK', 'HYLAK'])
        result = delete_double(np.array([400, 401]), np.array([7, 7]), loc, 'GRAND', id_src, np.array([10., 50.]))
        assert np.array_equal(result, loc)

    def test_double_id_after_nan_rows(self):
        """
        Regression test: the source check must use the position of the ID in the
        original table, not its position among the non-NaN entries (NaNs before the
        duplicated ID would otherwise make the check look at the wrong row).
        """
        loc = np.zeros((3, 3), dtype=np.int_)
        loc[0, 0] = 400
        loc[1, 1] = 401
        loc[2, 2] = 402
        id_oth = np.array([np.nan, 7, 7])
        id_src = np.array(['HYLAK', 'GRAND', 'GRAND'])
        result = delete_double(np.array([400, 401, 402]), id_oth, loc, 'GRAND', id_src,
                               np.array([np.nan, 10., 50.]))
        # the duplicated GRAND ID is used from GRAND rows: 401 (area 10) is removed
        # because 402 (area 50) is kept, while the HYLAK row (400) is untouched
        assert result[0, 0] == 400
        assert result[1, 1] == 0
        assert result[2, 2] == 402


class TestTwsMapsEfas:
    """
    End-to-end test of the EFAS domain: it also covers the delete_double workaround
    (two LISFLOOD IDs sharing a GRAND_ID, only the largest-area one is kept) and the
    GDW source branch.
    """

    data_dir = os.path.join(os.path.dirname(__file__), 'data', 'twsmaps', 'efas_double_ids')
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def test_twsmaps_efas_double_ids(self):
        out_file = os.path.join(self.data_dir, 'out.nc')
        ref_file = os.path.join(self.data_dir, 'reference.nc')

        cmd = [
            sys.executable, '-m', 'lisfloodutilities.twsmaps.twsmaps',
            '-d', 'EFAS', '-t', 'reservoir', '-g', '-o', out_file,
            '--file-shp1', os.path.join(self.data_dir, 'grand.shp'),
            '--file-shp2', os.path.join(self.data_dir, 'hylak.shp'),
            '--file-shp3', os.path.join(self.data_dir, 'gdw.shp'),
            '--file-tab', os.path.join(self.data_dir, 'table.xlsx'),
            '--file-loc', os.path.join(self.data_dir, 'locations.nc'),
        ]
        subprocess.run(cmd, check=False, cwd=self.root_dir)

        assert os.path.exists(out_file), f'twsmaps did not produce the output map "{out_file}".'

        generated = xr.open_dataset(out_file)
        polygon_id = generated['polygon_id'].values

        # the reservoir with the smallest area (300) is removed by delete_double
        assert not np.any(polygon_id == 300)
        # the kept reservoir (301) fills the cells covered by its GRanD polygon
        for row, col in [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]:
            assert polygon_id[row, col] == 301, f'unexpected id {polygon_id[row, col]} at (row {row}, col {col})'
        # the GDW-sourced reservoir (302) is assigned to its outlet location pixel
        assert polygon_id[0, 0] == 302

        if os.path.exists(ref_file):
            # compare the generated map with the reference one
            reference = xr.open_dataset(ref_file)
            all_equal = reference.equals(generated)
            reference.close()
            assert all_equal, (f'twsmaps EFAS extent map generation failed. Please check the differences '
                               f'between the generated map "{out_file}" and the reference map "{ref_file}".')

        generated.close()

        # clean up the generated output
        try:
            os.remove(out_file)
        except FileNotFoundError:
            pass
