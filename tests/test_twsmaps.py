import os
import sys
import subprocess
import xarray as xr


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
