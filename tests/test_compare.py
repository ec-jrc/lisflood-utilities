"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""
from nine import IS_PYTHON2
if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

import pytest
import shutil

from lisfloodutilities.compare.pcr import PCRComparator, TSSComparator
from lisfloodutilities.compare.nc import NetCDFComparator


class TestComparators:

    def test_netcdfcomp_files(self):
        comp = NetCDFComparator('tests/data/masks/area.nc', for_testing=True)
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc')
        assert 'have different variables names' in str(excinfo.value)

    def test_netcdfcomp_dirs(self):
        comp = NetCDFComparator('tests/data/masks/area.nc', for_testing=True)
        comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=True)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=False)
        assert 'is missing in tests/data/folder_b' in str(excinfo.value)

    def test_netcdfcomp_diffdirs(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', for_testing=False, save_diff_files=True, array_equal=True)
        assert comp.diff_folder is not None
        assert isinstance(comp.diff_folder, Path)
        assert comp.diff_folder.exists()
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/a.nc')
        assert not comp.diff_folder.joinpath('a.nc').exists()
        assert not comp.diff_folder.joinpath('b.nc').exists()
        assert not comp.diff_folder.joinpath('diff.nc').exists()
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        assert comp.diff_folder.joinpath('a_a.nc').exists(), 'Diff file A does not exist in %s' % comp.diff_folder
        assert comp.diff_folder.joinpath('a_b.nc').exists(), 'Diff file B does not exist'
        assert comp.diff_folder.joinpath('a_diff.nc').exists(), 'Diff file diff.nc does not exist'
        # Clean the diff folder after the test since it is no longer needed
        if comp.diff_folder is not None and comp.diff_folder.parent.name == 'diffs':
            shutil.rmtree(comp.diff_folder.parent)
        else:
            shutil.rmtree(comp.diff_folder)

    def test_netcdfcomp_submask(self):
        comp = NetCDFComparator(mask='tests/data/submask/subcatchment_mask.map')
        comp.compare_files('tests/data/submask/dis_subdomain.nc', 'tests/data/submask/dis.nc')

    def test_netcdfcomp_tol(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.05, rtol=0.1, for_testing=True)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.005, rtol=0.01, for_testing=True)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.0005, rtol=0.001, for_testing=True)
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        assert 'of different values - max diff:' in str(excinfo.value)

    def test_netcdfcomp_identical(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', array_equal=True, for_testing=True)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/a.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        assert 'Arrays are not equal' in str(excinfo.value)

    def test_netcdfcomp_no_asserts(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', array_equal=True, for_testing=False)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/a.nc')
        assert not comp.errors
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        assert comp.errors
        assert 'a.nc/field@0 is different' in comp.errors

    def test_netcdfcomp_no_mask(self):
        comp = NetCDFComparator(array_equal=True)
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/a.nc')
        comp = NetCDFComparator()
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/a.nc')
        with pytest.raises(ValueError) as excinfo:
            comp.compare_files('tests/data/submask/dis_subdomain.nc', 'tests/data/submask/dis.nc')
        assert 'operands could not be broadcast together with shapes (43,37) (57,80)' in str(excinfo.value)

    def test_netcdfcomp_timestep(self):
        # TODO
        pass

    def test_netcdfcomp_timestep_range(self):
        # TODO
        pass

    def test_tss(self):
        with pytest.raises(AssertionError) as excinfo:
            comp = TSSComparator(array_equal=True)
            comp.compare_files('tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_b/qLakeOut.tss')

        assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in str(excinfo.value)
        comp.compare_files('tests/data/folder_a/qLakeIn.tss', 'tests/data/folder_b/qLakeIn.tss')

        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=True)
        assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in str(excinfo.value)

    def test_tss_timestep(self):
        comp = TSSComparator(for_testing=True, array_equal=True)
        timestep = 62
        comp.compare_files('tests/data/folder_a/qLakeIn.tss', 'tests/data/folder_b/qLakeInSingleTimestep.tss', timestep)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_b/qLakeOut.tss', 'tests/data/folder_b/qLakeOut2.tss', timestep)
        assert '62 not found in tests/data/folder_b/qLakeOut2.tss' in str(excinfo.value)

    def test_tss_timestep_tol(self):
        comp = TSSComparator(for_testing=True, rtol=0.00001, atol=0.000001)
        timestep = 62
        comp.compare_files('tests/data/folder_b/qLakeIn.tss', 'tests/data/folder_b/qLakeInSingleTimestep.tss', timestep)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_b/qLakeOut.tss', 'tests/data/folder_b/qLakeOut3.tss', timestep)
        assert 'Not equal to tolerance rtol=1e-05, atol=1e-06\n\nMismatched elements: 1 / 2 (50%)' in str(excinfo.value)

    def test_tss_tolerance(self):
        comp = TSSComparator(for_testing=True, rtol=0.001, atol=0.0001)
        comp.compare_files('tests/data/folder_a/test_tol_ok_1.tss', 'tests/data/folder_b/test_tol_ok_2.tss')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/test_tol_fail_1.tss', 'tests/data/folder_b/test_tol_fail_2.tss')
        assert 'Not equal to tolerance rtol=0.001, atol=0.0001\n\nMismatched elements: 1 / 2 (50%)\nMax absolute difference: 0.007' in str(excinfo.value)
        comp = TSSComparator(for_testing=True, rtol=0.00001, atol=0.000001)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/test_tol_ok_1.tss', 'tests/data/folder_b/test_tol_ok_2.tss')
        assert 'Not equal to tolerance rtol=1e-05, atol=1e-06\n\nMismatched elements: 1 / 2 (50%)\nMax absolute difference: 5.e-05' in str(excinfo.value)

    def test_pcr(self):
        comp = PCRComparator(for_testing=True)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_b/1.map', 'tests/data/folder_c/1.map')
        assert 'tests/data/folder_b/1.map different from tests/data/folder_c/1.map' in str(excinfo.value)
        comp.compare_files('tests/data/folder_a/1.map', 'tests/data/folder_c/1.map')
        with pytest.raises(AssertionError, match=r"tests/data/folder_a/[1-4].map different from tests/data/folder_b/[1-4].map"):
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=True)
        comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/', skip_missing=True)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/', skip_missing=False)
        assert '.map is missing in tests/data/folder_c' in str(excinfo.value)
