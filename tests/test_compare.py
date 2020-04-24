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

import pytest

from lisfloodutilities.compare import TSSComparator, NetCDFComparator, PCRComparator


class TestComparators:

    def test_netcdfcomp_files(self):
        comp = NetCDFComparator('tests/data/area.nc', for_testing=True)
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc')
            assert 'have different variables names' in excinfo.value

    def test_netcdfcomp_dirs(self):
        comp = NetCDFComparator('tests/data/area.nc', for_testing=True)
        comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=False)
            assert 'is missing in tests/data/folder_b/' in excinfo.value

    def test_netcdfcomp_tol(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.05, rtol=0.1, for_testing=True)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.005, rtol=0.01, for_testing=True)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.0005, rtol=0.001, for_testing=True)
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
            assert 'of different values - max diff:' in excinfo.value

    def test_netcdfcomp_identical(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', array_equal=True, for_testing=True)
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
            assert 'Arrays are not equal' in excinfo.value

    def test_tss(self):
        comp = TSSComparator(for_testing=True, array_equal=True)
        comp.compare_files('tests/data/folder_a/qLakeIn.tss', 'tests/data/folder_b/qLakeIn.tss')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_b/qLakeOut.tss')
            assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in excinfo.value
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
            assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in excinfo.value

    def test_tss_tolerance(self):
        comp = TSSComparator(for_testing=True, rtol=0.001, atol=0.0001)
        comp.compare_files('tests/data/folder_a/test_tol_ok_1.tss', 'tests/data/folder_b/test_tol_ok_2.tss')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/test_tol_fail_1.tss', 'tests/data/folder_b/test_tol_fail_2.tss')
            assert 'tests/data/folder_a/test_tol_fail_1.tss different from tests/data/folder_b/test_tol_fail_2.tss (line: 4)' in excinfo.value
        comp = TSSComparator(for_testing=True, rtol=0.00001, atol=0.000001)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/test_tol_ok_1.tss', 'tests/data/folder_b/test_tol_ok_2.tss')
            assert 'tests/data/folder_a/test_tol_ok_1.tss different from tests/data/folder_b/test_tol_ok_2.tss (line: 4)' in excinfo.value

    def test_pcr(self):
        comp = PCRComparator(for_testing=True)
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_b/1.map', 'tests/data/folder_c/1.map')
            assert 'tests/data/folder_b/1.map different from tests/data/folder_c/1.map' in excinfo.value
        comp.compare_files('tests/data/folder_a/1.map', 'tests/data/folder_c/1.map')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
            assert 'tests/data/folder_b/1.map different from tests/data/folder_c/1.map' in excinfo.value
        comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/', skip_missing=False)
            assert '.map is missing in tests/data/folder_c' in excinfo.value

    def test_comp_no_asserts(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', array_equal=True, for_testing=False)
        errors = comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        assert not errors
        errors = comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        assert errors
        assert 'a.nc/field@0 is different' in errors
