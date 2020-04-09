import pytest

from lisfloodutilities.compare import TSSComparator, NetCDFComparator, PCRComparator


class TestComparators:

    def test_netcdfcomp_files(self):
        comp = NetCDFComparator('tests/data/area.nc')
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc')
            assert 'have different variables names' in excinfo.value

    def test_netcdfcomp_dirs(self):
        comp = NetCDFComparator('tests/data/area.nc')
        comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=False)
            assert 'is missing in tests/data/folder_b/' in excinfo.value

    def test_netcdfcomp_tol(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.05, rtol=0.1)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.005, rtol=0.01)
        comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.0005, rtol=0.001)
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
            assert 'of different values - max diff:' in excinfo.value

    def test_netcdfcomp_identical(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', array_equal=True)
        comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
            assert 'Arrays are not equal' in excinfo.value

    def test_tss(self):
        comp = TSSComparator()
        comp.compare_files('tests/data/folder_a/qLakeIn.tss', 'tests/data/folder_b/qLakeIn.tss')
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_files('tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_b/qLakeOut.tss')
            assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in excinfo.value
        with pytest.raises(AssertionError) as excinfo:
            comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
            assert 'tests/data/folder_a/qLakeOut.tss different from tests/data/folder_b/qLakeOut.tss' in excinfo.value

    def test_pcr(self):
        comp = PCRComparator()
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
