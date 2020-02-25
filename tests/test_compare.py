from lisfloodutilities.compare import TSSComparator, NetCDFComparator, PCRComparator

class TestComparators:

    def test_netcdfcomp(self):
        comp = NetCDFComparator('tests/data/areaOrigin.nc')
        assert not comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        assert not comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
        assert comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=False)
        assert comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc')

