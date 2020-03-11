from lisfloodutilities.compare import TSSComparator, NetCDFComparator, PCRComparator


class TestComparators:

    def test_netcdfcomp_files(self):
        comp = NetCDFComparator('tests/data/area.nc')
        assert not comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_b/ta.nc')
        assert comp.compare_files('tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc')

    def test_netcdfcomp_dirs(self):
        comp = NetCDFComparator('tests/data/area.nc')
        assert not comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
        assert comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/', skip_missing=False)

    def test_netcdfcomp_atol(self):
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.05, rtol=0.1)
        assert not comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.005, rtol=0.01)
        assert not comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')
        comp = NetCDFComparator('tests/data/folder_d/mask.nc', atol=0.0005, rtol=0.001)
        assert comp.compare_files('tests/data/folder_d/a.nc', 'tests/data/folder_d/b.nc')

    def test_tss(self):
        comp =TSSComparator()
        assert not comp.compare_files('tests/data/folder_a/qLakeIn.tss', 'tests/data/folder_b/qLakeIn.tss')
        assert comp.compare_files('tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_b/qLakeOut.tss')
        assert comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')

    def test_pcr(self):
        comp = PCRComparator()
        assert comp.compare_files('tests/data/folder_b/1.map', 'tests/data/folder_c/1.map')
        assert not comp.compare_files('tests/data/folder_a/1.map', 'tests/data/folder_c/1.map')
        assert comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_b/')
        assert not comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/')
        assert comp.compare_dirs('tests/data/folder_a/', 'tests/data/folder_c/', skip_missing=False)
