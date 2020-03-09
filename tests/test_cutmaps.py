from pathlib import Path

from lisfloodutilities.cutmaps.cutlib import get_filelist, get_cuts, cutmap


class TestCutlib:

    def test_getfiles_to_cut_filelist(self):
        l = get_filelist(filelist='tests/data/cutlist.txt')
        assert list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc',
                               'tests/data/folder_b/ta.nc', 'tests/data/folder_d/a.nc',
                               'tests/data/folder_d/b.nc'])) == l

    def test_getfiles_to_cut_folder(self):
        l = get_filelist(input_folder='tests/data/folder_a')
        assert list(map(Path, ['tests/data/folder_a/ta.nc', 'tests/data/folder_a/tp.nc'])) == l


    def test_getfiles_to_cut_glofas_setup(self):
        l = get_filelist(glofas_folder='tests/data/folder_a')
        assert list(
            map(Path, ['tests/data/folder_a/1.map', 'tests/data/folder_a/2.map',
                       'tests/data/folder_a/3.map', 'tests/data/folder_a/4.map',
                       'tests/data/folder_a/5.map', 'tests/data/folder_a/qLakeIn.tss',
                       'tests/data/folder_a/qLakeOut.tss', 'tests/data/folder_a/ta.nc',
                       'tests/data/folder_a/tp.nc'])) == l

    def test_get_cuts_withcuts(self):
        pass

    def test_get_cuts_withmaskfile(self):
        pass
    #
    # def test_cutmap(self):
    #     pass
