import os
from pathlib import Path
from lisfloodutilities.gridding.decumulate_daily_grids import main


cur_folder = os.path.dirname(os.path.realpath(__file__))


class TestDecumulate:

    def check_file_content(self, file1: Path, file2: Path):
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            content1 = f1.read()
            content2 = f2.read()
        file_name1 = str(file1).replace(self.base_data_folder, '')
        file_name2 = str(file2).replace(self.base_data_folder, '')
        assert content1 == content2, f'Files [{file_name1}] and [{file_name2}] do not have the same content.'
    
    def test_decumulate_kiwis(self):
        configuration_base_folder = Path(os.path.join(cur_folder, '../src/lisfloodutilities/gridding/configuration'))
        self.base_data_folder = os.path.join(cur_folder, 'data/decumulate/')
        input_folder = Path(os.path.join(self.base_data_folder, 'meteo_in'))
        reference_folder = Path(os.path.join(self.base_data_folder, 'reference/pr6'))
        output_folder = Path(os.path.join(self.base_data_folder, 'meteo_out/pr6'))

        output_folder.mkdir(parents=True, exist_ok=True)
        # Clean the output folder
        for filename_kiwis in output_folder.rglob('*.kiwis'):
            os.remove(filename_kiwis)
        
        argv = ['--conf', '1arcmin', '--pathconf', f'{configuration_base_folder}',
                '--var24h', 'pr', '--var6h', 'pr6', '--out', f'{output_folder}',
                '--pr24h', f'{input_folder}/pr/', '-s', '20051226000000', '-e', '20051226060001',
                '--pr6h', f'{input_folder}/pr6/']

        retval = main(argv)

        assert retval == 0, f'Got return value {retval}. Some error occurred while processing.'

        for filename_kiwis in reference_folder.rglob('*.kiwis'):
            out_file = str(filename_kiwis).replace(str(reference_folder), str(output_folder))
            assert os.path.exists(out_file), f'File [{out_file}] was not generated.'
            self.check_file_content(filename_kiwis, Path(out_file))
        
