import os
import xarray as xr
from lisfloodutilities.mctrivers.mctrivers import mct_mask


class TestMctMask:
    
    pardir = 'tests/data/mctrivers'
    
    def compare_masks(self, checked_case):
        'check if the generated mct mask is identical to the expected ones'
        gnr_mask = f'{pardir}/{checked_case}/mctmask.nc'
        if os.path.exists(gnr_mask):
            os.remove(gnr_mask)
    
        # command = 'python3 src/lisfloodutilities/mctrivers/mctrivers.py'
        # command += f' -i {pardir}/{checked_case}/changrad.nc -l {pardir}/{checked_case}/ldd.nc -u {pardir}/{checked_case}/upArea.nc'
        # command += f' -m {pardir}/{checked_case}/mask.nc -S 0.001 -N 5 -U 500000000 -O {gnr_mask}'
        # os.system(command)
        mct_mask(f'{pardir}/{checked_case}/changrad.nc', f'{pardir}/{checked_case}/ldd.nc', f'{pardir}/{checked_case}/upArea.nc', 
                 mask_file=f'{pardir}/{checked_case}/mask.nc', slp_threshold=0.001, nloops=5, minuparea=500000000, 
                 outputfile=gnr_mask)
    
        org_mask = f'{pardir}/{checked_case}/mctmask_original.nc'
        original = xr.open_dataset(org_mask)
        generated = xr.open_dataset(gnr_mask)
    
        all_equal = original.equals(generated)
        if all_equal:
            print(f'Test for mct river mask generation for {checked_case} passed. Generated mask is deleted.')
            os.remove(gnr_mask)
        else:
            fail_message = f'Test for mct river mask generation for {checked_case} failed.' 
            fail_message += f'Please check differences between the generated mask "{gnr_mask}" and the expected mask "{org_mask}".'
            assert all_equal, fail_message


    def test_mctmask(self):
        'check mct masks for all available domains'
        available_cases = os.listdir(pardir)
        for i_case in available_cases:
            self.compare_masks(i_case)