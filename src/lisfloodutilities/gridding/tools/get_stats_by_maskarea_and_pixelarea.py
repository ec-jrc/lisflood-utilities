import netCDF4 as nc
import numpy as np
from pathlib import Path

OUTPUT_PATH = '/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/var_statistics'

# For pr6, e0, es, et in mm/day
mm_per_day_to_m = 1000.0 / 4.0
# For pr6, e0, es, et in mm/6h
mm_to_m = 1000.0

# data_filename_pattern, vars_array, years, start_step, end_step, flip_updown, conversion_to_meters
data = {
    # 'KISTERS': ('/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/compare/{v}/KISTERS-{v}_{yyyy}*_monsum.nc', ['pr6'], ['2000'], 0, 0, False, conversion2meters),
    'KISTERS': ('/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/Kisters/decumulation/compare/{v}/KISTERS-{v}_{yyyy}*_monsum.nc', ['pr6'], ['2000', '2005', '2021'], 0, 0, False, mm_to_m),
    'EMO1': ('/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/compare/{v}/EMO-1arcmin-{v}_{yyyy}*_monsum.nc', ['pr6'], ['2000', '2005', '2021'], 0, 11, False, mm_per_day_to_m),
    # 'EFAS': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/EMO-1arcmin-{v}_{yyyy}_monsum.nc', ['pr6', 'e0', 'et'], ['2019', '2020', '2022'], 0, 11, False, mm_per_day_to_m),
    # 'EFAS_FRANCESCA_ATOS_2023': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_{yyyy}_from_francesca_and_ATOS_monsum.nc', ['pr', 'e0', 'et'], ['2023'], 0, 6, False, mm_per_day_to_m),
    # 'EFAS_FRANCESCA_ATOS_2024': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_{yyyy}_from_francesca_and_ATOS_monsum.nc', ['pr', 'e0', 'et'], ['2024'], 0, 2, False, mm_per_day_to_m),
    # 'EFAS_KISTERS_REFRESH_2023': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202307_202403_mm_per_day_monsum.nc', ['pr6'], ['2023'], 0, 5, False, mm_per_day_to_m),
    # 'EFAS_KISTERS_REFRESH_2024': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202307_202403_mm_per_day_monsum.nc', ['pr6'], ['2024'], 6, 8, False, mm_per_day_to_m),
    # 'EFAS_KISTERS_REFRESH_202404': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202404_mm_per_day_monsum.nc', ['pr6'], ['2024'], 0, 0, False, mm_per_day_to_m),
}

pixarea_filename = '/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/configuration/1arcmin/pixarea.nc'
flip_pixarea_updown = False
# pixarea_filename = '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/PixAREA_01degrees.tif.nc'
# flip_pixarea_updown = True
pixarea_var = 'Band1'

# masks_base_path = '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks'
masks_base_path = '/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/masks'

# mask_filename, mask_var, flip_mask_updown
mask_filenames = [
    (f'{masks_base_path}/mask_Elbe_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_Vistula_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_Europe_commonareas_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/Po_mask_1arcmin_LARGE.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_Oder_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_EmiliaRomagna_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_Piemonte_1arcmin.tif.nc', 'Band1', True),
    (f'{masks_base_path}/mask_Lombardia_1arcmin.tif.nc', 'Band1', True)
]


# mask_filenames = [
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Elbe_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Vistula_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Europe_commonareas_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Oder_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/test_maskmap_PoValley_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_EmiliaRomagna_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Piemonte_01degrees.tif.nc', 'Band1', True),
#     ('/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Lombardia_01degrees.tif.nc', 'Band1', True)
#     ]


for mask_filename, mask_var, flip_mask_updown in mask_filenames:
    f = Path(mask_filename)
    output_basename = f.stem

    # Open the mask NetCDF file
    mask_nc = nc.Dataset(mask_filename, 'r')
    # Read the mask data
    mask_data = mask_nc.variables[mask_var][:]
    if flip_mask_updown:
        mask_data = np.flipud(mask_data)

    # Open the Pixel Area NetCDF file
    pixarea_nc = nc.Dataset(pixarea_filename, 'r')
    # Read the pixel area data
    pixarea_data = pixarea_nc.variables[pixarea_var][:]
    if flip_pixarea_updown:
        pixarea_data = np.flipud(pixarea_data)
    masked_pixarea = np.ma.masked_where(mask_data !=1, pixarea_data)

    for dataset_name in data:
        print(f'Processing: {dataset_name} for {mask_filename}')
        data_filename_pattern, vars_array, years, start_step, end_step, flip_updown, conversion_to_meters = data[dataset_name]

        for v in vars_array:
            output_file_path = Path(OUTPUT_PATH, f'{v}_{output_basename}.tsv')
            with open(output_file_path, 'a') as out_file:
                for yyyy in years:
                    print(f'Processing: {v} {yyyy}')
                    current_filename_pattern = data_filename_pattern.format(yyyy=yyyy, v=v)
                    current_filename_pattern_path = Path(current_filename_pattern)
                    current_folder = current_filename_pattern_path.parent
                    current_filename_wildcard = current_filename_pattern_path.name

                    cur_idx = start_step
                    for data_filename in sorted(current_folder.rglob(current_filename_wildcard)):
                        print(f'Processing: {data_filename}')
                        # Open the data NetCDF file
                        data_nc = nc.Dataset(data_filename, 'r')

                        # Assuming the variable containing the monthly data is named 'monthly_data'
                        # and the mask variable is named 'mask'
                        monthly_data_var = v

                        # Initialize an empty list to store the yearly totals
                        yearly_totals = []

                        # Loop through each timeslice (assuming there are 12 for each year)
                        for i in range(start_step, end_step+1):
                            cur_idx += 1
                            # Read the data for the current timeslice
                            monthly_data = data_nc.variables[monthly_data_var][i, :, :]

                            if flip_updown:
                                monthly_data = np.flipud(monthly_data)

                            # Apply the mask to the data subset
                            # masked_monthly_data = np.where(mask_data, data_subset, 0)
                            masked_monthly_data = np.ma.masked_where(mask_data!=1, monthly_data)
                
                            monthly_sum_m = masked_monthly_data / conversion_to_meters
                
                            monthly_sum_m3 = monthly_sum_m * masked_pixarea
                
                            monthly_totals_m3 = np.sum(monthly_sum_m3)
                            # monthly_totals_m3 = np.mean(masked_monthly_data)
                            
                            out_file.write(f'{dataset_name}\t{yyyy}\t{cur_idx}\t{monthly_totals_m3}\n')

                            # If it's the first month, initialize the yearly sum with the masked data
                            if i == start_step:
                                yearly_sum = monthly_totals_m3
                            else:
                                # Otherwise, add the masked data to the running yearly sum
                                yearly_sum += monthly_totals_m3
                
                        # Add the yearly sum to the list of yearly totals
                        yearly_totals.append(yearly_sum)
                
                        # Convert the list to a numpy array if necessary
                        yearly_totals = np.array(yearly_totals)
                
                        total_values = np.sum(yearly_totals)
                
                        # Close the NetCDF files
                        data_nc.close()
                        data_nc = None
            
                    print(f"Total values of {v} year {yyyy}: {total_values}")
    
    mask_nc.close()
    pixarea_nc.close()
