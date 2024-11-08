import netCDF4 as nc
import numpy as np
from pathlib import Path

OUTPUT_PATH = '/mnt/nahaUsers/gomesgo/ERA5_vs_EFAS/var_statistics2'

data = {
    # 'EFAS': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/EMO-1arcmin-{v}_{yyyy}_monsum.nc', ['pr6', 'e0', 'et'], ['2019', '2020', '2022'], 0, 12),
    # 'EFAS_FRANCESCA_ATOS_2023': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_{yyyy}_from_francesca_and_ATOS_monsum.nc', ['pr', 'e0', 'et'], ['2023'], 0, 7),
    # 'EFAS_FRANCESCA_ATOS_2024': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_{yyyy}_from_francesca_and_ATOS_monsum.nc', ['pr', 'e0', 'et'], ['2024'], 0, 3),
    # 'EFAS_KISTERS_REFRESH_2023': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202307_202403_mm_per_day_monsum.nc', ['pr6'], ['2023'], 0, 6),
    # 'EFAS_KISTERS_REFRESH_2024': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202307_202403_mm_per_day_monsum.nc', ['pr6'], ['2024'], 6, 9),
    'EFAS_KISTERS_REFRESH_202404': ('/mnt/nahaUsers/grimast/meteo_template/meteo_from_francesca/{v}_EFASv5_Kisters_refresh_202404_mm_per_day_monsum.nc', ['pr6'], ['2024'], 0, 1),
    }

pixarea_filename = '/mnt/nahaUsers/gomesgo/ERA5_vs_EFAS/pixarea.nc'
flip_pixarea_updown = False
# pixarea_filename = '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/PixAREA_01degrees.tif.nc'
# flip_pixarea_updown = True

pixarea_var = 'Band1'


mask_var = 'Band1'



mask_filenames = [
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Elbe_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Vistula_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Europe_commonareas_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/Po_mask_1arcmin_LARGE.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Oder_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_EmiliaRomagna_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Piemonte_1arcmin.tif.nc',
    '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Lombardia_1arcmin.tif.nc',
    ]


# mask_filenames = [
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Elbe_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Vistula_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Europe_commonareas_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Oder_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/test_maskmap_PoValley_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_EmiliaRomagna_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Piemonte_01degrees.tif.nc',
#     '/mnt/nahaUsers/grimast/meteo_template/Goncalo_monthly_statistics/masks/mask_Lombardia_01degrees.tif.nc',
#     ]


for mask_filename in mask_filenames:
    f=Path(mask_filename)
    output_basename=f.stem

    # Open the mask NetCDF file
    mask_nc = nc.Dataset(mask_filename, 'r')
    # Read the mask data
    mask_data = mask_nc.variables[mask_var][:]
    mask_data = np.flipud(mask_data)

    # Open the Pixel Area NetCDF file
    pixarea_nc = nc.Dataset(pixarea_filename, 'r')
    # Read the pixel area data
    pixarea_data = pixarea_nc.variables[pixarea_var][:]
    if flip_pixarea_updown:
        pixarea_data = np.flipud(pixarea_data)
    masked_pixarea = np.ma.masked_where(mask_data !=1, pixarea_data)

    for dataset_name in data:
        data_filename_pattern, vars_array, years, start_step, end_step = data[dataset_name]
    
        for v in vars_array:
            output_file_path = Path(OUTPUT_PATH, f'{v}_{output_basename}.csv')
            with open(output_file_path, 'a') as out_file:
                for yyyy in years:
                    data_filename = data_filename_pattern.format(yyyy=yyyy, v=v)
            
                    # Open the data NetCDF file
                    data_nc = nc.Dataset(data_filename, 'r')
            
                    # Assuming the variable containing the monthly data is named 'monthly_data'
                    # and the mask variable is named 'mask'
                    monthly_data_var = v
            
                    # data_lat = data_nc.variables['lat'][:]
                    # data_lon = data_nc.variables['lon'][:]
                    # mask_lat = mask_nc.variables['lat'][:]
                    # mask_lon = mask_nc.variables['lon'][:]
            #        print('data_lat', data_lat)
            #        print('data_lon', data_lon)
            #        print('mask_lat', mask_lat)
            #        print('mask_lon', mask_lon)
                    # pixarea_lat = pixarea_nc.variables['lat'][:]
                    # pixarea_lon = pixarea_nc.variables['lon'][:]
            #        print('pixarea_lat', pixarea_lat)
            #        print('pixarea_lon', pixarea_lon)
            
            
                    # Initialize an empty list to store the yearly totals
                    yearly_totals = []
            
                    # Loop through each timeslice (assuming there are 12 for each year)
                    for i in range(start_step, end_step):
                        # Read the data for the current timeslice
                        monthly_data = data_nc.variables[monthly_data_var][i, :, :]
            
                        # Apply the mask to the data subset
                        # masked_monthly_data = np.where(mask_data, data_subset, 0)
                        masked_monthly_data = np.ma.masked_where(mask_data !=1, monthly_data)
            
                        monthly_sum_m = masked_monthly_data / 1000 / 4
            
                        monthly_sum_m3 = monthly_sum_m * masked_pixarea
            
                        monthly_totals_m3 = np.sum(monthly_sum_m3)
                        # monthly_totals_m3 = np.mean(masked_monthly_data)
                        
                        # print(f"Total values of {v} yyyy-mm {yyyy}-{i}: {monthly_totals_m3}")
                        # out_file.write(f"{dataset_name};{yyyy};{i+1};{monthly_totals_m3}\n")
                        
                        out_file.write(f"EFAS;{yyyy};{i+1+3};{monthly_totals_m3}\n")
                        # if yyyy == '2023':
                        #     out_file.write(f"EFAS;{yyyy};{i+1+6};{monthly_totals_m3}\n")
                        # elif yyyy == '2024':
                        #     out_file.write(f"EFAS;{yyyy};{i+1-6};{monthly_totals_m3}\n")
                                    
                        # If it's the first month, initialize the yearly sum with the masked data
                        if i == 0:
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
            
                    # print(f"Total values of {v} year {yyyy}: {total_values}")
    
    mask_nc.close()
    pixarea_nc.close()
