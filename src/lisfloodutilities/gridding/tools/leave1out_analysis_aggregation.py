import pandas as pd


BASE_PATH = '/BGFS/DISASTER/grids/leave_one_out'

variable_name = 'pr6'

is_test_running=False

# interpolation_algorithms = ['spheremap', 'idw', 'adw', 'cdd']
interpolation_algorithms = ['cdd']

if not is_test_running:
    algorithm_folder_name = {
         'spheremap': 'sph_all',
         'idw': 'idw_all',
         'adw': 'adw_all',
         'cdd': 'cdd_all'
         }
else:
    algorithm_folder_name = {
        'spheremap': 'spheremap',
        'idw': 'idw',
        'adw': 'adw',
        'cdd': 'cdd'
    }

timesteps_file = f'{BASE_PATH}/all-time-steps.txt'
series = pd.read_csv(timesteps_file, delimiter='\t', header=None, squeeze=True)
timesteps = list(series.T)
series = None

df_dic = {
    'mae': None,
    'mbe': None,
    'mse': None,
    'pearson_r': None,
    'csi': None,
    'values_sum': None,
    'pixels_with_values_0<X<=0.1': None,
    'pixels_with_values_0.1<X<1': None
    }

df_file_name = {
    'mae': f'{variable_name}_mae.stats',
    'mbe': f'{variable_name}_mbe.stats',
    'mse': f'{variable_name}_mse.stats',
    'pearson_r': f'{variable_name}_pearson_r.stats',
    'csi': f'{variable_name}_csi.stats',
    'values_sum': f'{variable_name}_values_sum.stats',
    'pixels_with_values_0<X<=0.1': f'{variable_name}_pixels_with_values_between_0_and_0.1.stats',
    'pixels_with_values_0.1<X<1': f'{variable_name}_pixels_with_values_between_0.1_and_1.stats'
    }

print('Start reading...')
for timestep in timesteps:
    for interpolation_algorithm in interpolation_algorithms:
        # /BGFS/DISASTER/grids/leave_one_out/spheremap/pr6/pr6202110041800_spheremap_10mm.stats
        file_in = f'{BASE_PATH}/{algorithm_folder_name[interpolation_algorithm]}/{variable_name}/{variable_name}{timestep}_{interpolation_algorithm}.stats'
        print(f'reading: {file_in}')
        df = pd.read_csv(file_in, delimiter='\t')
        for df_column in df_dic:
            df_out = df_dic[df_column]
            if df_out is None:
                column_names = ['timestep.interpolation']
                column_names.extend(df['test'].tolist())
                df_out = pd.DataFrame(columns=column_names)
            new_row = [f'{timestep}.{interpolation_algorithm}']
            new_row.extend(df[df_column].tolist())
            df_out.loc[len(df_out)] = new_row
            df_dic[df_column] = df_out
        df = None

print('Writing results...')
# Write results
for df_column in df_dic:
    df_out = df_dic[df_column]
    file_out = df_file_name[df_column]
    df_out = df_out.set_index('timestep.interpolation')
    df_out = df_out.replace(-9999.0, float('nan'))
    df_out['mean'] = df_out.loc[:, 't01':].mean(axis=1)
    df_out['std'] = df_out.loc[:, 't01':].std(axis=1)
    df_out['median'] = df_out.loc[:, 't01':].median(axis=1)
    quantiles = df_out.apply(lambda x: x.quantile([0.1, 0.3, 0.5, 0.7, 0.9]), axis=1)
    df_out = pd.concat([df_out, quantiles], axis=1)
    if not is_test_running:
        file_out = f'{BASE_PATH}/{algorithm_folder_name[interpolation_algorithm]}/{file_out}'
    else:
        file_out = f'{BASE_PATH}/{algorithm_folder_name[interpolation_algorithm]}/test_dataset_{file_out}'
    df_out.to_csv(file_out, sep='\t', index=True)
    print(f'wrote: {file_out}')

print('Finished...')

