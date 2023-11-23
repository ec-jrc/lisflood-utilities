import pandas as pd

variable_name = 'pr'

BASE_PATH = '/BGFS/DISASTER/grids/leave_one_out'

timesteps = ['199106170600', '199106180600']
# timesteps = ['199106180600']
interpolation_algorithms = ['spheremap', 'idw', 'adw', 'cdd']

for timestep in timesteps:
    for interpolation_algorithm in interpolation_algorithms:
        file_in = f'{BASE_PATH}/{interpolation_algorithm}/{variable_name}/{variable_name}{timestep}_{interpolation_algorithm}_BY_CATCHMENT.stats'
        file_out = f'{BASE_PATH}/{interpolation_algorithm}/{variable_name}/{variable_name}{timestep}_{interpolation_algorithm}_BY_CATCHMENT_out.stats'
        df = pd.read_csv(file_in, delimiter='\t')
        df = df.set_index('catchments')
        df['mean'] = df.loc[:, 't01':].mean(axis=1)
        df['std'] = df.loc[:, 't01':].std(axis=1)
        df.to_csv(file_out, sep='\t', index=True)

