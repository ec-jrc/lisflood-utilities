from netCDF4 import Dataset, num2date
import numpy as np
import sys
import argparse

VALUE_NAN = -9999.0
min_values = {'pd': 0.0, 'tn': -55.0, 'tx': -55.0, 'ws': 0.0, 'rg': 0.0, 'pr': 0.0, 'pr6': 0.0, 'ta6': -55.0, 'ta': -55.0, 'e0': 0.0, 'es': 0.0, 'et': 0.0}
max_values = {'pd': 50.0, 'tn': 50.0, 'tx': 55.0, 'ws': 50.0, 'rg': 115000000.0, 'pr': 600.0, 'pr6': 600.0, 'ta6': 55.0, 'ta': 55.0, 'e0': 60.0, 'es': 60.0, 'et': 60.0}

def setNaN(value, var_id, defaultNaN=np.nan):
    try:
        value[value==1e31] = defaultNaN
    except Exception as e:
        print(f"value==1e31 : {str(e)}")
    try:
        value[value==VALUE_NAN] = defaultNaN
    except Exception as e:
        print(f"value=={VALUE_NAN} : {str(e)}")
    try:
        value[value==-32768.0] = defaultNaN
    except Exception as e:
        print(f"value==-32768.0 : {str(e)}")
    try:
        value[value==31082] = defaultNaN
    except Exception as e:
        print(f"value==31082 : {str(e)}")
    value[value < min_values[var_id]] = defaultNaN
    value[value > max_values[var_id]] = defaultNaN
    return value

def print_grid_statistics(current_timestamp, var_code, grid: np.ndarray):
    grid_min = np.nanmin(grid)
    grid_max = np.nanmax(grid)
    grid_mean = np.nanmean(grid)
    grid_percentile_10 = np.nanpercentile(grid, 10)
    grid_percentile_90 = np.nanpercentile(grid, 90)
    stats_string = (
        f'#APP_STATS: {{"TIMESTAMP": "{current_timestamp}", "VAR_CODE": "{var_code}", '
        f'"MINIMUM_VALUE": {grid_min:.2f}, "MAXIMUM_VALUE": {grid_max:.2f}, '
        f'"MEAN_VALUE": {grid_mean:.2f}, "PERCENTILE_10": {grid_percentile_10:.2f}, '
        f'"PERCENTILE_90": {grid_percentile_90:.2f}}}'
    )
    print(stats_string)

def check_netcdf(filename: str, var_id: str):
    """
    Extract statistics from a NetCDF file and print them to the screen.

    Parameters:
    - filename (str): The name of the NetCDF file.
    - var_id (str): The current ID of the variable.
    """
    try:
        # Open the NetCDF file in read mode
        with Dataset(filename, 'r') as ds:
            # print(f"File: {filename}")
            time_var = ds.variables['time']
            time_units = time_var.units
            time_calendar = time_var.calendar
            for t in range(0, len(time_var), 1):
                # print(f'time: {t}')
                cur_grid = ds.variables[var_id][t, :, :]
                time_value = time_var[t]
                time_datetime = num2date(time_value, time_units, calendar=time_calendar)
                current_timestamp = time_datetime.strftime('%Y-%m-%d %H:%M:%S')
                print_grid_statistics(current_timestamp, var_id, setNaN(cur_grid, var_id))
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred in file {filename}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Check a variable in a NetCDF file.")
    parser.add_argument("filename", help="The name of the NetCDF file.")
    parser.add_argument("-v", "--var", required=True, help="The ID of the variable.")
    args = parser.parse_args()

    check_netcdf(args.filename, args.var)

if __name__ == "__main__":
    main()
