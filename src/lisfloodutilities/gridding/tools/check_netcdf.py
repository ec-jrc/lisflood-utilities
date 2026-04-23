from netCDF4 import Dataset
import sys
import argparse

def check_netcdf(filename: str, var_id: str):
    """
    Check if a NetCDF file can be opened and contains grids in all timesteps.

    Parameters:
    - filename (str): The name of the NetCDF file.
    - var_id (str): The current ID of the variable.
    """
    try:
        # Open the NetCDF file in read mode
        with Dataset(filename, 'r') as ds:
            # print(f"File: {filename}")
            for t in range(0, len(ds.variables['time']), 1):
                # print(f'time: {t}')
                cur_grid = ds.variables[var_id][t, :, :]
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
