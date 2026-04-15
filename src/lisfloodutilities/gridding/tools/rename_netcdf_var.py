from netCDF4 import Dataset
import sys
import argparse

def rename_netcdf_variable(filename: str, var_id: str, new_var_id: str):
    """
    Rename a variable in a NetCDF file.

    Parameters:
    - filename (str): The name of the NetCDF file.
    - var_id (str): The current ID of the variable.
    - new_var_id (str): The new ID for the variable.
    """
    try:
        # Open the NetCDF file in write mode
        with Dataset(filename, 'r+') as ds:
            # Rename the variable
            ds.renameVariable(var_id, new_var_id)
            print(f"Variable '{var_id}' renamed to '{new_var_id}' in {filename}")
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    parser = argparse.ArgumentParser(description="Rename a variable in a NetCDF file.")
    parser.add_argument("filename", help="The name of the NetCDF file.")
    parser.add_argument("-o", "--old-id", required=True, help="The current ID of the variable.")
    parser.add_argument("-n", "--new-id", required=True, help="The new ID for the variable.")
    args = parser.parse_args()

    rename_netcdf_variable(args.filename, args.old_id, args.new_id)

if __name__ == "__main__":
    main()
