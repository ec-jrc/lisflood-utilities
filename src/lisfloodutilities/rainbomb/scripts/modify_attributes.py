from netCDF4 import Dataset

def getarg():
    """
    Get program arguments.

    :return: args:  namespace of program arguments
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--inputfile",
        type=str,
        required=True,
        help="Path of input file for amending attributes",
    )
    parser.add_argument(
        "-y",
        "--year",
        type=str,
        required=True,
        help="Year for defining start of time",
    )
    parser.add_argument(
        "-v",
        "--variable",
        type=str,
        required=True,
        help="Variable name",
    )
    parser.add_argument(
        "-u",
        "--units",
        type=str,
        required=True,
        help="Units of the variable",
    )
    args = parser.parse_args()  # assign namespace to args
    return args


def main():
    """Function for running the script as main"""

    # ----------- Read arguments -----------
    args = getarg()
    input_file = args.inputfile
    year = args.year
    variable_name = args.variable
    units = args.units


    nf2 = Dataset(input_file, 'r+', format='NETCDF4_CLASSIC') 
    nf2.variables['time'].units = f'days since {year}-01-01 00:00:00'
    tp = nf2.variables[variable_name]
    tp.units = units
    
    nf2.close()


if __name__ == "__main__":
    main()