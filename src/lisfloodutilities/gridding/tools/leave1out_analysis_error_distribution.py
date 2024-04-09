__author__="Goncalo Gomes"
__date__="$Jul 12, 2023 12:01:00$"
__version__="0.1"
__updated__="$Jul 13, 2023 10:41:00$"

"""
Copyright 2019-2020 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import sys
import os
import csv
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError
from lisfloodutilities.gridding.lib.utils import FileUtils
import pandas as pd
import numpy as np
import rioxarray as rxr
from osgeo import gdal
import netCDF4 as nc
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import mean_absolute_error
from sklearn.feature_selection import r_regression


def ratio(y_true: np.array, y_pred: np.array) -> np.array:
    '''
    Parameters:
        y_true (array): Array of observed values
        y_pred (array): Array of prediction values

    Returns:
        Ratio array between observed and predicted values
    '''
    return np.abs(1 - y_pred / y_true)


def write_results(outfile: Path, input_df: pd.DataFrame):
    if outfile.is_file():
        output_df = pd.read_csv(outfile, delimiter='\t')
        df = pd.merge(output_df, input_df, on=['x', 'y'], suffixes=('_1', '_2'), how='outer')
        df['ratio'] = df['ratio_1'].fillna(0) + df['ratio_2'].fillna(0)
        df.drop(['ratio_1', 'ratio_2'], axis=1, inplace=True)
    else:
        df = pd.DataFrame({
            'x': input_df['x'],
            'y': input_df['y'],
            'ratio': input_df['ratio'],
        })
    df.to_csv(outfile, sep='\t', index=False)


def run(file_true_values: str, outfile: str):
    outfile_path = Path(outfile)
    file_true_values_path = Path(file_true_values)
    # Load the CSV file into a pandas dataframe
    df_true_and_predicted_values = pd.read_csv(file_true_values_path.with_suffix('.tab'), delimiter='\t')
    # remove all duplicates by x,y since we do not know which to keep and they give issues writing the results
    df_true_and_predicted_values = df_true_and_predicted_values.drop_duplicates(['x', 'y'], keep=False)

    limit_value = 0.0
    df_true_and_predicted_values = df_true_and_predicted_values[(df_true_and_predicted_values['true_value'] > limit_value) & (df_true_and_predicted_values['interpolated_value'] > limit_value)]

    true_values = df_true_and_predicted_values['true_value']
    predicted_values = df_true_and_predicted_values['interpolated_value']

    df_true_and_predicted_values['ratio'] = ratio(true_values, predicted_values)
    # leaving only x, y, ratio
    df_true_and_predicted_values.drop(['true_value', 'interpolated_value'], axis=1, inplace=True)
    # write the results to the output file
    write_results(outfile_path, df_true_and_predicted_values)


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script uses RepeatedStratifiedKFold to generate several input files as training dataset where a small percentage
    of evenly distributed stations are removed to constitute the test cases. Later the training dataset will be interpolated
    with several algorithms: Spheremap, Inverse distance and Angular Distance Weighting. Then the error at the coordinates of
    the test cases is calculated by the difference between the interpolated value and the real value in the test case.
    '''
    program_license = """
    Copyright 2019-2020 European Union
    Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
    You may not use this work except in compliance with the Licence.
    You may obtain a copy of the Licence at:
    https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
    Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the Licence for the specific language governing permissions and limitations under the Licence.
    """

    try:
        # setup option parser
        parser = ArgumentParser(epilog=program_license, description=program_version_string+program_longdesc)

        # set defaults
        parser.set_defaults(is_compressed_data=False)

        parser.add_argument("-t", "--test", dest="file_true_values", required=True, type=FileUtils.file_type,
                            help="Set file path for true values form the test dataset in txt format.",
                            metavar="txt_file")
        parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_or_folder,
                            help="Set output file containing the statistics.",
                            metavar="output_file")

        # process options
        args = parser.parse_args(argv)

        print(f"True values: {args.file_true_values}")
        print(f"Output File: {args.outfile}")

        run(args.file_true_values, args.outfile)
        print("Finished.")
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()