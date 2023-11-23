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
import pandas as pd
from lisfloodutilities.gridding.lib.utils import FileUtils
from sklearn.model_selection import RepeatedStratifiedKFold


def write_datasets(i: int, orig_outfilepath: Path, X_df: pd.DataFrame, y_df: pd.DataFrame, write_train_dataset=False):
    if i > 9:
        test_folder = Path(orig_outfilepath.parent).joinpath(f't{i}')
    else:
        test_folder = Path(orig_outfilepath.parent).joinpath(f't0{i}')
    output_kiwis_filename = orig_outfilepath.name
    kiwis_filepath = Path(test_folder).joinpath(output_kiwis_filename)
    Path(kiwis_filepath.parent).mkdir(parents=True, exist_ok=True)
    # Create an empty kiwis file once to be used later by the spheremap program to detect the files to be processed
    if write_train_dataset:
        with open(kiwis_filepath, mode='a'): pass
        output_filename = output_kiwis_filename.replace('all.kiwis', '20230714101901.txt')
    else:
        # write test datset
        output_filename = output_kiwis_filename.replace('all.kiwis', 'TEST_DATASET.txt')
    output_file = test_folder.joinpath(output_filename)
    df2 = pd.DataFrame({
        'x': X_df['station_local_x'].values,
        'y': X_df['station_local_y'].values,
        'value': y_df['ts_value'].values,
    })
    df2.to_csv(output_file, sep='\t', header=False, index=False)


def run(infolder: str, outfolder: str):
    inwildcard = '*_all.kiwis'

    for filename in sorted(Path(infolder).rglob(inwildcard)):
        print(f'Processing file: {filename}')
        outfile = str(filename).replace(infolder, outfolder)
        outfilepath = Path(outfile)
        # Create the output parent folders if not exist yet
        Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)

        # Load the CSV file into a pandas dataframe
        df = pd.read_csv(filename, delimiter='\t')

        # Clean the dataset from unwanted rows
        df = df.astype(str)
        df = df[df['EFAS-ADDATTR-NOGRIDDING'].isin(['no'])]
        df = df[df['EFAS_ADDATTR_ISINNEWDOMAIN'].isin(['yes'])]
        df = df[~df['EXCLUDE'].isin(['yes'])]
        df = df[df['q_code'].isin(['40', '120'])]

        # Use the x and y coordinate columns as our input features (X),
        # and the value column as our target variable (y)
        X = df[['station_local_x', 'station_local_y']]
        # y = df[['site_no', 'ts_value']]
        y = df[['ts_value']]

        # Instantiate the RepeatedStratifiedKFold object
        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=47)

        i = 0
        # Iterate over the folds created by RepeatedStratifiedKFold
        for train_index, test_index in rskf.split(X, y):
            i += 1
            # Split the data into training and testing sets based on the current fold indices
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]

            # Write the train dataset
            X_train_df = X_train[['station_local_x', 'station_local_y']]
            y_train_df = y_train[['ts_value']]
            write_datasets(i, outfilepath, X_train_df, y_train_df, write_train_dataset=True)

            # Write the test dataset
            X_test_df = X_test[['station_local_x', 'station_local_y']]
            y_test_df = y_test[['ts_value']]
            write_datasets(i, outfilepath, X_test_df, y_test_df, write_train_dataset=False)


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

#        # set defaults
#        parser.set_defaults(quiet=False,
#                            out_tiff=False,
#                            overwrite_output=False,
#                            start_date='',
#                            end_date=END_DATE_DEFAULT)

        parser.add_argument("-i", "--in", dest="infolder", required=True, type=FileUtils.folder_type,
                            help="Set input folder path with kiwis files",
                            metavar="input_folder")
        parser.add_argument("-o", "--out", dest="outfolder", required=True, type=FileUtils.folder_type,
                            help="Set output folder base path.",
                            metavar="output_folder")

        # process options
        args = parser.parse_args(argv)

        print(f"Input Folder: {args.infolder}")
        print(f"Output Folder: {args.outfolder}")

        run(args.infolder, args.outfolder)
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
