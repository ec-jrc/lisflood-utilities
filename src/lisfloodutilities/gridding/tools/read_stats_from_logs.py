from dask.dataframe.io.tests.test_json import df
__author__="Goncalo Gomes"
__date__="$Jun 06, 2024 10:45:00$"
__version__="0.1"
__updated__="$Jun 06, 2024 10:45:00$"

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
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import json
from lisfloodutilities.gridding.lib.utils import FileUtils


def have_the_same_columns(df1: pd.DataFrame, df2: pd.DataFrame):
    return all(df1.columns.isin(df2.columns))

def run(infolder: str, outfile: str, search_string: str, inwildcard: str = '*.log'):
    out_df = None
    outfilepath = Path(outfile)
    # Create the output parent folders if not exist yet
    Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)

    for filename in sorted(Path(infolder).rglob(inwildcard)):
        print(f'Processing file: {filename}')    
        with open(filename, 'r') as file:
            lines = file.readlines()
        filtered_lines = [line for line in lines if line.startswith(search_string)]
        stats_dictionaries = [json.loads(line.strip()[len(search_string):]) for line in filtered_lines]
        if len(stats_dictionaries) > 0:
            df = pd.DataFrame(stats_dictionaries)
            if out_df is None:
                if not df.empty:
                    out_df = df
            elif not df.empty and have_the_same_columns(out_df, df):
                out_df = pd.concat([out_df], df)
            else:
                print('ERROR: Both datasets do not have the same columns')
    if out_df is None or out_df.empty:
        print('WARNING: No lines containing statistics where found.')
    else:
        out_df = out_df.drop_duplicates()
        out_df.to_csv(outfilepath, index=False, header=True, sep='\t')
        print(f'Wrote file: {outfilepath}')
        print(out_df)


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script parses a list of log files containing statistics in the form of dictionary and converts the logged statistics into one tab separated file.
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
        parser.set_defaults(search_string='#APP_STATS: ',
                            input_wildcard='*.log')

        parser.add_argument("-i", "--in", dest="infolder", required=True, type=FileUtils.folder_type,
                            help="Set input folder path with log files (*.log)",
                            metavar="/input/folder/logfiles/")
        parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_type,
                            help="Set output file name (*.tsv).",
                            metavar="/path/to/output_file.tsv")
        parser.add_argument("-w", "--wildcard", dest="input_wildcard", required=False, type=str,
                            help=('Set the input wildcard to filter out the log files to be processed.'),
                            metavar='*.log')
        parser.add_argument("-s", "--search", dest="search_string", required=False, type=str,
                            help=('Set line tag that identifies the statistics dictionary. '
                                  'It will be used to parse the line, so the space at the end is necessary.'),
                            metavar='"#APP_STATS: "')

        # process options
        args = parser.parse_args(argv)

        print(f"Input Folder: {args.infolder}")
        print(f"Output Filer: {args.outfile}")
        print(f"Input Wildcard: {args.input_wildcard}")
        print(f'Search String: [{args.search_string}]')

        run(args.infolder, args.outfile, args.search_string, args.input_wildcard)
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
