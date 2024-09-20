
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
import csv
from datetime import datetime, timedelta
from lisfloodutilities.gridding.lib.utils import FileUtils


COL_OUTPUT_QUALITY_CODE_WRONG = 'QUALITY_CODE_WRONG'
COL_OUTPUT_TOTAL_OBSERVATIONS = 'TOTAL_OBSERVATIONS'
COL_OUTPUT_TIMESTAMP = 'TIMESTAMP'
COL_OUTPUT_VAR_CODE = 'VAR_CODE'
COL_OUTPUT_PROVIDER_ID = 'PROVIDER_ID'


def run(statfile: str, outfile: str):
    
    outfilepath = Path(outfile)
    # Create the output parent folders if not exist yet
    Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)
    
    statfilepath = Path(statfile)
    print(f'Reading statistics file: {statfilepath}')
    df_stats = pd.read_csv(statfilepath, sep="\t")
    provider_ids = sorted(df_stats[COL_OUTPUT_PROVIDER_ID].unique())
    
    print('provider_ids:', provider_ids)

    first_timestamp_cell = df_stats[COL_OUTPUT_TIMESTAMP].iloc[0]
    yyyy = first_timestamp_cell[:4]
    
    ncols = len(provider_ids)

    out_row1 = [yyyy, 'DP']
    out_row1.extend(provider_ids)
    out_row2 = [yyyy, 'average stations with data']
    out_row2.extend([''] * ncols)
    out_row3 = [yyyy, 'average error']
    out_row3.extend([''] * ncols)
    out_row4 = [yyyy, 'max number of errors in a day']
    out_row4.extend([''] * ncols)

    i = 0
    for provider_id in provider_ids:
        average_stations = df_stats.loc[df_stats[COL_OUTPUT_PROVIDER_ID] == provider_id, COL_OUTPUT_TOTAL_OBSERVATIONS].mean()
        out_row2[2 + i] = round(average_stations,0)
        average_error = df_stats.loc[df_stats[COL_OUTPUT_PROVIDER_ID] == provider_id, COL_OUTPUT_QUALITY_CODE_WRONG].mean()
        out_row3[2 + i] = round(average_error)
        max_error = df_stats.loc[df_stats[COL_OUTPUT_PROVIDER_ID] == provider_id, COL_OUTPUT_QUALITY_CODE_WRONG].max()
        out_row4[2 + i] = round(max_error)
        i += 1

    with open(outfilepath, 'a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(out_row1)
        writer.writerow(out_row2)
        writer.writerow(out_row3)
        writer.writerow(out_row4)

    print(f'Wrote file: {outfilepath}')


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script extracts kiwis logged statistics into another tab separated file.
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

    # try:
    if True:
        # setup option parser
        parser = ArgumentParser(epilog=program_license, description=program_version_string+program_longdesc)

        # set defaults
        # parser.set_defaults(input_wildcard='*.tsv')

        parser.add_argument("-s", "--stat", dest="statfile", required=True, type=FileUtils.file_type,
                            help="Set input file containing kiwis statistics name (*.tsv).",
                            metavar="/path/to/kiwis_stats_ws_2001.tsv")
        parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_type,
                            help="Set output file name (*.tsv).",
                            metavar="/path/to/output_file.tsv")

        # process options
        args = parser.parse_args(argv)

        print(f"Statistics File: {args.statfile}")
        print(f"Output File: {args.outfile}")

        run(args.statfile, args.outfile)
        print("Finished.")
    # except Exception as e:
    #     indent = len(program_name) * " "
    #     sys.stderr.write(program_name + ": " + repr(e) + "\n")
    #     sys.stderr.write(indent + "  for help use --help")
    #     return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()
