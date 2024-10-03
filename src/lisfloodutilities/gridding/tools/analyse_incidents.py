
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

DELIMITER_INPUT = ';'
DELIMITER_OUTPUT = '\t'

COL_PROVIDER_ID = 'SITE'
COL_STATION_NUM = 'STATION'
COL_PARAMETER = 'PARAM'
COL_TIMESERIES = 'TS'
COL_STATION_NAME = 'NAME'
COL_STATUS = 'Status'
COL_LAT = 'latitude'
COL_LON = 'longitude'
COL_HEIGHT = 'elevation'
COL_RELATED_PATHS = 'RelatedPaths'
COL_RESULTS_TYPE = 'resultsType'
COL_MSG = 'message'
COL_TIME_INTERVAL_START = 'from'
COL_TIME_INTERVAL_END = 'until'
COL_NUM_INVALID = 'nInvalid'
COL_NUM_VALID = 'nValid'
COL_PERCENTAGE_INVALID = 'percInvalid'
COL_INCIDENTS = 'incidents'
COL_TOTAL_INCIDENTS = 'totalIncidents'


incident_type_columns = {}


def get_total_incidents(row: pd.Series) -> int:
    global incident_type_columns
    incidents = row[COL_INCIDENTS]
    if incidents is None:
        return 0
    incidents = incidents.strip()
    if len(incidents) == 0 or incidents == 'nan' or incidents == '{}':
        return 0
    incidents_dic = eval(incidents)
    if not isinstance(incidents_dic, dict):
        return 0
    total_incidents = 0
    for incident_key in incidents_dic:
        # Used to store all the incident types available
        incident_type_columns[incident_key] = ''
        try:
            total_incidents += int(incidents_dic[incident_key])
        except Exception as e:
            print(f'ERROR evaluating row: {row}')
    return total_incidents


def get_total_incident_by_type(row: pd.Series, incident_type: str) -> int:
    incidents = row[COL_INCIDENTS]
    if incidents is None:
        return 0
    incidents = incidents.strip()
    if len(incidents) == 0 or incidents == 'nan' or incidents == '{}':
        return 0
    incidents_dic = eval(incidents)
    if not isinstance(incidents_dic, dict):
        return 0
    total_incidents = 0
    for incident_key in incidents_dic:
        if incident_key == incident_type:
            try:
                total_incidents = int(incidents_dic[incident_key])
            except Exception as e:
                print(f'ERROR evaluating row: {row}')
    return total_incidents


def run(infolder: str, outfolder: str):
    global incident_type_columns
    inwildcard = '*.csv'
    
    for filename in sorted(Path(infolder).rglob(inwildcard)):
        print(f'Processing file: {filename}')
        outfile = f'{filename}_out.tsv'.replace(infolder, outfolder)
        outfilepath = Path(outfile)
        # Create the output parent folders if not exist yet
        Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(filename, delimiter=DELIMITER_INPUT, low_memory=False)
        df = df.astype({COL_PROVIDER_ID: 'str', COL_INCIDENTS: 'str'})
        incident_type_columns = {}
        df[COL_TOTAL_INCIDENTS] = df.apply(get_total_incidents, axis=1)

        # define aggregation functions
        agg_funcs = {COL_TOTAL_INCIDENTS: [('Total Incidents', 'sum'), ('Number of Stations', 'count')]}

        incident_types = list(incident_type_columns.keys())
        for incident_type in incident_types:
            agg_funcs[incident_type] = [(incident_type, 'sum')]
            df[incident_type] = None
            df[incident_type] = df.apply(get_total_incident_by_type, axis=1, incident_type=incident_type)

        groupy_cols = [COL_PROVIDER_ID, COL_PARAMETER, COL_TIMESERIES]
        df = df.groupby(groupy_cols).agg(agg_funcs).reset_index()

        # Eliminate the top level of column names since the new names are
        # written on the bottom level and insert the name of the the
        # aggregation column on the bottom level
        df.columns = df.columns.droplevel(0)
        df.reset_index(drop=True, inplace=True)
        columns = list(df.columns)
        columns[0] = 'Provider'
        columns[1] = COL_PARAMETER
        columns[2] = COL_TIMESERIES
        df.columns = columns
        
        df['Incidents per Station'] = df['Total Incidents'].div(df['Number of Stations'])
        df['Incidents per Station'] = df['Incidents per Station'].round(2)
        
        df.sort_values(by=['Incidents per Station', 'Total Incidents', 'Number of Stations'], ascending=[False, False, False], inplace=True)

        if df.empty:
            print(f'WARNING: No data was found in file {filename}')
        else:
            df.to_csv(outfilepath, index=False, header=True, sep=DELIMITER_OUTPUT)
            print(f'Wrote file: {outfilepath}')
            # print(out_df)


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script parses a list of CSV files containing KIWIS incidents and analyses them to produce a report into tab separated file for each.
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

        # # set defaults
        # parser.set_defaults(search_string='#APP_STATS: ')

        parser.add_argument("-i", "--in", dest="infolder", required=True, type=FileUtils.folder_type,
                            help="Set input folder path with log files (*.csv)",
                            metavar="/input/folder/logfiles/")
        parser.add_argument("-o", "--out", dest="outfolder", required=True, type=FileUtils.folder_type,
                            help="Set output folder where the analysis files will be stored (*_out.tsv).",
                            metavar="/path/to/output_file.tsv")

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
