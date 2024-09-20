
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
from datetime import datetime, timedelta
from lisfloodutilities.gridding.lib.utils import Config, FileUtils


COL_PROVIDER_ID = 'site_no'
COL_QUALITY_CODE = 'q_code'
QUALITY_CODE_WRONG = 160
COL_OUTPUT_QUALITY_CODE_WRONG = 'QUALITY_CODE_WRONG'
COL_OUTPUT_TOTAL_OBSERVATIONS = 'TOTAL_OBSERVATIONS'
COL_OUTPUT_TIMESTAMP = 'TIMESTAMP'
COL_OUTPUT_VAR_CODE = 'VAR_CODE'
COL_OUTPUT_PROVIDER_ID = 'PROVIDER_ID'


def set_values(row: pd.Series, timestamp: str, var_code: str, provider_id: str, total_wrong_values: int) -> pd.Series:
    if str(row[COL_OUTPUT_TIMESTAMP]) == timestamp and str(row[COL_OUTPUT_VAR_CODE]) == var_code and str(row[COL_OUTPUT_PROVIDER_ID]) == provider_id:
        row[COL_OUTPUT_QUALITY_CODE_WRONG] += total_wrong_values
        row[COL_OUTPUT_TOTAL_OBSERVATIONS] += total_wrong_values
    return row

def get_timestamp_from_filename(conf: Config, filename: Path) -> datetime:
    file_timestamp = datetime.strptime(filename.name, conf.input_timestamp_pattern)
    if conf.force_time is not None:
        new_time = datetime.strptime(conf.force_time, "%H%M").time()
        file_timestamp = datetime.combine(file_timestamp.date(), new_time)
    return file_timestamp

def run(config_filename: str, infolder: str, variable_code: str, statfile: str, outfile: str):
    
    conf = Config(config_filename)

    INPUT_TIMESTAMP_PATTERN = conf.input_timestamp_pattern
    inwildcard = conf.input_wildcard
    
    netcdf_offset_file_date = int(conf.get_config_field('VAR_TIME','OFFSET_FILE_DATE'))

    outfilepath = Path(outfile)
    # Create the output parent folders if not exist yet
    Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)
    
    statfilepath = Path(statfile)
    print(f'Reading statistics file: {statfilepath}')
    df_stats = pd.read_csv(statfilepath, sep="\t")

    for filename in sorted(Path(infolder).rglob(inwildcard)):
        print(f'Processing file: {filename}')
        df_kiwis = pd.read_csv(filename, sep="\t")
        file_timestamp = get_timestamp_from_filename(conf, filename)
        file_timestamp = file_timestamp + timedelta(days=netcdf_offset_file_date)
        timestamp = file_timestamp.strftime('%Y-%m-%d %H:%M:%S')
        new_df = df_kiwis.groupby([COL_PROVIDER_ID, COL_QUALITY_CODE]).size().reset_index(name='count')
        new_df = new_df.loc[new_df[COL_QUALITY_CODE] == QUALITY_CODE_WRONG]
        if not new_df.empty:
            for index, row in new_df.iterrows():
                provider_id = str(row[COL_PROVIDER_ID])
                total_wrong_values = int(row['count'])
                df_stats = df_stats.apply(set_values, axis=1, timestamp=timestamp, var_code=variable_code,
                                          provider_id=provider_id, total_wrong_values=total_wrong_values)
    df_stats.to_csv(outfilepath, index=False, header=True, sep='\t')
    print(f'Wrote file: {outfilepath}')
    print(df_stats)


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
        # parser.set_defaults(input_wildcard='*.kiwi')

        parser.add_argument("-i", "--in", dest="infolder", required=True, type=FileUtils.folder_type,
                            help="Set input folder path with kiwis files (*.kiwi)",
                            metavar="/input/folder/kiwis/")
        parser.add_argument("-s", "--stat", dest="statfile", required=True, type=FileUtils.file_type,
                            help="Set input file containing kiwis statistics name (*.tsv).",
                            metavar="/path/to/kiwis_stats_file.tsv")
        parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_type,
                            help="Set output file name (*.tsv).",
                            metavar="/path/to/output_file.tsv")
        parser.add_argument("-v", "--var", dest="variable_code", required=True,
                            help="Set the variable to be processed.",
                            metavar="{pr,pd,tn,tx,ws,rg,...}")
        parser.add_argument("-c", "--conf", dest="config_type", required=True,
                            help="Set the grid configuration type to use.",
                            metavar="{5x5km, 1arcmin,...}")
        parser.add_argument("-p", "--pathconf", dest="config_base_path", required=False, type=FileUtils.folder_type,
                            help="Overrides the base path where the configurations are stored.",
                            metavar="/path/to/config")

        # process options
        args = parser.parse_args(argv)

        configuration_base_folder = os.path.join(program_path, '../../src/lisfloodutilities/gridding/configuration')
        if args.config_base_path is not None and len(args.config_base_path) > 0:
            configuration_base_folder = args.config_base_path

        file_utils = FileUtils(args.variable_code)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, args.config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        print(f"Variable: {args.variable_code}")
        print(f"Input Folder: {args.infolder}")
        print(f"Statistics File: {args.statfile}")
        print(f"Output File: {args.outfile}")
        print(f"Config File: {config_filename}")

        run(config_filename, args.infolder, args.variable_code, args.statfile, args.outfile)
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
