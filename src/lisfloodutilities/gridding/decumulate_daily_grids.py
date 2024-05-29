__author__="Goncalo Gomes"
__date__="$May 14, 2024 12:01:00$"
__version__="0.1"
__updated__="$Mar 14, 2024 16:01:00$"

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
from typing import List, Tuple
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError
from datetime import datetime, timedelta
import numpy as np
import numpy.ma as ma
import pandas as pd
from collections import OrderedDict
from scipy.spatial import cKDTree
from lisfloodutilities.gridding.lib.utils import Config, FileUtils
from lisfloodutilities.gridding.lib.filters import KiwisFilter, DowgradedDailyTo6HourlyObservationsKiwisFilter

# Disable the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

quiet_mode = False


def print_msg(msg: str = ''):
    global quiet_mode
    if not quiet_mode:
        print(msg)


def get_filter_class(conf: Config) -> KiwisFilter:
    # Allow us to get the output columns out of the dataframe
    plugins_columns_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
    plugins_columns_dic = eval(plugins_columns_def)
    plugins_columns = OrderedDict(plugins_columns_dic)
    filter_class = KiwisFilter(filter_columns=plugins_columns)
    return filter_class


def get_timestamp_from_filename(conf: Config, filename: Path) -> datetime:
    file_timestamp = datetime.strptime(filename.name, conf.input_timestamp_pattern)
    if conf.force_time is not None:
        new_time = datetime.strptime(conf.force_time, "%H%M").time()
        file_timestamp = datetime.combine(file_timestamp.date(), new_time)
    return file_timestamp


def get_dataframes(conf: Config, kiwis_files: List[Path]) -> Tuple[List[pd.DataFrame], List[datetime], List[str]]:
    filter_class = get_filter_class(conf)
    print_msg(f'Filter class: {filter_class.get_class_description()}')
    filepath_kiwis = []
    kiwis_timestamps = []
    kiwis_str_timestamps = []
    for filename_kiwis in kiwis_files:
        kiwis_timestamp = get_timestamp_from_filename(conf, filename_kiwis)
        kiwis_timestamp_str = kiwis_timestamp.strftime(FileUtils.DATE_PATTERN_CONDENSED_SHORT)
        filepath_kiwis.append(filename_kiwis)
        kiwis_str_timestamps.append(kiwis_timestamp_str)
        kiwis_timestamps.append(kiwis_timestamp)
    df_kiwis_array = filter_class.filter(filepath_kiwis, kiwis_str_timestamps, [])
    return df_kiwis_array, kiwis_timestamps, kiwis_str_timestamps

def get_providers_radius(conf: Config, kiwis_timestamp_24h: datetime) -> dict:
    """
    Get the configuration of providers and radius correspondent to the current daily date
    """
    TIMESTAMP_PATTERN = '%Y-%m-%d %H:%M:%S'
    filter_args = {}
    plugin_decumulation_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_DECUMULATION_CONFIG')
    plugin_decumulation_config = eval(plugin_decumulation_def)
    
    for cur_config in plugin_decumulation_config:
        try:
            start_timestamp = datetime.strptime(cur_config['START'], TIMESTAMP_PATTERN)
        except Exception as e:
            start_timestamp = None
        try:
            end_timestamp = datetime.strptime(cur_config['END'], TIMESTAMP_PATTERN)
        except Exception as e:
            end_timestamp = None
        providers_radius = cur_config['PROVIDERS_RADIUS']
        if ((start_timestamp is None or kiwis_timestamp_24h >= start_timestamp) and
            (end_timestamp is None or kiwis_timestamp_24h <= end_timestamp)):
            filter_args = providers_radius
    return filter_args

def get_decumulated_dataframes(conf: Config, kiwis_filepaths: List[Path], kiwis_str_timestamps: List[str],
                               kiwis_dataframes: List[pd.DataFrame], kiwis_timestamp_24h: datetime) -> List[pd.DataFrame]:
    """
    Applies the decumulation filter to the dataframes. The first dataframe in the array should be the Daily
    one and then the remaining 4 dataframes will be the 6 hourly ones.
    """

    filter_args = get_providers_radius(conf, kiwis_timestamp_24h)

    # If there was no data interval defined for the current daily date, then no need to decumulate
    if not filter_args:
        return kiwis_dataframes
    
    plugins_columns_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
    plugins_columns_dic = eval(plugins_columns_def)
    filter_columns = OrderedDict(plugins_columns_dic)
    filter_class = DowgradedDailyTo6HourlyObservationsKiwisFilter(filter_columns, filter_args)
    df_kiwis_array = filter_class.filter(kiwis_filepaths, kiwis_str_timestamps, kiwis_dataframes)
    return df_kiwis_array


def run(config_filename_24h: str, config_filename_6h: str, file_utils_24h: FileUtils, file_utils_6h: FileUtils,
        kiwis_24h_06am_path: Path, kiwis_6h_12pm_path: Path, kiwis_6h_18pm_path: Path, kiwis_6h_12am_path: Path, kiwis_6h_06am_path: Path,
        output_path: Path = None):
    """
    Interpolate text files containing (x, y, value) using inverse distance interpolation.
    Produces as output, either a netCDF file containing all the grids or one TIFF file per grid.
    1. Read the config files
    2. Get the ordered list of files
    3. Interpolate the grids performing height correction on some variables
    4. Write the resulting grids
    """
    global quiet_mode

    interpolation_mode = 'adw'
    memory_save_mode = 0
    start_date = None
    end_date = None
    conf_24h = Config(config_filename_24h, start_date, end_date, quiet_mode, interpolation_mode, memory_save_mode)
    conf_6h = Config(config_filename_6h, start_date, end_date, quiet_mode, interpolation_mode, memory_save_mode)

    print_msg('Start reading files')
    df_kiwis_array_24h, kiwis_timestamps_24h, kiwis_str_timestamps_24h = get_dataframes(conf_24h, [kiwis_24h_06am_path])
    df_kiwis_array_6h, kiwis_timestamps_6h, kiwis_str_timestamps_6h = get_dataframes(conf_6h, [kiwis_6h_12pm_path, kiwis_6h_18pm_path,
                                                                                               kiwis_6h_12am_path, kiwis_6h_06am_path])
    # Check timestamps are correct
    if not (kiwis_timestamps_24h[0] == kiwis_timestamps_6h[3] and
            (kiwis_timestamps_24h[0] - timedelta(hours=6)) == kiwis_timestamps_6h[2] and
            (kiwis_timestamps_24h[0] - timedelta(hours=12)) == kiwis_timestamps_6h[1] and
            (kiwis_timestamps_24h[0] - timedelta(hours=18)) == kiwis_timestamps_6h[0]):
        raise ArgumentTypeError("The input kiwis do not respect the expected timestamps.")
    
    kiwis_dataframes = df_kiwis_array_24h
    kiwis_dataframes.extend(df_kiwis_array_6h)
    kiwis_timestamps = kiwis_str_timestamps_24h
    kiwis_timestamps.extend(kiwis_str_timestamps_6h)
    kiwis_filepaths = [kiwis_24h_06am_path, kiwis_6h_12pm_path, kiwis_6h_18pm_path, kiwis_6h_12am_path, kiwis_6h_06am_path]

    df_kiwis_array = get_decumulated_dataframes(conf_24h, kiwis_filepaths, kiwis_timestamps, kiwis_dataframes, kiwis_timestamps_24h[0])

    i = 0
    for kiwis_filepath in kiwis_filepaths[1:]:
        i += 1
        if output_path is not None:
            filepath = Path.joinpath(output_path, kiwis_filepath.name)
        else:
            filepath = kiwis_filepath
        df_kiwis_array[i].to_csv(filepath, index=False, header=True, sep="\t")
        print_msg(f'Wrote file: {filepath}')
    print_msg('Finished writing files')

def get_existing_file_path(parser: ArgumentParser, input_file_path: str) -> Path:
    file_path = Path(input_file_path)
    if not file_path.is_file():
        parser.error(f'Input file {file_path} does not exist or cannot be opened.')
    return file_path

def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script performs the decumulation of daily kiwis files into the respective 6hourly precipitation kiwis files in order to increase the station density or to fill in 6hourly gaps.
    JIRA Issue: EMDCC-1484
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
        parser.set_defaults(quiet=False,
                            out_folder=None)

        parser.add_argument("-d", "--pr24h", dest="kiwis_24h_06am_path", required=True, type=FileUtils.file_type,
                            help="Set the input kiwis file containing daily precipitation for a given day.",
                            metavar="/path/to/pr200102150600_all.kiwis")
        parser.add_argument("-1", "--pr6h12pm", dest="kiwis_6h_12pm_path", required=True, type=FileUtils.file_type,
                            help="Set the input kiwis file containing the first 6hourly timestep of precipitation 12:00 of the previous day.",
                            metavar="/path/to/pr6200102141200_all.kiwis")
        parser.add_argument("-2", "--pr6h18pm", dest="kiwis_6h_18pm_path", required=True, type=FileUtils.file_type,
                            help="Set the input kiwis file containing the second 6hourly timestep of precipitation 18:00 of the previous day.",
                            metavar="/path/to/pr6200102141800_all.kiwis")
        parser.add_argument("-3", "--pr6h12am", dest="kiwis_6h_12am_path", required=True, type=FileUtils.file_type,
                            help="Set the input kiwis file containing the first 6hourly timestep of precipitation 00:00 of the daily precipitation day.",
                            metavar="/path/to/pr6200102150000_all.kiwis")
        parser.add_argument("-4", "--pr6h06am", dest="kiwis_6h_06am_path", required=True, type=FileUtils.file_type,
                            help="Set the input kiwis file containing the first 6hourly timestep of precipitation 06:00 of the daily precipitation day.",
                            metavar="/path/to/pr6200102150600_all.kiwis")
        parser.add_argument("-o", "--out", dest="output_folder", required=False, type=FileUtils.folder_type,
                            help="Set the output folder where the resulting 6h kiwis will be stored. If this folder is not set, it will change the input kiwis.",
                            metavar="/path/to/output/folder")
        parser.add_argument("-c", "--conf", dest="config_type", required=True,
                            help="Set the grid configuration type to use.",
                            metavar="{5x5km, 1arcmin,...}")
        parser.add_argument("-p", "--pathconf", dest="config_base_path", required=False, type=FileUtils.folder_type,
                            help="Overrides the base path where the configurations are stored.",
                            metavar="/path/to/config")
        parser.add_argument("-v", "--var24h", dest="variable_code_24h", required=True,
                            help="Set the daily variable to be processed.",
                            metavar="{pr,ta,...}")
        parser.add_argument("-6", "--var6h", dest="variable_code_6h", required=True,
                            help="Set the 6hourly variable to be processed.",
                            metavar="{pr6,ta6,...}")
        parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Set script output into quiet mode [default: %(default)s]")

        # process options
        args = parser.parse_args(argv)
        quiet_mode = args.quiet

        configuration_base_folder = os.path.join(program_path, '../src/lisfloodutilities/gridding/configuration')

        if args.config_base_path is not None and len(args.config_base_path) > 0:
            configuration_base_folder = args.config_base_path

        file_utils_24h = FileUtils(args.variable_code_24h, quiet_mode)
        config_type_path_24h = file_utils_24h.get_config_type_path(configuration_base_folder, args.config_type)
        config_filename_24h = file_utils_24h.get_config_file(config_type_path_24h)

        file_utils_6h = FileUtils(args.variable_code_6h, quiet_mode)
        config_type_path_6h = file_utils_6h.get_config_type_path(configuration_base_folder, args.config_type)
        config_filename_6h = file_utils_6h.get_config_file(config_type_path_6h)

        if args.output_folder is not None and len(args.output_folder) > 0:
            output_path = Path(args.output_folder)
            print_msg(f"Output Folder: {output_path}")
        else:
            output_path = None
            print_msg(f"Output Folder: 6hourly kiwis files will be overwritten")

        kiwis_24h_06am_path = get_existing_file_path(parser, args.kiwis_24h_06am_path)
        kiwis_6h_12pm_path = get_existing_file_path(parser, args.kiwis_6h_12pm_path)
        kiwis_6h_18pm_path = get_existing_file_path(parser, args.kiwis_6h_18pm_path)
        kiwis_6h_12am_path = get_existing_file_path(parser, args.kiwis_6h_12am_path)
        kiwis_6h_06am_path = get_existing_file_path(parser, args.kiwis_6h_06am_path)
        
        print_msg(f"Daily PR kiwis file:  {args.kiwis_24h_06am_path}")
        print_msg(f"6hourly PR kiwis file 12:00:  {args.kiwis_6h_12pm_path}")
        print_msg(f"6hourly PR kiwis file 18:00:  {args.kiwis_6h_18pm_path}")
        print_msg(f"6hourly PR kiwis file 00:00:  {args.kiwis_6h_12am_path}")
        print_msg(f"6hourly PR kiwis file 06:00:  {args.kiwis_6h_06am_path}")
        print_msg(f"Config File 24h: {config_filename_24h}")
        print_msg(f"Config File 6h: {config_filename_6h}")

        run(config_filename_24h, config_filename_6h, file_utils_24h, file_utils_6h,
            kiwis_24h_06am_path, kiwis_6h_12pm_path, kiwis_6h_18pm_path,
            kiwis_6h_12am_path, kiwis_6h_06am_path, output_path=output_path)
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()

