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
import importlib
from lisfloodutilities.gridding.lib.utils import Config, FileUtils
from lisfloodutilities.gridding.lib.filters import KiwisFilter, DowgradedDailyTo6HourlyObservationsKiwisFilter

# Disable the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

quiet_mode = False

TIMESTAMP_PATTERN = '%Y-%m-%d %H:%M:%S'


def print_msg(msg: str = ''):
    global quiet_mode
    if not quiet_mode:
        print(msg)


def get_filter_class(conf: Config) -> KiwisFilter:
    # Allow us to get the output columns out of the dataframe
    plugins_columns_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
    plugins_columns_dic = eval(plugins_columns_def)
    plugins_columns = OrderedDict(plugins_columns_dic)
    filter_class = KiwisFilter(filter_columns=plugins_columns, var_code=conf.var_code, quiet_mode=quiet_mode)
    return filter_class


def get_filter_classes(conf: Config) -> list:
    '''
    Create and instantiate the Filter classes to process the kiwis files
    '''
    global quiet_mode
    plugins_array = []
    # Load the class dynamically
    plugins_config = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_PLUGIN_CLASSES')
    plugins_dic = eval(plugins_config)
    plugins = OrderedDict(plugins_dic)
    plugins_columns_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
    plugins_columns_dic = eval(plugins_columns_def)
    plugins_columns = OrderedDict(plugins_columns_dic)
    module_name = 'lisfloodutilities.gridding.lib.filters'
    try:
        for plugin in plugins:
            class_name = plugin
            class_args = plugins[plugin]
            module = importlib.import_module(module_name)
            class_instance = getattr(module, class_name)(plugins_columns, class_args, conf.var_code, quiet_mode)
            plugins_array.append(class_instance)
    except ImportError:
        print(f"Error: Could not import module '{module_name}'")
    except AttributeError:
        print(f"Error: Could not find class '{class_name}' in module '{module_name}'")
    return plugins_array


def get_timestamp_from_filename(conf: Config, filename: Path) -> datetime:
    file_timestamp = datetime.strptime(filename.name, conf.input_timestamp_pattern)
    if conf.force_time is not None:
        new_time = datetime.strptime(conf.force_time, "%H%M").time()
        file_timestamp = datetime.combine(file_timestamp.date(), new_time)
    return file_timestamp


def get_dataframes(conf: Config, kiwis_files: List[Path]) -> Tuple[List[pd.DataFrame], List[datetime], List[str], str]:
    filter_classes = get_filter_classes(conf)
    filepath_kiwis = []
    kiwis_timestamps = []
    kiwis_str_timestamps = []
    for filename_kiwis in kiwis_files:
        kiwis_timestamp = get_timestamp_from_filename(conf, filename_kiwis)
        kiwis_timestamp_str = kiwis_timestamp.strftime(FileUtils.DATE_PATTERN_CONDENSED_SHORT)
        filepath_kiwis.append(filename_kiwis)
        kiwis_str_timestamps.append(kiwis_timestamp_str)
        kiwis_timestamps.append(kiwis_timestamp)
    df_kiwis_array = []
    last_filter_class = KiwisFilter()
    for filter_class in filter_classes:
        last_filter_class = filter_class
        print_msg(f'{filter_class.get_class_description()}')
        df_kiwis_array = filter_class.filter(filepath_kiwis, kiwis_str_timestamps, df_kiwis_array)
    return df_kiwis_array, kiwis_timestamps, kiwis_str_timestamps, last_filter_class.COL_PROVIDER_ID

def get_providers_radius(conf: Config, kiwis_timestamp_24h: datetime) -> dict:
    """
    Get the configuration of providers and radius correspondent to the current daily date
    """
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
                               kiwis_dataframes: List[pd.DataFrame], filter_args: dict) -> List[pd.DataFrame]:
    """
    Applies the decumulation filter to the dataframes. The first dataframe in the array should be the Daily
    one and then the remaining 4 dataframes will be the 6 hourly ones.
    """

    # If there was no data interval defined for the current daily date, then no need to decumulate
    if not filter_args:
        return kiwis_dataframes
    
    plugins_columns_def = conf.get_config_field('PROPERTIES', 'KIWIS_FILTER_COLUMNS')
    plugins_columns_dic = eval(plugins_columns_def)
    filter_columns = OrderedDict(plugins_columns_dic)
    filter_class = DowgradedDailyTo6HourlyObservationsKiwisFilter(filter_columns, filter_args, conf.var_code, quiet_mode)
    df_kiwis_array = filter_class.filter(kiwis_filepaths, kiwis_str_timestamps, kiwis_dataframes)
    return df_kiwis_array

def processable_file(file_timestamp: datetime, start_date: datetime = None, end_date: datetime = None) -> bool:
    return (start_date is not None and start_date <= file_timestamp and
            end_date is not None and file_timestamp <= end_date)

def get_timestamp_from_filename(conf: Config, filename: Path) -> datetime:
    file_timestamp = datetime.strptime(filename.name, conf.input_timestamp_pattern)
    if conf.force_time is not None:
        new_time = datetime.strptime(conf.force_time, "%H%M").time()
        file_timestamp = datetime.combine(file_timestamp.date(), new_time)
    return file_timestamp

def get_24h_kiwis_paths(conf: Config, infolder: Path):
    kiwis_paths = []
    for filename_kiwis in sorted(infolder.rglob(conf.input_wildcard)):
        kiwis_timestamp = get_timestamp_from_filename(conf, filename_kiwis)
        if processable_file(kiwis_timestamp, conf.start_date, conf.end_date):
            kiwis_paths.append((filename_kiwis, kiwis_timestamp))
    return kiwis_paths

def print_statistics(provider_ids: List[str], df_kiwis_24h: pd.DataFrame, df_kiwis_array_6h_before: List[pd.DataFrame],
                     df_kiwis_array_6h_after: List[pd.DataFrame], column_provider_id_24h: str, column_provider_id_6h: str,
                     kiwis_timestamps_6h: List[datetime]):
    # Count the number of rows by provider_id daily data
    df_24h_count = df_kiwis_24h.groupby(column_provider_id_24h).size().reset_index(name='count')
    df_24h_count.rename(columns={column_provider_id_24h: 'provider_id'}, inplace=True)

    i = 0
    for df_kiwis_6h_before in df_kiwis_array_6h_before:
        df_kiwis_6h_after = df_kiwis_array_6h_after[i]
        # Count the number of rows by provider_id 6hourly data before decumulation
        df_6h_before_count = df_kiwis_6h_before.groupby(column_provider_id_6h).size().reset_index(name='count')
        df_6h_before_count.rename(columns={column_provider_id_6h: 'provider_id'}, inplace=True)
        filtered_before_df = df_6h_before_count[df_6h_before_count['provider_id'].isin(provider_ids)]
        # Count the number of rows by provider_id 6hourly data
        df_6h_after_count = df_kiwis_6h_after.groupby(column_provider_id_6h).size().reset_index(name='count')
        df_6h_after_count.rename(columns={column_provider_id_6h: 'provider_id'}, inplace=True)
        # Merge both dataframes by provider_id
        merged_df = df_24h_count.merge(df_6h_after_count, on='provider_id', how='outer')
        merged_df.rename(columns={'count_x': 'count_24h', 'count_y': 'count_6h_after'}, inplace=True)
        # Fill NaN values with 0 for the non matching count column
        merged_df['count_24h'] = merged_df['count_24h'].fillna(0)
        merged_df['count_6h_after'] = merged_df['count_6h_after'].fillna(0)
        # Filter rows that are in the decumulated provider IDs
        filtered_df = merged_df[merged_df['provider_id'].isin(provider_ids)]
        TOTAL_6H_STATIONS = len(df_kiwis_6h_before)
        TIMESTAMP = kiwis_timestamps_6h[i].strftime(TIMESTAMP_PATTERN)
        for index, row in filtered_df.iterrows():
            PROVIDER_ID = row['provider_id']
            TOTAL_24H_STATIONS = row['count_24h']
            count_before = 0
            count_before_series = filtered_before_df.loc[filtered_before_df['provider_id'] == PROVIDER_ID, 'count']
            if not count_before_series.empty:
                count_before = count_before_series.iloc[0]
            count_after = row['count_6h_after']
            DECUMULATED_24H_STATIONS = count_after - count_before
            DECUMULATED_STATIONS_24H_PERCENT = 0.0
            if TOTAL_24H_STATIONS > 0.0:
                DECUMULATED_STATIONS_24H_PERCENT = 100.0 * DECUMULATED_24H_STATIONS / TOTAL_24H_STATIONS
            DECUMULATED_STATIONS_RELATIVE_TO_6H_PERCENT = 0.0
            if TOTAL_6H_STATIONS > 0.0:
                DECUMULATED_STATIONS_RELATIVE_TO_6H_PERCENT = 100.0 * DECUMULATED_24H_STATIONS / TOTAL_6H_STATIONS
            stats_string = (
                f'#APP_STATS: {{"TIMESTAMP": "{TIMESTAMP}", "PROVIDER_ID": {PROVIDER_ID}, '
                f'"TOTAL_24H_STATIONS": {TOTAL_24H_STATIONS}, "DECUMULATED_24H_STATIONS": {DECUMULATED_24H_STATIONS}, '
                f'"DECUMULATED_STATIONS_24H_PERCENT": {DECUMULATED_STATIONS_24H_PERCENT:.2f}, "TOTAL_6H_STATIONS": {TOTAL_6H_STATIONS}, '
                f'"DECUMULATED_STATIONS_RELATIVE_TO_6H_PERCENT": {DECUMULATED_STATIONS_RELATIVE_TO_6H_PERCENT:.2f}}}'
            )
            print_msg(stats_string)
        i += 1

def run(conf_24h: Config, conf_6h: Config, beginning_of_interval_offset: int, kiwis_24h_06am_path: Path, kiwis_6h_12pm_path: Path,
        kiwis_6h_18pm_path: Path, kiwis_6h_12am_path: Path, kiwis_6h_06am_path: Path, input_path_6h: Path, output_path: Path = None):
    """
    While processing the 4 grids of 6hourly precipitation Day1 12:00, 18:00 and Day2 00:00, 06:00 we will use the daily precipitation
    of the corresponding period meaning Day2 06:00 to decumulate its values where there is missing 6hourly precipitation.
    For each daily precipitation observation from one of the providers above, if there is no 6h observation in a radius of 1.5 km,
    insert the daily value divided by 4 in all 4 6hourly datasets.
    If there is one or more 6hourly stations in the radius and the station have less then 4 values, insert the missing values in the
    corresponding 6hourly dataset by and changing its value using the formula (PR - Sum(PR6)) / (number of missing values),
    if and only if the resulting value is positive (>=0).
    Daily stations that have a real 6h station in the radius (1.5km) that is complete, meaning it has all four 6hourly values. The decumulation is NOT done.

    - read 1 daily KiWIS-file and 4 corresponding 6-hourly KiWIS-files
    - calculate decumulated values
    - write 4 new 6-hourly KiWIS-files including decumulated values
    """

    global quiet_mode

    print_msg('Start reading files')
    df_kiwis_array_24h, kiwis_timestamps_24h, kiwis_str_timestamps_24h, column_provider_id_24h = get_dataframes(conf_24h, [kiwis_24h_06am_path])
    df_kiwis_array_6h, kiwis_timestamps_6h, kiwis_str_timestamps_6h, column_provider_id_6h = get_dataframes(conf_6h, [kiwis_6h_12pm_path,
                                                                                                                      kiwis_6h_18pm_path,
                                                                                                                      kiwis_6h_12am_path,
                                                                                                                      kiwis_6h_06am_path])
    # Check timestamps are correct
    if not ((kiwis_timestamps_24h[0] + timedelta(hours=beginning_of_interval_offset)) == kiwis_timestamps_6h[3] and
            (kiwis_timestamps_24h[0] + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=6)) == kiwis_timestamps_6h[2] and
            (kiwis_timestamps_24h[0] + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=12)) == kiwis_timestamps_6h[1] and
            (kiwis_timestamps_24h[0] + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=18)) == kiwis_timestamps_6h[0]):
        raise ArgumentTypeError("The input kiwis do not respect the expected timestamps.")

    kiwis_dataframes = df_kiwis_array_24h
    kiwis_dataframes.extend(df_kiwis_array_6h)
    kiwis_timestamps = kiwis_str_timestamps_24h
    kiwis_timestamps.extend(kiwis_str_timestamps_6h)
    kiwis_filepaths = [kiwis_24h_06am_path, kiwis_6h_12pm_path, kiwis_6h_18pm_path, kiwis_6h_12am_path, kiwis_6h_06am_path]

    filter_args = get_providers_radius(conf_24h, kiwis_timestamps_24h[0])
    df_kiwis_array = get_decumulated_dataframes(conf_24h, kiwis_filepaths, kiwis_timestamps, kiwis_dataframes, filter_args)

    provider_ids = list(filter_args.keys())
    print_statistics(provider_ids, df_kiwis_array_24h[0], df_kiwis_array_6h, df_kiwis_array[1:],
                     column_provider_id_24h, column_provider_id_6h, kiwis_timestamps_6h)

    # Write output files
    i = 0
    for kiwis_filepath in kiwis_filepaths[1:]:
        i += 1
        if output_path is not None:
            outfile = str(kiwis_filepath).replace(str(input_path_6h), str(output_path))
            filepath = Path(outfile)
            # Create the output parent folders if not exist yet
            Path(filepath.parent).mkdir(parents=True, exist_ok=True)
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

def get_6hourly_filepath(parser: ArgumentParser, conf: Config, kiwis_6h_folder_path: Path, kiwis_6h_timestamp: datetime) -> Path:
    kiwis_6h_filename = kiwis_6h_timestamp.strftime(conf.input_timestamp_pattern)
    for filename_kiwis in sorted(kiwis_6h_folder_path.rglob(kiwis_6h_filename)):
        return filename_kiwis
    parser.error(f'Input file {kiwis_6h_filename} does not exist in folder {kiwis_6h_folder_path}.')
    return None

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
                            output_folder=None,
                            use_beginning_of_interval=False)

        parser.add_argument("-d", "--pr24h", dest="kiwis_24h_folder_path", required=True, type=FileUtils.folder_type,
                            help="Set the input kiwis file folder containing daily precipitation.",
                            metavar="/path/to/pr")
        parser.add_argument("-g", "--pr6h", dest="kiwis_6h_folder_path", required=True, type=FileUtils.folder_type,
                            help="Set the input kiwis file folder containing 6 hourly precipitation.",
                            metavar="/path/to/pr6")
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
        parser.add_argument("-s", "--start", dest="start_date",
                            help="Set the start date and time from which data is imported [default: date defining the time units inside the config file]",
                            metavar="YYYYMMDDHHMISS")
        parser.add_argument("-e", "--end", dest="end_date",
                            help="Set the end date and time until which data is imported [default: %(default)s]",
                            metavar="YYYYMMDDHHMISS")
        parser.add_argument("-b", "--boi", dest="use_beginning_of_interval", action="store_true",
                            help="Indicate that the daily timesteps are at the beginning of the interval [default: %(default)s]")
        parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Set script output into quiet mode [default: %(default)s]")

        # process options
        args = parser.parse_args(argv)
        quiet_mode = args.quiet

        configuration_base_folder = os.path.join(program_path, '../src/lisfloodutilities/gridding/configuration')

        start_date = None
        try:
            start_date = datetime.strptime(args.start_date, FileUtils.DATE_PATTERN_CONDENSED)
            start_date_str =  start_date.strftime(FileUtils.DATE_PATTERN_SEPARATED)
            print_msg(f"Start Date:  {start_date_str}")
        except Exception as t:
            print_msg(f"Start Date: DEFAULT")

        end_date = None
        try:
            end_date = datetime.strptime(args.end_date, FileUtils.DATE_PATTERN_CONDENSED)
        except Exception as t:
            print_msg(f"Warning: {args.end_date} does not seem a valid date. Setting up default end date")
            end_date = datetime.strptime(END_DATE_DEFAULT, FileUtils.DATE_PATTERN_CONDENSED)
        end_date_str =  end_date.strftime(FileUtils.DATE_PATTERN_SEPARATED)
        print_msg(f"End Date:  {end_date_str}")

        if args.config_base_path is not None and len(args.config_base_path) > 0:
            configuration_base_folder = args.config_base_path

        file_utils_24h = FileUtils(args.variable_code_24h, quiet_mode)
        config_type_path_24h = file_utils_24h.get_config_type_path(configuration_base_folder, args.config_type)
        config_filename_24h = file_utils_24h.get_config_file(config_type_path_24h)

        file_utils_6h = FileUtils(args.variable_code_6h, quiet_mode)
        config_type_path_6h = file_utils_6h.get_config_type_path(configuration_base_folder, args.config_type)
        config_filename_6h = file_utils_6h.get_config_file(config_type_path_6h)

        interpolation_mode = 'adw'
        memory_save_mode = 5
        conf_24h = Config(config_filename_24h, start_date, end_date, quiet_mode, interpolation_mode, memory_save_mode)
        conf_6h = Config(config_filename_6h, start_date, end_date, quiet_mode, interpolation_mode, memory_save_mode)
        print_msg(f"Config File 24h: {config_filename_24h}")
        print_msg(f"Config File 6h: {config_filename_6h}")

        if args.output_folder is not None and len(args.output_folder) > 0:
            output_path = Path(args.output_folder)
            print_msg(f"Output Folder: {output_path}")
        else:
            output_path = None
            print_msg(f"Output Folder: 6hourly kiwis files will be overwritten")

        print_msg(f"Timesteps are beginning of the interval: {args.use_beginning_of_interval}")

        kiwis_24h_paths = get_24h_kiwis_paths(conf_24h, Path(args.kiwis_24h_folder_path))
        kiwis_6h_folder_path = Path(args.kiwis_6h_folder_path)
        
        # Defines the offset for the 6hourly data when the daily data is at
        # the beginning of the interval while the 6hourly files should be at the end of the interval
        beginning_of_interval_offset = 0
        if args.use_beginning_of_interval:
            beginning_of_interval_offset = 24

        for filename_kiwis, kiwis_timestamp in kiwis_24h_paths:
            kiwis_24h_06am_path = get_existing_file_path(parser, str(filename_kiwis))
            print_msg(f"Daily {args.variable_code_24h} kiwis file:  {kiwis_24h_06am_path}")
            kiwis_6h_06am_timestamp = kiwis_timestamp + timedelta(hours=beginning_of_interval_offset)
            kiwis_6h_12am_timestamp = kiwis_timestamp + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=6)
            kiwis_6h_18pm_timestamp = kiwis_timestamp + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=12)
            kiwis_6h_12pm_timestamp = kiwis_timestamp + timedelta(hours=beginning_of_interval_offset) - timedelta(hours=18)

            kiwis_6h_06am_path = get_6hourly_filepath(parser, conf_6h, kiwis_6h_folder_path, kiwis_6h_06am_timestamp)
            kiwis_6h_12am_path = get_6hourly_filepath(parser, conf_6h, kiwis_6h_folder_path, kiwis_6h_12am_timestamp)
            kiwis_6h_18pm_path = get_6hourly_filepath(parser, conf_6h, kiwis_6h_folder_path, kiwis_6h_18pm_timestamp)
            kiwis_6h_12pm_path = get_6hourly_filepath(parser, conf_6h, kiwis_6h_folder_path, kiwis_6h_12pm_timestamp)

            print_msg(f"6hourly {args.variable_code_6h} kiwis file 06:00:  {kiwis_6h_06am_path}")
            print_msg(f"6hourly {args.variable_code_6h} kiwis file 00:00:  {kiwis_6h_12am_path}")
            print_msg(f"6hourly {args.variable_code_6h} kiwis file 18:00:  {kiwis_6h_18pm_path}")
            print_msg(f"6hourly {args.variable_code_6h} kiwis file 12:00:  {kiwis_6h_12pm_path}")

            run(conf_24h, conf_6h, beginning_of_interval_offset, kiwis_24h_06am_path, kiwis_6h_12pm_path, kiwis_6h_18pm_path,
                kiwis_6h_12am_path, kiwis_6h_06am_path, input_path_6h=kiwis_6h_folder_path, output_path=output_path)
        return 0
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()

