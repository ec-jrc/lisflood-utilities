__author__="Goncalo Gomes"
__date__="$Mar 26, 2023 12:01:00$"
__version__="0.1"
__updated__="$Mar 28, 2023 16:01:00$"

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
from datetime import datetime, timedelta
from lisfloodutilities.gridding.lib.utils import Printable, Dem, Config, FileUtils, GriddingUtils, KiwisLoader
from lisfloodutilities.gridding.lib.writers import NetCDFWriter, GDALWriter


END_DATE_DEFAULT = datetime.now().strftime('%Y%m%d060000')

quiet_mode = False


def print_msg(msg: str = ''):
    global quiet_mode
    if not quiet_mode:
        print(msg)

def interpolation_mode_type(mode: str) -> str:
    if not mode or mode not in Config.INTERPOLATION_MODES:
        raise ArgumentTypeError(f'You must select a mode out of {list(Config.INTERPOLATION_MODES.keys())}.')
    return mode

def memory_save_mode_type(mode: str) -> str:
    if mode not in Config.MEMORY_SAVE_MODES:
        raise ArgumentTypeError(f'You must select a RAM save mode out of {list(Config.MEMORY_SAVE_MODES.keys())}.')
    return mode

def run(config_filename: str, infolder: str, output_file: str, processing_dates_file: str, file_utils: FileUtils,
        output_tiff: bool, output_netcdf: bool, overwrite_output: bool, start_date: datetime = None, end_date: datetime = None,
        interpolation_mode: str = 'adw', use_broadcasting: bool = False, memory_save_mode: str = None):
    """
    Interpolate text files containing (x, y, value) using inverse distance interpolation.
    Produces as output, either a netCDF file containing all the grids or one TIFF file per grid.
    1. Read the config files
    2. Get the ordered list of files
    3. Interpolate the grids performing height correction on some variables
    4. Write the resulting grids
    """
    global quiet_mode

    conf = Config(config_filename, start_date, end_date, quiet_mode, interpolation_mode, memory_save_mode)

    if conf.start_date > conf.end_date:
        raise ArgumentTypeError("Start date is greater than End date.")

    dates_to_process = file_utils.read_processing_dates_file(processing_dates_file)
    grid_utils = GriddingUtils(conf, quiet_mode, use_broadcasting)

    size_lons = len(conf.dem_lons)
    size_lats = len(conf.dem_lats)
    print_msg("cellSizeX: %s cellSizeY: %s dem_lons: %s dem_lats: %s minX: %s maxY: %s" % (
        str(conf.dem_cell_size_x),
        str(conf.dem_cell_size_y),
        str(size_lons), str(size_lats),
        str(conf.dem_min_x),
        str(conf.dem_max_y)))

    print_msg('Start reading files')
    inwildcard = conf.var_code + FileUtils.FILES_WILDCARD

    netcdf_offset_file_date = int(conf.get_config_field('VAR_TIME','OFFSET_FILE_DATE'))

    outfile = output_file
    if output_tiff:
        output_writer_tiff = GDALWriter(conf, overwrite_output, quiet_mode)
    if output_netcdf:
        output_writer_netcdf = NetCDFWriter(conf, overwrite_output, quiet_mode)
        output_writer_netcdf.open(Path(outfile))
    file_loader = KiwisLoader(conf, Path(infolder), dates_to_process, overwrite_output, quiet_mode)
    for filename in file_loader:
        file_timestamp = file_utils.get_timestamp_from_filename(filename) + timedelta(days=netcdf_offset_file_date)
        print_msg(f'Processing file: {filename}')
        if output_tiff:
            outfilepath = filename.with_suffix('.tiff')
            output_writer_tiff.open(outfilepath)
        grid_data = grid_utils.generate_grid(filename)
        if output_netcdf:
            output_writer_netcdf.write(grid_data, file_timestamp)
        if output_tiff:
            output_writer_tiff.write(grid_data, file_timestamp)
            output_writer_tiff.close()
    if output_netcdf:
        output_writer_netcdf.close()
    print_msg('Finished writing files')


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script interpolates meteo input variables data into a single NETCDF4 file and, if selected, generates also a GEOTIFF file per timestep.
    The resulting netCDF is CF-1.6 compliant.
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
                            out_tiff=False,
                            out_netcdf=False,
                            overwrite_output=False,
                            start_date='',
                            end_date=END_DATE_DEFAULT,
                            interpolation_mode='adw',
                            use_broadcasting=False,
                            memory_save_mode='0')

        parser.add_argument("-i", "--in", dest="in_file_or_folder", required=True, type=FileUtils.file_or_folder,
                            help="Set a single input kiwis file or folder path containing all the kiwis files.",
                            metavar="/path/to/pr200102150600_all.kiwis")
        parser.add_argument("-o", "--out", dest="output_file", required=False, type=FileUtils.file_type,
                            help="Set the output netCDF file path containing all the timesteps between start and end dates.",
                            metavar="/path/to/pr2001.nc")
        parser.add_argument("-c", "--conf", dest="config_type", required=True,
                            help="Set the grid configuration type to use.",
                            metavar="{5x5km, 1arcmin,...}")
        parser.add_argument("-p", "--pathconf", dest="config_base_path", required=False, type=FileUtils.folder_type,
                            help="Overrides the base path where the configurations are stored.",
                            metavar="/path/to/config")
        parser.add_argument("-v", "--var", dest="variable_code", required=True,
                            help="Set the variable to be processed.",
                            metavar="{pr,pd,tn,tx,ws,rg,...}")
        parser.add_argument("-d", "--dates", dest="processing_dates_file", required=False, type=FileUtils.file_type,
                            help="Set file containing a list of filenames to be processed in the form of <var>YYYYMMDDHHMI_YYYYMMDDHHMISS.txt",
                            metavar="files2process.txt")
        parser.add_argument("-s", "--start", dest="start_date",
                            help="Set the start date and time from which data is imported [default: date defining the time units inside the config file]",
                            metavar="YYYYMMDDHHMISS")
        parser.add_argument("-e", "--end", dest="end_date",
                            help="Set the end date and time until which data is imported [default: %(default)s]",
                            metavar="YYYYMMDDHHMISS")
        parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Set script output into quiet mode [default: %(default)s]")
        parser.add_argument("-t", "--tiff", dest="out_tiff", action="store_true",
                            help="Outputs a tiff file per timestep [default: %(default)s]")
        parser.add_argument("-n", "--netcdf", dest="out_netcdf", action="store_true",
                            help="Outputs a single netCDF with all the timesteps [default: %(default)s]")
        parser.add_argument("-f", "--force", dest="overwrite_output", action="store_true",
                            help="Force write to existing file. TIFF files will be overwritten and netCDF file will be appended. [default: %(default)s]")
        parser.add_argument("-m", "--mode", dest="interpolation_mode", required=False, type=interpolation_mode_type,
                            help="Set interpolation mode. [default: %(default)s]",
                            metavar=f"{list(Config.INTERPOLATION_MODES.keys())}")
        parser.add_argument("-r", "--ramsavemode", dest="memory_save_mode", required=False, type=memory_save_mode_type,
                            help="Set memory save mode level. Used to reduce memory usage [default: %(default)s]",
                            metavar=f"{list(Config.MEMORY_SAVE_MODES.keys())}")
        parser.add_argument("-b", "--broadcast", dest="use_broadcasting", action="store_true",
                            help="When set, computations will run faster in full broadcasting mode but require more memory. [default: %(default)s]")

        # process options
        args = parser.parse_args(argv)
        quiet_mode = args.quiet

        configuration_base_folder = os.path.join(program_path, '../src/lisfloodutilities/gridding/configuration')
        if args.config_base_path is not None and len(args.config_base_path) > 0:
            configuration_base_folder = args.config_base_path

        file_utils = FileUtils(args.variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, args.config_type)

        config_filename = file_utils.get_config_file(config_type_path)

        if not args.out_tiff and not args.out_netcdf:
            parser.error(f'You must choose at least one output format, TIFF (--tiff) and/or netCDF (--netcdf)')
        if args.out_netcdf:
            if args.output_file is None:
                parser.error("--netcdf requires defining the output file with --out.")
            else:
                print_msg("Output Type: netCDF")
                print_msg(f"Output File: {args.output_file}")
        # TIFF output is written to the folder where each of the kiwis files exist
        if args.out_tiff:
            print_msg("Output Type: TIFF")
            output_path = Path(args.in_file_or_folder)
            if output_path.is_file():
                output_path = output_path.parent
            print_msg(f"Output Folder: {output_path}")

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

        print_msg(f"Input Folder:  {args.in_file_or_folder}")
        print_msg(f"Overwrite Output: {args.overwrite_output}")
        print_msg(f"Interpolation Mode: {args.interpolation_mode}")
        print_msg(f"RAM Save Mode: {args.memory_save_mode}")
        print_msg(f"Broadcasting: {args.use_broadcasting}")
        print_msg(f"Processing Dates File: {args.processing_dates_file}")
        print_msg(f"Config File: {config_filename}")

        run(config_filename, args.in_file_or_folder, args.output_file, args.processing_dates_file, file_utils, args.out_tiff,
            args.out_netcdf, args.overwrite_output, start_date, end_date, args.interpolation_mode, args.use_broadcasting,
            args.memory_save_mode)
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()

