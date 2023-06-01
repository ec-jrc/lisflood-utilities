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
from lisfloodutilities.gridding.lib.utils import Printable, Dem, Config, FileUtils, GriddingUtils
from lisfloodutilities.gridding.lib.writers import NetCDFWriter, GDALWriter


END_DATE_DEFAULT = datetime.now().strftime('%Y%m%d060000')

quiet_mode = False


def print_msg(msg: str = ''):
    global quiet_mode
    if not quiet_mode:
        print(msg)

def run(config_filename: str, infolder: str, outfolder_or_file: str, processing_dates_file: str, file_utils: FileUtils,
        output_tiff: bool, overwrite_output: bool, start_date: datetime = None, end_date: datetime = None):
    """
    Interpolate text files containing (x, y, value) using inverse distance interpolation.
    Produces as output, either a netCDF file containing all the grids or one TIFF file per grid.
    1. Read the config files
    2. Get the ordered list of files
    3. Interpolate the grids performing height correction on some variables
    4. Write the resulting grids
    """
    global quiet_mode

    conf = Config(config_filename, start_date, end_date, quiet_mode)

    if conf.start_date > conf.end_date:
        raise ArgumentTypeError("Start date is greater than End date.")

    dates_to_process = file_utils.read_processing_dates_file(processing_dates_file)
    grid_utils = GriddingUtils(conf, quiet_mode)

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

    if output_tiff:
        outfolder = outfolder_or_file
        output_writer = GDALWriter(conf, overwrite_output, quiet_mode)
        for filename in sorted(Path(infolder).rglob(inwildcard)):
            file_timestamp = file_utils.get_timestamp_from_filename(filename) + timedelta(days=netcdf_offset_file_date)
            if not file_utils.processable_file(file_timestamp, dates_to_process, conf.start_date, conf.end_date):
                continue  # Skip processing file
            print_msg(f'Processing file: {filename}')
            outfile = str(filename).replace(infolder, outfolder)
            outfilepath = Path(outfile).with_suffix('.tiff')
            # Create the output parent folders if not exist yet
            Path(outfilepath.parent).mkdir(parents=True, exist_ok=True)
            output_writer.open(Path(outfilepath))
            grid_data = grid_utils.generate_grid(filename)
            output_writer.write(grid_data, file_timestamp)
            output_writer.close()
    else: # NetCDF
        outfile = outfolder_or_file
        output_writer = NetCDFWriter(conf, overwrite_output, quiet_mode)
        output_writer.open(Path(outfile))
        for filename in sorted(Path(infolder).rglob(inwildcard)):
            file_timestamp = file_utils.get_timestamp_from_filename(filename) + timedelta(days=netcdf_offset_file_date)
            if not file_utils.processable_file(file_timestamp, dates_to_process, conf.start_date, conf.end_date):
                continue  # Skip processing file
            print_msg(f'Processing file: {filename}')
            grid_data = grid_utils.generate_grid(filename)
            output_writer.write(grid_data, file_timestamp)
        output_writer.close()
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
    This script interpolates meteo input variables data into either a single NETCDF4 file or one GEOTIFF file per timestep.
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
                            overwrite_output=False,
                            start_date='',
                            end_date=END_DATE_DEFAULT)

        parser.add_argument("-i", "--in", dest="infolder", required=True, type=FileUtils.folder_type,
                            help="Set input folder path with kiwis/point files",
                            metavar="input_folder")
        parser.add_argument("-o", "--out", dest="outfolder_or_file", required=True, type=FileUtils.file_or_folder,
                            help="Set output folder base path for the tiff files or the netCDF file path.",
                            metavar="{output_folder, netcdf_file}")
        parser.add_argument("-c", "--conf", dest="config_type", required=True,
                            help="Set the grid configuration type to use.",
                            metavar="{5x5km, 1arcmin,...}")
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
                            help="Outputs a tiff file per timestep instead of the default single netCDF [default: %(default)s]")
        parser.add_argument("-f", "--force", dest="overwrite_output", action="store_true",
                            help="Force write to existing file. TIFF files will be overwritten and netCDF file will be appended. [default: %(default)s]")

        # process options
        args = parser.parse_args(argv)
        quiet_mode = args.quiet

        configuration_base_folder = os.path.join(program_path, '../src/lisfloodutilities/gridding/configuration')

        file_utils = FileUtils(args.variable_code, quiet_mode)

        config_type_path = file_utils.get_config_type_path(configuration_base_folder, args.config_type)

        config_filename = file_utils.get_config_file(config_type_path)

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

        print_msg(f"Input Folder:  {args.infolder}")
        print_msg(f"Overwrite Output: {args.overwrite_output}")
        if args.out_tiff:
            print_msg("Output Type: TIFF")
            print_msg(f"Output Folder: {args.outfolder_or_file}")
        else:
            print_msg("Output Type: netCDF")
            print_msg(f"Output File: {args.outfolder_or_file}")
        print_msg(f"Processing Dates File: {args.processing_dates_file}")
        print_msg(f"Config File: {config_filename}")

        run(config_filename, args.infolder, args.outfolder_or_file, args.processing_dates_file, 
            file_utils, args.out_tiff, args.overwrite_output, start_date, end_date)
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()

