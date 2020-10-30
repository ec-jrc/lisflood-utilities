"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

cutmaps
A tool to cut netcdf files
"""

import argparse
import os
import shutil
import sys

from .. import version, logger
from .cutlib import mask_from_ldd, get_filelist, get_cuts, cutmap
from ..nc2pcr import convert


def parse_and_check_args(parser, cliargs):
    args = parser.parse_args(cliargs)
    if args.mask and args.cuts:
        parser.error('[--mask | --cuts] arguments are mutually exclusive')
    if not (args.mask or args.cuts) and not (args.ldd and args.stations):
        parser.error('(--mask | --cuts | [--ldd, --stations]) You need to pass mask path or cuts coordinates '
                     'or a list of stations along with LDD path')
    if (args.mask or args.cuts) and (args.ldd or args.stations):
        parser.error('(--mask | --cuts | [--ldd, --stations]) '
                     '--mask, --cuts and --ldd and --stations arguments are mutually exclusive')
    return args


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):
        group_mask = self.add_argument_group(title='Cut with a provided mask or a bounding box or '
                                                   'create mask cookie-cutter on-fly from stations list and ldd map')
        group_filelist = self.add_mutually_exclusive_group(required=True)
        group_mask.add_argument("-m", "--mask", help='mask file cookie-cutter, .map if pcraster, .nc if netcdf')
        group_mask.add_argument("-c", "--cuts", help='Cut coordinates in the form lonmin_lonmax:latmin_latmax')
        group_mask.add_argument("-l", "--ldd", help='Path to LDD file')
        group_mask.add_argument("-N", "--stations",
                                help='Path to stations.txt file.'
                                     'Read documentation to know about the format')
        group_mask.add_argument("-C", "--clonemap",
                                help='Path to PCRaster clonemap; used to convert ldd.nc to ldd.map')

        group_filelist.add_argument("-f", "--folder", help='Directory with netCDF files to be cut')
        group_filelist.add_argument("-S", "--static-data",
                                    help='Directory with EFAS/GloFAS static maps. '
                                         'Output files will have same directories structure')

        self.add_argument("-o", "--outpath", help='path where to save cut files',
                          default='./cutmaps_out', required=True)
        self.add_argument("-W", "--overwrite", help='Set flag to overwrite existing files',
                          default=False, required=False, action='store_true')


def main(cliargs):
    parser = ParserHelpOnError(description='Cut netCDF file: {}'.format(version))
    parser.add_args()
    args = parse_and_check_args(parser, cliargs)
    mask = args.mask
    cuts = args.cuts

    ldd = args.ldd
    stations = args.stations

    input_folder = args.folder
    static_data_folder = args.static_data
    overwrite = args.overwrite
    pathout = args.outpath
    if not os.path.exists(pathout):
        logger.warning('\nOutput folder %s not existing. Creating it...', pathout)
        os.mkdir(pathout)
    if ldd and stations:
        logger.info('\nTry to produce a mask from LDD and stations points: %s %s', ldd, stations)
        if ldd.endswith('.nc'):
            # convert ldd.nc to a pcraster map as we are going to use pcraster commands
            clonemap = args.clonemap
            ldd = convert(ldd, clonemap, '.map')
        mask, outlets_nc = mask_from_ldd(ldd, stations)
        # copy outlets.nc (produced from stations txt file) and the new mask to output folder
        shutil.copy(outlets_nc, os.path.join(pathout, 'my_outlets.nc'))
        shutil.copy(mask, os.path.join(pathout, 'my_mask.map'))

    logger.info('\n\nCutting using: %s\n Files to cut from: %s\n Output: %s\n Overwrite existing: %s\n\n',
                mask or cuts,
                input_folder or static_data_folder,
                pathout, overwrite)

    list_to_cut = get_filelist(input_folder, static_data_folder)
    x_min, x_max, y_min, y_max = get_cuts(cuts=cuts, mask=mask)

    # walk through list_to_cut
    for file_to_cut in list_to_cut:

        filename, ext = os.path.splitext(file_to_cut)

        # localdir used only with static_data_folder.
        # It will track folder structures in a EFAS/GloFAS like setup and replicate it in output folder
        localdir = os.path.dirname(file_to_cut)\
            .replace(os.path.dirname(static_data_folder), '')\
            .lstrip('/') if static_data_folder else ''

        fileout = os.path.join(pathout, localdir, os.path.basename(file_to_cut))
        if os.path.isdir(file_to_cut) and static_data_folder:
            # just create folder
            os.makedirs(fileout, exist_ok=True)
            continue
        if ext != '.nc':
            if static_data_folder:
                logger.warning('%s is not in netcdf format, just copying to ouput folder', file_to_cut)
                shutil.copy(file_to_cut, fileout)
            else:
                logger.warning('%s is not in netcdf format, skipping...', file_to_cut)
            continue
        elif os.path.isfile(fileout) and os.path.exists(fileout) and not overwrite:
            logger.warning('%s already existing. This file will not be overwritten', fileout)
            continue

        cutmap(file_to_cut, fileout, x_min, x_max, y_min, y_max)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
