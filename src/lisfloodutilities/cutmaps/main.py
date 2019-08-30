"""
Copyright 2019 European Union

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
import sys

from .cutlib import get_filelist, get_cuts, cutmap, logger


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

    def add_args(self):
        group_mask = self.add_mutually_exclusive_group()
        group_filelist = self.add_mutually_exclusive_group()

        group_mask.add_argument("-m", "--mask", help='mask file cookie-cutter, .map if pcraster, .nc if netcdf')
        group_mask.add_argument("-c", "--cuts", help='Cut coordinates in the form Xmin-Xmax:Ymin-Ymax')
        group_filelist.add_argument("-l", "--list", help='list of files to be cut, in a text file, one file per line, '
                                                         'main variable i.e. pr/e0/tx/tavg must be '
                                                         'last variable in the input file')
        group_filelist.add_argument("-f", "--folder", help='Directory with netCDF files to be cut')

        self.add_argument("-o", "--outpath", help='path where to save cut files', default='./cutmaps_out', required=True)


def main(cliargs):
    parser = ParserHelpOnError(description='Convert PCRaster maps to a NetCDF map stack')

    parser.add_args()

    args = parser.parse_args(cliargs)

    mask = args.mask
    cuts = args.cuts
    filelist = args.list
    pathout = args.outpath
    input_folder = args.folder
    logger.info('Cutting using: %s\n Files to cut from: %s\n Output: %s\n ', mask or cuts, filelist or input_folder, pathout)

    list_to_cut = get_filelist(filelist, input_folder)
    x_max, x_min, y_max, y_min = get_cuts(cuts, mask)

    # walk through list_to_cut
    for file_to_cut in list_to_cut:
        filename, ext = os.path.splitext(file_to_cut)
        fileout = os.path.join(pathout, os.path.basename(file_to_cut))
        if ext != '.nc':
            logger.warning('%s is not in netcdf format, skipping...', file_to_cut)
            continue
        if os.path.isfile(fileout) and os.path.exists(fileout):
            logger.warning('%s already existing. This file will not be overwritten', fileout)
            continue

        cutmap(file_to_cut, fileout, x_max, x_min, y_max, y_min)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
