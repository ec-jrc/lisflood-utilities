import argparse
import logging
import os
import sys

from lisflood.utilities.cutmaps.cutlib import get_filelist, get_cuts, cutmap

logger = logging.getLogger()


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
    logger.info('{} {} {} {}'.format(x_min, x_max, y_min, y_max))

    # walk through list_to_cut
    for file_to_cut in list_to_cut:
        filename, ext = os.path.splitext(file_to_cut)
        fileout = os.path.join(pathout, os.path.basename(file_to_cut))
        if ext != '.nc':
            logger.warning('{} is not in netcdf format, skipping...'.format(file_to_cut))
            continue
        if os.path.isfile(fileout) and os.path.exists(fileout):
            logger.warning(fileout, 'Alreay existing, be sure you are not swiping original file with cut version.'
                           'This file will not be overwritten')
            continue

        cutmap(file_to_cut, fileout, x_max, x_min, y_max, y_min)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
