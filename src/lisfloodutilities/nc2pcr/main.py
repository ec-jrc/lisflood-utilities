"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

nc2prc
A tool that convert netCDF files to PCRaster maps
"""

import sys
import argparse

from . import convert
from .. import version


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):
        self.add_argument('-i', '--input', required=True,
                          help='Path to input NetCDF. Current version does not work with mapstacks',
                          metavar='input')
        self.add_argument('-o', '--output', help='Path to PCRaster output map.', default='./out.map',
                          metavar='output')
        self.add_argument('-c', '--clonemap', help='Path to PCRaster clone map.',
                          metavar='clonemap')
        self.add_argument("-l", "--is_ldd", help='Set flag if input file is a LDD map', default=False,
                          required=False, action='store_true')


def main_script():
    sys.exit(main(sys.argv[1:]))


def main(args):
    parser = ParserHelpOnError(description='Convert a netCDF file to a PCRaster map {}'.format(version))

    parser.add_args()
    parsed_args = parser.parse_args(args)
    configuration = {'inf': parsed_args.input,
                     'clonemap': parsed_args.clonemap,
                     'outf': parsed_args.output,
                     'is_ldd': parsed_args.is_ldd}
    # MAIN METHOD
    convert(**configuration)


if __name__ == '__main__':
    main_script()
