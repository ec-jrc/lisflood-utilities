"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

pcr2nc
A tool that convert PCRaster maps to a netCDF4 mapstack CF-1.7
"""

import sys
import argparse
import json

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from . import convert
from .. import version


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):
        required_group = self.add_argument_group('required arguments')
        required_group.add_argument('-m', '--metadata', help='Path to json/yaml metadata file for NetCDF',
                                    metavar='metadata', required=True)
        required_group.add_argument('-i', '--input', required=True,
                                    help='Path to input dataset. It can be a single PCRaster map, '
                                         'a folder containing PCRaster maps or a path with wildcards '
                                         'to use with glob.glob python function',
                                    metavar='input')
        self.add_argument('-o', '--output_file', help='Path to netcdf output file.', default='./file.nc',
                          metavar='output_file')


def parse_metadata(metadata_file):
    if metadata_file.endswith('.yaml') or metadata_file.endswith('.yml'):
        with open(metadata_file) as f:
            metadata = yaml.load(f, Loader=Loader)
    else:
        # suppose json format
        with open(metadata_file) as f:
            metadata = json.load(f)
    return metadata


def main_script():
    sys.exit(main(sys.argv[1:]))


def main(args):
    parser = ParserHelpOnError(description='Convert PCRaster to netCDF {}'.format(version))

    parser.add_args()
    parsed_args = parser.parse_args(args)
    configuration = {'dataset': parsed_args.input,
                     'output': parsed_args.output_file,
                     'metadata': parse_metadata(parsed_args.metadata)}
    # MAIN METHOD
    convert(**configuration)


if __name__ == '__main__':
    main_script()
