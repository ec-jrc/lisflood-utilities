"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

Usage: compare -a /workarea/datatests/modela_results/ -b /workarea/datatests/modelb_results/ -m /workarea/GLOFAS/maps/areamodel.nc
"""

import argparse
import sys
from nine import IS_PYTHON2
if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

import numpy as np

from lisfloodutilities.compare.nc import NetCDFComparator
from .. import version, logger

np.set_printoptions(precision=4, linewidth=300, suppress=True)


def main(cliargs):
    parser = ParserHelpOnError(description='Compare netCDF outputs: {}'.format(version))
    parser.add_args()
    args = parser.parse_args(cliargs)

    dataset_a = args.dataset_a.as_posix()
    dataset_b = args.dataset_b.as_posix()
    maskfile = args.maskarea.as_posix()

    array_equal = args.array_equal
    skip_missing = args.skip_missing
    atol = args.atol
    rtol = args.rtol
    max_diff_perc = args.max_diff_percentage
    max_large_diff_perc = args.max_largediff_percentage
    save_diff_files = args.save_diffs
    comparator = NetCDFComparator(maskfile, array_equal=array_equal, for_testing=False, atol=atol, rtol=rtol,
                                  max_perc_diff=max_diff_perc, max_perc_large_diff=max_large_diff_perc,
                                  save_diff_files=save_diff_files)
    logger.info('\n\nComparing %s and %s\n\n ', dataset_a, dataset_b)
    errors = comparator.compare_dirs(dataset_a, dataset_b, skip_missing=skip_missing)

    for i, e in enumerate(errors):
        logger.error('%d - %s', i, e)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    main_script()


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):

        self.add_argument("-a", "--dataset_a",
                          type=lambda p: Path(p),
                          help='path to dataset version A',
                          required=True)
        self.add_argument("-b", "--dataset_b",
                          type=lambda p: Path(p),
                          help='path to dataset version B',
                          required=True)
        self.add_argument("-m", "--maskarea",
                          type=lambda p: Path(p),
                          help='path to mask',
                          required=True)
        self.add_argument("-e", "--array-equal",
                          help='flag to compare files to be identical',
                          default=False, action='store_true')
        self.add_argument("-s", "--skip-missing",
                          help='flag to skip missing files in comparison',
                          default=False, action='store_true')
        self.add_argument("-D", "--save-diffs",
                          help='flag to save diffs in netcdf files for visual comparisons. '
                               'Files are saved in ./diffs folder of current directory.'
                               'For each file presenting differences, you will find files diffs, '
                               'original A and B (only for timesteps where differences are found).',
                          required=False, default=False, action='store_true')
        self.add_argument('-r', '--rtol', help='rtol', type=float, default=0.001)
        self.add_argument('-t', '--atol', help='atol', type=float, default=0.0001)
        self.add_argument('-p', '--max-diff-percentage', help='threshold for diffs percentage', type=float, default=0.2)
        self.add_argument('-l', '--max-largediff-percentage', help='threshold for large diffs percentage', type=float, default=0.1)
