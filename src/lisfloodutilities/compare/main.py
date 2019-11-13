"""
Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

Usage: compare -a /workarea/datatests/modela_results/ -b /workarea/datatests/modelb_results/ -m /workarea/GLOFAS/maps/areamodel.nc
"""

import sys

import numpy as np

from . import NetCDFComparator, logger
from .helpers import ParserHelpOnError

np.set_printoptions(precision=4, linewidth=300, suppress=True)


def main(cliargs):
    parser = ParserHelpOnError(description='Compare netCDF outputs from two different LisFlood versions')
    parser.add_args()
    args = parser.parse_args(cliargs)

    dataset_a = args.dataset_a
    dataset_b = args.dataset_b
    maskfile = args.maskarea
    comparator = NetCDFComparator(maskfile)
    logger.info('\n\nComparing %s and %s\n\n ', dataset_a, dataset_b)
    errors = comparator.compare_dirs(dataset_a, dataset_b, )

    for i, e in enumerate(errors):
        logger.error('%d - %s', i, e)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    main_script()
