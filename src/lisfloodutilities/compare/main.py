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
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from .helpers import logger, ParserHelpOnError, check_vars

np.set_printoptions(precision=4, linewidth=300, suppress=True)


def main(cliargs):
    parser = ParserHelpOnError(description='Compare netCDF outputs from two different LisFlood versions')
    parser.add_args()
    args = parser.parse_args(cliargs)

    dataset_a = args.dataset_a
    dataset_b = args.dataset_b
    maskfile = args.maskarea
    masknc = Dataset(maskfile)
    maskvar = masknc.variables['areaModel']
    maskarea = np.logical_not(maskvar[:, :])

    logger.info('\n\nComparing %s and %s\n\n ', dataset_a, dataset_b)
    errors = []
    for fa in Path(dataset_a).glob('**/*.nc'):
        fb = os.path.join(dataset_b, os.path.basename(fa))
        if not os.path.exists(fb):
            logger.info('skipping %s as it is not in %s', fb, dataset_b)
            continue
        nca = Dataset(fa)
        ncb = Dataset(fb)

        num_dims = 3 if 'time' in nca.variables else 2
        var_name = [k for k in nca.variables if len(nca.variables[k].dimensions) == num_dims][0]
        vara = nca.variables[var_name]
        varb = ncb.variables[var_name]
        if 'time' in nca.variables:
            stepsa = nca.variables['time'][:]
            stepsb = ncb.variables['time'][:]
            if len(stepsa) != len(stepsb):
                logger.error('Different number of timesteps')
                errors.append('File: {}: different size in time axis A:{} B:{}'.format(os.path.basename(fa), len(stepsa), len(stepsb)))
                continue
            for i, step in enumerate(stepsa):
                vara_step = vara[i][:, :]
                varb_step = varb[i][:, :]
                check_vars(maskarea, i, vara_step, varb_step, var_name, fa, errors)
        else:
            vara_step = vara[:, :]
            varb_step = varb[:, :]
            check_vars(maskarea, None, vara_step, varb_step, var_name, fa, errors)

    for i, e in enumerate(errors):
        logger.error('%d - %s', i, e)


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
