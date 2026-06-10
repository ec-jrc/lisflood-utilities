"""
Copyright 2019-2026 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import argparse
import sys
import numpy as np
import os
import time
import logging
from netCDF4 import Dataset
from pathlib import Path
from typing import Optional

from ..cutmaps.helpers import verify_existing_netcdf

logging.basicConfig(format='[%(asctime)s] - %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
logger = logging.getLogger()


def main(cliargs):
    parser = ParserHelpOnError(description='Verify that the Water Regions map is consistent with the map of the calibration catchments')
    parser.add_args()
    args = parser.parse_args(cliargs)
    calib_catchments = args.calib_catchments 
    waterregions = args.waterregions 
    output_message = verify_waterregions(calib_catchments, waterregions, print_progress=True)
    logger.info('\n Verify the consistency between %s and %s ', calib_catchments, waterregions)
    logger.info(output_message)


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)
        
    def add_args(self):
        self.add_argument("-cc", "--calib_catchments",
                          help='map of calibration catchments, netcdf format',
                          required=True, type=Path)
        self.add_argument("-wr", "--waterregions",
                          help='map of water regions, netcdf format',
                          required=True, type=Path)


def verify_waterregions(calib_catchments: Path , waterregions: Path, print_progress: bool = False) -> str:
    # Check the input files
    msg = verify_existing_netcdf(input_file=calib_catchments, file_id='calibration catchments')
    if len(msg) > 0:
        return msg
    msg = verify_existing_netcdf(input_file=waterregions, file_id='water regions')
    if len(msg) > 0:
        return msg

    waterregions_map = Dataset(waterregions, 'r', format='NETCDF4_CLASSIC')
    ywr = []
    wr_yresolution = 0
    xwr = []
    wr_xresolution = 0
    wr = []
    for v in waterregions_map.variables:
        if (v == 'y' or v == 'lat'):
           ywr = waterregions_map.variables[v][:]
           wr_yresolution = np.abs((ywr[-1] - ywr[0])/(len(ywr[:])-1.0))
        if (v == 'x' or v == 'lon'):
           xwr = waterregions_map.variables[v][:]
           wr_xresolution = np.abs((xwr[-1] - xwr[0])/(len(xwr[:])-1.0))
        if len(waterregions_map.variables[v].dimensions) == 2:
           wr = waterregions_map.variables[v][:]

    cal_catch = np.array([])
    flipud = False
    calibration_catchments_map = Dataset(calib_catchments, 'r', format='NETCDF4_CLASSIC')
    for v in calibration_catchments_map.variables:
        if (v == 'y' or v == 'lat'):
           ycc = calibration_catchments_map.variables[v][:]
           cc_resolution = np.abs((ycc[-1] - ycc[0])/(len(ycc[:])-1.0))
           if round(ycc[0],1) != round(ywr[0],1):
              flipud = True
              check_grid_y = np.amax(np.abs(np.flip(ycc)-ywr))
              if (check_grid_y > wr_yresolution):
                   return 'Error: the two maps do not have the same coordinates. Please, check the input maps.'
        if (v == 'x' or v == 'lon'):
           xcc = calibration_catchments_map.variables[v][:]
           check_grid_x = np.amax(np.abs(xcc-xwr))
           if (check_grid_x > wr_xresolution):
                   return 'Error: the two maps do not have the same coordinates. Please, check the input maps.'
        if len(calibration_catchments_map.variables[v].dimensions) == 2:
           cal_catch = calibration_catchments_map.variables[v][:]
           if flipud:
              print('Warning: one of the maps has inverted y-axis, is this an intended feature?')
              cal_catch = np.flipud(cal_catch)

    cal_catch[cal_catch < 1] = -9999
    cal_catch[np.isnan(cal_catch) == 1] = -9999
    wr_ids = np.unique(wr)

    id_error_wr = []
    cal_catch_error_wr = []
    output_message = ''

    for idx, wr_id in enumerate(wr_ids):
        if print_progress:
            progress = (idx + 1) / len(wr_ids) * 100
            print(f'\rProgress: {progress:.1f}%', end='')
        extract_wr = np.where(wr == wr_id, cal_catch, -9999)
        num_cal_catch = np.unique(extract_wr)
        num_cal_catch_check = np.extract(num_cal_catch != -9999, num_cal_catch)
        if len(num_cal_catch_check) > 1:
            id_error_wr.append(int(wr_id))
            num_cal_catch_write = num_cal_catch[num_cal_catch != -9999].astype(int)
            cal_catch_error_wr.append(num_cal_catch_write)
            msg = 'The water regions WR are included in more than one calibration catchment'
            msg_catchments = [list(c) for c in cal_catch_error_wr]
            output_message = f'ERROR: {msg}\nWR={id_error_wr}\nCalibration catchments={msg_catchments}'
    if len(id_error_wr) == 0:
        output_message = 'PASSED: Each calibration catchment contains only a finite number of water regions.'

    waterregions_map.close()
    calibration_catchments_map.close()

    return output_message


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    main_script()                         

