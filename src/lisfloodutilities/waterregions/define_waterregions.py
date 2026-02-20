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

import argparse
import shutil
import sys
import numpy as np
import pcraster as pcr
import os
import time
import logging
from lisfloodutilities.nc2pcr import convert as convnc2pcr
from lisfloodutilities.pcr2nc import convert as convpcr2nc
import tempfile
from pathlib import Path
from osgeo import gdal
from netCDF4 import Dataset
from typing import List, Optional, Tuple, Union

from ..cutmaps.helpers import (verify_existing_netcdf, col2netcdf,
                               array_to_nc_from_clone, get_river_network_from_map)
from earthkit.hydro import catchments

import yaml
import json

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

logging.basicConfig(format='[%(asctime)s] - %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
logger = logging.getLogger()


def main(cliargs):
    parser = ParserHelpOnError(description='Define Water Regions consistent with calibration points: {}')
    parser.add_args()
    args = parser.parse_args(cliargs)

    calib_points = Path(args.calib_points)
    countries_id = Path(args.countries_id)
    ldd = Path(args.ldd)
    waterregions_initial = Path(args.waterregions_initial)
    output_wr = Path(args.output_wr)
    metadata_file = Path(args.metadata_file)
    
    if output_wr.suffix == '.nc':
       metadata_parsed = parse_metadata(metadata_file)
    else:
       metadata_parsed = {}

    if output_wr is not None and output_wr.parent.exists():
        tmp_folder_path = Path(output_wr.parent, "tmp")
    else:
        tmp_folder_path = Path(tempfile.gettempdir(), "lisflood_tmp")
    tmpToRemove=False
    if not os.path.exists(tmp_folder_path):
         tmpToRemove=True
         os.mkdir(tmp_folder_path)

    waterregion_nc, waterregion_unique, subcatchment = define_waterregions(calib_points, countries_id,
                                                                           ldd, waterregions_initial,
                                                                           output_wr, metadata_parsed,
                                                                           tmpdir=tmp_folder_path)    
    if tmpToRemove:
         if os.path.exists(tmp_folder_path):
            try:
               shutil.rmtree(tmp_folder_path)
            except:
                pass
    
    if waterregion_nc is None or waterregion_unique is None or subcatchment is None:
        logger.error('Water regions definition failed due to input errors. Please, check the error messages above.')
        sys.exit(1)

    logger.info('\nUsing %s and %s to define the water regions\n ', calib_points, ldd)
    
    # check the consistency between the water regions and the calibration catchments
    cal_catch = subcatchment
    cal_catch[cal_catch < 1] = -9999
    cal_catch[np.isnan(cal_catch)] = -9999

    wr_ids = np.unique(waterregion_unique)

    id_error_wr = []
    cal_catch_error_wr = []
    output_message = ''

    for wr_id in wr_ids:
           extract_wr = np.where(waterregion_unique == wr_id, cal_catch, -9999)
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
             output_message = 'OK! Each water region is completely included inside one calibration catchment.'
    
    print(output_message)


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)
        
    def add_args(self):

        self.add_argument("-p", "--calib_points",
                          help='list of calibration points: lon or x, lat or y, point id',
                          required=True, type=Path)
        self.add_argument("-l", "--ldd",
                          help='LDD map netCDF format',
                          required=True, type=Path)
        self.add_argument("-C", "--countries_id",
                          help='map of Countries ID netCDF format',
                          required=True, type=Path)
        self.add_argument("-w", "--waterregions_initial",
                          help='initial map of water regions netCDF format',
                          required=True, type=Path)
        self.add_argument("-o", "--output_wr",
                          help='output map of water regions netCDF format',
                          required=True, type=Path)
        self.add_argument("-m","--metadata_file",                        
                          help='Path to metadata file for netCDF, .yaml or .json format',
                          required=False, type=Path)


def parse_metadata(metadata_file: Path) -> dict:
    if metadata_file.suffix == '.yaml' or metadata_file.suffix == '.yml':
        with open(metadata_file) as f:
            metadata = yaml.load(f, Loader=Loader)
    else:
        # suppose json format
        with open(metadata_file) as f:
            metadata = json.load(f)
    return metadata                            


def verify_inputs(calib_points: Union[Path, str], countries_id: Union[Path, str],
                  ldd: Union[Path, str], waterregions_initial: Union[Path, str],
                  output_wr: Union[Path, str]) -> List[str]:
    msgs = []
    calib_points_path = calib_points if isinstance(calib_points, Path) else Path(calib_points)
    # Guarantee the calib_points exists
    if not calib_points_path.exists() or not calib_points_path.is_file():
        msgs.append('Error: The calibration points must be an existing file')
    # Guarantee the countries_id exists and is netCDF
    msg = verify_existing_netcdf(input_file=countries_id, file_id='countries ids')
    if len(msg) > 0:
        msgs.append(msg)
    # Guarantee the ldd exists and is netCDF
    msg = verify_existing_netcdf(input_file=ldd, file_id='LDD')
    if len(msg) > 0:
        msgs.append(msg)
    # Guarantee the waterregions_initial exists and is netCDF
    msg = verify_existing_netcdf(input_file=waterregions_initial, file_id='water regions initial')
    if len(msg) > 0:
        msgs.append(msg)
    # Check the output water regions file
    output_wr_path = output_wr if isinstance(output_wr, Path) else Path(output_wr)
    if (output_wr_path is None or output_wr_path.suffix != '.nc'):
        msgs.append('Error: The output file must be a netCDF file.')
    return msgs


def get_data_from_timeless_nc(nc_file: Union[Path, str], var_name: Optional[str] = None) -> np.ndarray:
    nc_file_path = nc_file if isinstance(nc_file, Path) else Path(nc_file)
    with Dataset(nc_file_path, 'r', 'format=NETCDF4_classic') as nc:
        if var_name is None:
            var_name = str([v for v in nc.variables if len(nc.variables[v].dimensions) == 2][0])
        data = nc[var_name][:]
    return data

def define_waterregions(calib_points: Path, countries_id: Path, ldd: Path,
                        waterregions_initial: Path, output_wr: Path, metadata_parsed: dict = {},
                        tmpdir: Union[Path, None] = None) ->  Tuple[Union[Path, None],
                                                                    Union[np.ndarray, None],
                                                                    Union[np.ndarray, None]]:
    if tmpdir is None:
        tmpdir = Path(tempfile.gettempdir())

    # #0. Guarantee the input files exist and are netCDF
    msg = verify_inputs(calib_points, countries_id, ldd, waterregions_initial, output_wr)
    if len(msg) > 0:
        print('\n'.join(msg))
        return None, None, None
          
    #1. The calibration points are converted into a map
    pointmap_file = Path(tmpdir, 'points.nc')
    try:
       os.remove(pointmap_file)
    except:
       pass
    points = col2netcdf(calib_points, pointmap_file, countries_id, quiet=False)

    country_ids_map = get_data_from_timeless_nc(countries_id)
    waterregions_initial_map = get_data_from_timeless_nc(waterregions_initial)

    #2. Catchment map1 derived from calibration points
    network = get_river_network_from_map(ldd)
    subcatchments = catchments.find(network, points)
    try:
        os.remove(pointmap_file)
    except:
        pass

    #3. Making map with all valid ldd cells
    #4. Taking only the catchments from calibration points
    land = np.where(subcatchments > 0, subcatchments, np.nan)

    #5. Splitting riverbasins if part of different countries, so if cross-border catchments
    # Guaranty the resulting ids are unique by multiplying by 2 different prime numbers and then adding
    subcatchments_unique = (5 * land) + (7 * country_ids_map)

    #6. Covering the missing areas with the old waterregions
    #8. Each water region is given a unique ID number
    waterregions_unique = np.where(np.isnan(subcatchments_unique),
                                   subcatchments_unique + (11 * waterregions_initial_map),
                                   subcatchments_unique)

    #9. Save the water region map
    array_to_nc_from_clone(out_path=output_wr, clone_path=ldd, grid=waterregions_unique, metadata=metadata_parsed)

    print(waterregions_unique)

    return output_wr, waterregions_unique, subcatchments
                          
def main_script():
    gdal.UseExceptions()
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    main_script()                         

