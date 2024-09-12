"""
Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import os
import sys
import argparse
import pandas as pd
import logging
from datetime import datetime

from lisfloodpreprocessing import Config, read_input_files
from lisfloodpreprocessing.utils import find_conflicts
from lisfloodpreprocessing.finer_grid import coordinates_fine
from lisfloodpreprocessing.coarser_grid import coordinates_coarse

def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Correct the coordinates of a set of points to match the river network in the
        LISFLOOD static map.
        First, it uses a reference value of catchment area to find the most accurate
        pixel in a high-resolution map.
        Second, it finds the pixel in the low-resolution map that better matches the
        catchment shape derived from the high-resolution map.
        """,
        prog=prog
    )
    parser.add_argument('-c', '--config-file', type=str, required=True,
                        help='Path to the YML configuration file')
    parser.add_argument('-r', '--reservoirs', action='store_true', default=False,
                        help='The input points are reservoirs')
    args = parser.parse_args()

    # create logger
    logger = logging.getLogger('lfcoords')
    logger.setLevel(logging.INFO)
    logger.propagate = False
    log_format = logging.Formatter('%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    # console handler
    c_handler = logging.StreamHandler()
    c_handler.setFormatter(log_format)
    c_handler.setLevel(logging.INFO)
    logger.addHandler(c_handler)
    # File handler
    f_handler = logging.FileHandler(f'lfcoords_{datetime.now():%Y%m%d%H%M}.log')
    f_handler.setFormatter(log_format)
    f_handler.setLevel(logging.INFO)
    logger.addHandler(f_handler)

    # read configuration
    try:
        cfg = Config(args.config_file)
    except Exception as e:
        logger.error(f'Reading the config file: {e}')
        sys.exit(1)
    
    # read input files
    try:
        inputs = read_input_files(cfg)
    except Exception as e:
        logger.error(f'Reading the input files: {e}')
        sys.exit(2)       

    # find coordinates in high resolution
    try:
        points_HR, polygons_HR = coordinates_fine(cfg,
                                                  points=inputs['points'],
                                                  ldd_fine=inputs['ldd_fine'],
                                                  upstream_fine=inputs['upstream_fine'],
                                                  save=True)
    except Exception as e:
        logger.error(f'Locating the points in the finer grid: {e}')
        sys.exit(3)
    
    # find conflicts in high resolution
    try:
        find_conflicts(points_HR,
                       columns=[f'{var}_{cfg.FINE_RESOLUTION}' for var in ['lat', 'lon']],
                       save=cfg.OUTPUT_FOLDER / f'conflicts_{cfg.FINE_RESOLUTION}.shp')
    except Exception as e:
        logger.error(f'Finding conflicts in the finer grid: {e}')
        sys.exit(4)

    # find coordinates in LISFLOOD
    try:
        points_LR, polygons_LR = coordinates_coarse(cfg,
                                                    points_fine=points_HR,
                                                    polygons_fine=polygons_HR,
                                                    ldd_coarse=inputs['ldd_coarse'],
                                                    upstream_coarse=inputs['upstream_coarse'],
                                                    reservoirs=args.reservoirs,
                                                    save=True)
    except Exception as e:
        logger.error(f'Locating the points in the finer grid: {e}')
        sys.exit(5)
        
    # find conflicts in LISFLOOD
    try:
        find_conflicts(points_LR,
                       columns=[f'{var}_{cfg.COARSE_RESOLUTION}' for var in ['lat', 'lon']],
                       save=cfg.OUTPUT_FOLDER / f'conflicts_{cfg.COARSE_RESOLUTION}.shp')
    except Exception as e:
        logger.error(f'Finding conflicts in the coarser grid: {e}')
        sys.exit(6)

def main_script():
    sys.exit(main())

if __name__ == "__main__":
    main_script()
