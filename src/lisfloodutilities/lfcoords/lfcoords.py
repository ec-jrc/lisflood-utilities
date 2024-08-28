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

from lisfloodutilities.lfcoords import Config
from lilisfloodutilities.lfcoords.finer_grid import coordinates_fine
from llisfloodutilities.lfcoords.coarser_grid import coordinates_coarse

def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Correct the coordinates of a set of stations to match the river network in the
        LISFLOOD static map.
        First, it uses a reference value of catchment area to find the most accurate
        pixel in a high-resolution map.
        Second, it finds the pixel in the low-resolution map that better matches the
        catchment shape derived from the high-resolution map.
        """,
        prog=prog
    )
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Path to the YML configuration file')
    args = parser.parse_args()

    # create logger
    logger = logging.getLogger('correct-coordinates')
    logger.setLevel(logging.INFO)
    logger.propagate = False
    log_format = logging.Formatter('%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    c_handler = logging.StreamHandler()
    c_handler.setFormatter(log_format)
    c_handler.setLevel(logging.INFO)
    logger.addHandler(c_handler)

    # read configuration
    try:
        cfg = Config(args.config_file)
    except Exception as e:
        logger.error(f'Reading the config file: {e}')
        sys.exit(1)

    # find coordinates in high resolution
    try:
        stations_HR = coordinates_fine(cfg, save=False)
    except Exception as e:
        logger.error(f'Locating the points in the finer grid: {e}')
        sys.exit(2)

    # find coordinates in LISFLOOD
    try:
        coordinates_coarse(cfg, stations_HR, save=True)
    except Exception as e:
        logger.error(f'Locating the points in the finer grid: {e}')
        sys.exit(3)

def main_script():
    sys.exist(main())

if __name__ == "__main__":
    main_script()
