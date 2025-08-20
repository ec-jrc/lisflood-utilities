import sys
import argparse
import logging
from datetime import datetime

import numpy as np
import pandas as pd
import geopandas as gpd

from lisfloodpreprocessing import Config, read_input_files
from lisfloodpreprocessing.utils import find_conflicts
from lisfloodpreprocessing.finer_grid import coordinates_fine
from lisfloodpreprocessing.coarser_grid import coordinates_coarse

logging.getLogger('pyogrio').propagate = False

def main():
    """
    Main function to correct point coordinates to match the river network in
    LISFLOOD static maps.
    """
    parser = argparse.ArgumentParser(
        description="""
        Correct the coordinates of a set of points to match the river network in the
        LISFLOOD static map.
        First, it uses a reference value of catchment area to find the most accurate 
        pixel in a high-resolution map.
        Second, it finds the pixel in the low-resolution map that better matches the 
        catchment shape derived from the high-resolution map.
        """
    )
    parser.add_argument(
        '-c', '--config-file', type=str, required=True, 
        help='Path to the configuration file'
    )
    args = parser.parse_args()

    # create the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # define a shared log format
    log_format = logging.Formatter(
        '%(asctime)s | %(levelname)s | %(name)s | %(message)s', 
        datefmt='%Y-%m-%d %H:%M:%S'
    )
        
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

    # main script logic
    success = False

    try:
        logger.info('Starting coordinate correction process...')
        
        # read configuration
        logger.info(f"Reading configuration from {args.config_file}")
        cfg = Config(args.config_file)
    
        # read input files
        logger.info('Reading input files...')
        inputs = read_input_files(cfg)      
    
        # find coordinates in high resolution
        logger.info('Processing points in the high-resolution grid...')
        points_HR, polygons_HR = coordinates_fine(
            cfg,
            points=inputs['points'],
            ldd_fine=inputs['ldd_fine'],
            upstream_fine=inputs['upstream_fine'],
            save=True
        )
    
        # find conflicts in high resolution
        logger.info('Finding conflicts in the high-resolution grid...')
        conflicts_fine = find_conflicts(
            points_HR,
            resolution=cfg.fine_resolution,
            pct_error=cfg.pct_error,
            save=cfg.output_folder / f'conflicts_{cfg.fine_resolution}.shp'
        )
        if conflicts_fine is not None:
            points_HR.drop(conflicts_fine.index, axis=0, inplace=True)
    
        # find coordinates in LISFLOOD
        logger.info('Processing points in the LISFLOOD grid...')
        points_LR, polygons_LR = coordinates_coarse(
            cfg,
            points_fine=points_HR,
            polygons_fine=polygons_HR,
            ldd_coarse=inputs['ldd_coarse'],
            upstream_coarse=inputs['upstream_coarse'],
            save=True
        )
    
        # find conflicts in LISFLOOD
        logger.info('Finding conflicts in the LISFLOOD grid...')
        conflicts_coarse = find_conflicts(
            points_LR,
            resolution=cfg.coarse_resolution,
            pct_error=cfg.pct_error,
            save=cfg.output_folder / f'conflicts_{cfg.coarse_resolution}.shp'
        )

        logger.info('Process completed successfully')
        success = True

    except FileNotFoundError as e:
        logger.error(f"A required file was not found: {e}")
    except OSError as e:
        logger.error(f"An I/O error occurred: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()