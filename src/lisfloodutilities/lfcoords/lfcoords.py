import pandas as pd
import argparse
import logging

from lisfloodpreprocessing import Config
from lisfloodpreprocessing.finer_grid import coordinates_fine
from lisfloodpreprocessing.coarser_grid import coordinates_coarse

def main():
    
    parser = argparse.ArgumentParser(
        description="""
        Correct the coordinates of a set of stations to match the river network in the
        LISFLOOD static map.
        First, it uses a reference value of catchment area to find the most accurate 
        pixel in a high-resolution map.
        Second, it finds the pixel in the low-resolution map that better matches the 
        catchment shape derived from the high-resolution map.
        """
    )
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Path to the configuration file')
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
    cfg = Config(args.config_file)
    
    # find coordinates in high resolution
    stations_HR = coordinates_fine(cfg, save=False)
    
    # find coordinates in LISFLOOD
    coordinates_coarse(cfg, stations_HR, save=True)

if __name__ == "__main__":
    main()