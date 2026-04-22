import sys
import argparse
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'basin-delineation_{datetime.now():%Y%m%d%H%M}.log')
    ]
)
for noisy_lib in ['pyogrio', 'rasterio', 'numba']:
    logging.getLogger(noisy_lib).setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

from . import Config
from .io import read_shapefiles
from .utils import find_conflicts
from .finer_grid import coordinates_fine
from coarser_grid import coordinates_coarse


def main():
    """
    Main function to correct point coordinates to match the river network in
    LISFLOOD static maps.
    """
    parser = argparse.ArgumentParser(
        description="""
        Correct the coordinates of a set of points to match the river network in a
        low-resolution DEM (digital elevation model).
            1. It uses a reference value of catchment area to find the most accurate
            pixel in a high-resolution map.
            2. It finds the pixel in the low-resolution map that better matches the
            catchment shape derived from the high-resolution map.
        If the low-resolution maps are not provided, the tool only does the first step.
        If the points and basins in the high-resolution are provided, together with 
        the low-resolution maps, the tool only does the second step.
        """
    )
    parser.add_argument(
        '-c', '--config-file', type=str, required=True, 
        help='Path to the configuration file'
    )
    args = parser.parse_args()

    # main script logic
    success = False

    try:
        logger.info('Starting coordinate correction process...')
        
        # read configuration
        logger.info(f"Reading configuration from {args.config_file}")
        cfg = Config(args.config_file)

        if cfg.run_fine:
            # find coordinates in high resolution
            logger.info('Processing points in the high-resolution grid...')
            points_fine, basins_fine = coordinates_fine(cfg, save=True)

            # find conflicts in high resolution
            logger.info('Finding conflicts in the high-resolution grid...')
            conflicts_fine = find_conflicts(
                points_fine,
                resolution=cfg.fine_resolution,
                pct_error=cfg.pct_error,
                save=cfg.output_folder / f'conflicts_{cfg.fine_resolution}.geojson'
            )
            if len(conflicts_fine) > 0:
                logger.warning(f"{len(conflicts_fine)} conflicts in the high-resolution grid")
                points_fine.drop(conflicts_fine.index, axis=0, inplace=True)
        else:
            # read points and basins in the finer resolution grid
            points_fine = read_shapefiles(cfg.points_fine)
            basins_fine = read_shapefiles(cfg.basins_fine)

        if cfg.run_coarse:
            # find coordinates in low resolution
            logger.info('Processing points in the low-resolution grid...')
            points_LR, polygons_LR = coordinates_coarse(cfg,
                                                        points_fine=points_fine,
                                                        polygons_fine=basins_fine,
                                                        save=True
                                                        )
        
            # find conflicts in low resolution
            logger.info('Finding conflicts in the low-resolution grid...')
            conflicts_coarse = find_conflicts(
                points_LR,
                resolution=cfg.coarse_resolution,
                pct_error=cfg.pct_error,
                save=cfg.output_folder / f'conflicts_{cfg.coarse_resolution}.geojson'
            )
            if len(conflicts_coarse) > 0:
                logger.warning(f"{len(conflicts_coarse)} conflicts in the low-resolution grid")
        
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