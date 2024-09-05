import numpy as np
import pandas as pd
import rioxarray
import pyflwdir
from tqdm import tqdm
from pathlib import Path
from typing import Optional
import logging
import warnings
warnings.filterwarnings("ignore")

from lisfloodutilities.lfcoords import Config
from lisfloodutilities.lfcoords.utils import find_pixel, catchment_polygon 

# set logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def coordinates_fine(
    cfg: Config,
    save: bool = True
) -> Optional[pd.DataFrame]:
    """
    Processes point coordinates to find the most accurate pixel in a high-resolution map, 
    based on a reference value of catchment area. It updates the point coordinates and 
    exports the catchment areas as shapefiles.

    The function reads the upstream area map and local drainage direction (LDD) map in 
    fine resolution, as well as the point coordinates with their reference upstream areas.
    It then iteratively searches for the best-matching pixel within an increasing search 
    range and applies penalties and factors to determine the closest match. The function 
    creates a boolean map of each point's catchment area and vectorizes it into a 
    GeoDataFrame for export as a shapefile.

    Parameters:
    -----------
    cfg: Config
        Configuration object containing file paths and parameters specified in  the configuration file.
    save: bool
        If True, the updated points table is saved to a CSV file. If False, the updated points DataFrame is returned without saving.
        
    Returns:
    --------
    points: pandas.DataFrame
        If save is False, returns a pandas DataFrame with updated point coordinates and upstream areas in the finer grid. Otherwise, the function returns None, and the results are saved directly to a CSV file.
    """
    
    # READ INPUTS
    
    # read upstream map with fine resolution
    upstream_fine = rioxarray.open_rasterio(cfg.UPSTREAM_FINE).squeeze(dim='band')
    logger.info(f'Map of upstream area corretly read: {cfg.UPSTREAM_FINE}')

    # read local drainage direction map
    ldd_fine = rioxarray.open_rasterio(cfg.LDD_FINE).squeeze(dim='band')
    logger.info(f'Map of local drainage directions corretly read: {cfg.LDD_FINE}')

    # read points text file
    points = pd.read_csv(cfg.POINTS, index_col='ID')
    logger.info(f'Table of points correctly read: {cfg.POINTS}')
    

    # PROCESSING
    
    # resolution of the input map
    cellsize = np.mean(np.diff(upstream_fine.x)) # degrees
    cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
    suffix_fine = f'{cellsize_arcsec}sec'
    logger.info(f'The resolution of the finer grid is {cellsize_arcsec} arcseconds')
    
    # add columns to the table of points
    new_cols = sorted([f'{col}_{suffix_fine}' for col in points.columns])
    points[new_cols] = np.nan

    # create river network
    fdir_fine = pyflwdir.from_array(ldd_fine.data,
                                    ftype='d8', 
                                    transform=ldd_fine.rio.transform(),
                                    check_ftype=False,
                                    latlon=True)

    # output path
    SHAPE_FOLDER_FINE = cfg.SHAPE_FOLDER / suffix_fine
    SHAPE_FOLDER_FINE.mkdir(parents=True, exist_ok=True)

    for ID, attrs in tqdm(points.iterrows(), total=points.shape[0], desc='points'):  

        # reference coordinates and upstream area
        lat_ref, lon_ref, area_ref = attrs[['lat', 'lon', 'area']]

        # search new coordinates in an increasing range
        ranges = [55, 101, 151]
        penalties = [500, 500, 1000]
        factors = [2, .5, .25]
        acceptable_errors = [50, 80, np.nan]
        for rangexy, penalty, factor, max_error in zip(ranges, penalties, factors, acceptable_errors):
            logger.debug(f'Set range to {rangexy}')
            lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=penalty, factor=factor)
            if error <= max_error:
                break

        # update new columns in 'points'
        points.loc[ID, new_cols] = [int(upstream_fine.sel(y=lat, x=lon).item()), round(lat, 6), round(lon, 6)]

        # boolean map of the catchment associated to the corrected coordinates
        basin_arr = fdir_fine.basins(xy=(lon, lat)).astype(np.int32)

        # vectorize the boolean map into geopandas
        basin_gdf = catchment_polygon(basin_arr.astype(np.int32),
                                      transform=ldd_fine.rio.transform(),
                                      crs=ldd_fine.rio.crs,
                                      name='ID')
        basin_gdf['ID'] = ID
        basin_gdf.set_index('ID', inplace=True)
        basin_gdf[attrs.index] = attrs.values

        # export shape file
        output_file = SHAPE_FOLDER_FINE / f'{ID}.shp'
        basin_gdf.to_file(output_file)
        logger.info(f'Catchment {ID} exported as shapefile: {output_file}')
    
    # return/save
    points.sort_index(axis=1, inplace=True)
    if save:
        output_csv = cfg.POINTS.parent / f'{cfg.POINTS.stem}_{suffix_fine}.csv'
        points.to_csv(output_csv)
        logger.info(f'The updated points table in the finer grid has been exported to: {output_csv}')
    else: 
        return points