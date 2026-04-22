import logging
from pathlib import Path
from typing import Literal, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
import pyflwdir
import rioxarray as rxr

from . import Config

logger = logging.getLogger(__name__)


def read_flowdir(
    path: Union[str, Path], 
    ftype: Optional[Literal['d8', 'ldd', 'nextxy', 'infer']] = 'infer'
) -> Tuple[pyflwdir.FlwdirRaster, str]:
    """
    Reads a flow direction raster and converts it into a PyFlwdir object.

    This function utilizes rioxarray to open a geospatial raster, ensures it is 
    a 2D array by squeezing the band dimension, and initializes a 
    pyflwdir.FlwdirRaster object with the associated affine transform.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the flow direction raster file (e.g., .tif, .asc).
    ftype : {'d8', 'ldd', 'nextxy', 'infer'}, optional
        The flow direction type/convention of the input raster. 
        If 'infer' (default), the function attempts to determine the type 
        based on the data values.

    Returns
    -------
    flwdir : pyflwdir.FlwdirRaster
        The initialized PyFlwdir flow direction object, ready for 
        topological operations.
    resolution : str
        A string representing the grid resolution (e.g., '30sec' or '1min').
        The unit is 'min' if the cell size is greater than 60 arcseconds, 
        otherwise 'sec'.

    Notes
    -----
    The function assumes the input raster is in a geographic coordinate 
    system (`latlon=True`) and disables `check_ftype` to speed up 
    initialization if the type is already known.
    """

    # read flow direction map
    # Pylance incorrectly infers open_rasterio returns list[Dataset]; it actually returns RasterioArrayAdapter (DataArray)
    raster = rxr.open_rasterio(path).squeeze(dim='band')  # type: ignore[union-attr]

    # Convert to float32 to ensure compatibility with pyflwdir
    # (pyflwdir requires one of: int16, int32, uint8, uint16, float32, float64, int8)
    raster_data = raster.data.astype(np.float32)

    # convert into a pyflwdir object
    flwdir = pyflwdir.from_array(
            raster_data,
            ftype=ftype,
            transform=raster.rio.transform(),
            check_ftype=False,
            latlon=True
        )
    
    # raster resolution
    cellsize_deg = abs(flwdir.transform[0]) # degrees
    cellsize = int(np.round(cellsize_deg * 3600, 0)) # arcsec
    unit = 'sec'
    if cellsize >= 60:
        cellsize = int(np.round(cellsize_deg * 60, 0)) # arcmin
        unit = 'min'
    resolution = f'{cellsize}{unit}'

    return flwdir, resolution


def read_shapefiles(path: Union[str, Path]) -> gpd.GeoDataFrame:
    """
    Reads a geospatial file into a GeoDataFrame and standardizes the format.

    This function loads a vector file (e.g., Shapefile, GeoJSON), converts all 
    column names to lowercase for consistency, sets the 'id' column as the 
    index, and ensures the index is of integer type.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the geospatial file to be read.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with a lowercase schema and an integer-based 'id' index.
    """

    points = gpd.read_file(path)
    points.columns = points.columns.str.lower()
    points.set_index('id', inplace=True)
    points.index = points.index.astype(int)

    return points


def save_results(
        cfg: Config,
        points: gpd.GeoDataFrame,
        polygons: gpd.GeoDataFrame,
        resolution: str
) -> None:
    """
    Exports the resulting outlet points and catchment polygons to GeoJSON files.

    The filenames are constructed using a base name (stem) derived from the 
    original input points, suffixed with the specific grid resolution. The 
    function logs the successful export of both layers to the output folder 
    defined in the configuration.

    Parameters
    ----------
    cfg : Config
        Configuration object containing output folder paths and input file stems.
    points : gpd.GeoDataFrame
        The GeoDataFrame containing the updated outlet point coordinates.
    polygons : gpd.GeoDataFrame
        The GeoDataFrame containing the delineated catchment boundaries.
    resolution : str
        A string indicating the grid resolution (e.g., '1min', '30sec') 
        used for file naming and logging.

    Returns
    -------
    None
        The function saves files directly to the disk and does not return a value.
    """

    if cfg.points is not None:
        stem = cfg.points.stem
    else:
        # cfg.points_fine is guaranteed to be not None when cfg.points is None
        # (validated in Config.__init__)
        assert cfg.points_fine is not None
        stem = cfg.points_fine.stem.split('_')[0]

    if cfg.output_folder is None:
        raise ValueError("Output folder path is not defined in the configuration.")

    # polygons
    polygon_file = cfg.output_folder / f'{stem}_basins_{resolution}.geojson'
    polygons.to_file(polygon_file, driver='GeoJSON')
    logger.info(f'Basins at {resolution} resolution exported to: {polygon_file}')
    
    # points
    point_file = cfg.output_folder / f'{stem}_outlets_{resolution}.geojson'
    points.to_file(point_file, driver='GeoJSON')
    logger.info(f'Outlets at {resolution} resolution exported to: {point_file}')