import pytest
from pathlib import Path
import geopandas as gpd
from geopandas.testing import assert_geodataframe_equal

from lisfloodutilities.lfcoords import Config
from lisfloodutilities.lfcoords.io import read_shapefiles
from lisfloodutilities.lfcoords.finer_grid import coordinates_fine
from lisfloodutilities.lfcoords.coarser_grid import coordinates_coarse

@pytest.fixture
def base_path():
    """Provides the path to the test data directory."""
    return Path(__file__).parent / 'data' / 'lfcoords' / 'basin-delineation'

@pytest.fixture
def cfg(base_path):
    """Initializes the Config object once for the test session."""
    config_path = base_path / 'config.yml'
    return Config(config_path)

def test_basin_delineation(base_path, cfg):
    """
    Validates that both fine and coarse delineation match the 
    pre-computed GeoJSONs.
    """

    # Compute test values
    test_fine, basins_fine = coordinates_fine(cfg, save=False)
    test_coarse, _ = coordinates_coarse(cfg,
                                        test_fine,
                                        basins_fine,
                                        save=False)

    # Load expected values
    stem = cfg.points.stem
    expected_fine_path = base_path / f'{stem}_outlets_{cfg.fine_resolution}.geojson'
    expected_fine = read_shapefiles(expected_fine_path)
    expected_coarse_path = base_path / f'{stem}_outlets_{cfg.coarse_resolution}.geojson'
    expected_coarse = read_shapefiles(expected_coarse_path)

    # Spatial Assertions
    assert_geodataframe_equal(
        test_fine,
        expected_fine,
        check_less_precise=True,
        check_like=True,
        check_dtype=False
    )
    assert_geodataframe_equal(
        test_coarse,
        expected_coarse,
        check_less_precise=True,
        check_like=True,
        check_dtype=False
    )
