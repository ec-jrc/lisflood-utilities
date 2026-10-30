"""Integration tests for the ldd2vec tool.

Each test case converts an input LDD netCDF into a vector line network and
compares the produced ESRI Shapefile against a reference shapefile. The
comparison is order-independent (features may be written in any order) and
checks the attribute value, the geometry, and the CRS of every feature.

Adding coverage for another LDD input is intentionally cheap: drop the input
netCDF under ``tests/data/ldd2vec/input/`` and its expected shapefile set under
``tests/data/ldd2vec/reference/``, then append one ``LddVecCase`` entry to
``CASES`` below. No new test function is required.
"""

import os
from collections import Counter
from dataclasses import dataclass, field
from typing import List, Optional

import geopandas as gpd
import pytest
from shapely import wkt

from lisfloodutilities.ldd2vec.ldd2vec import ldd_to_vector


DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "ldd2vec")
INPUT_DIR = os.path.join(DATA_DIR, "input")
REFERENCE_DIR = os.path.join(DATA_DIR, "reference")

# Geometry coordinates are compared after rounding to this many decimal places
# so that writing/reading a shapefile (which stores doubles) does not introduce
# spurious mismatches. Coordinates are grid cell centres/corners, so this is far
# finer than the grid resolution.
_COORD_PRECISION = 9


@dataclass
class LddVecCase:
    """A single input-to-reference conversion case.

    :ivar name: human-readable id used for the parametrized test id.
    :ivar input_name: file name of the input LDD netCDF under ``input/``.
    :ivar reference_name: file name of the reference ``.shp`` under ``reference/``.
    :ivar river_network_format: river network format passed to ``ldd_to_vector``.
    :ivar coords_names: optional explicit ``[y, x]`` coordinate names; empty to
        auto-detect.
    """

    name: str
    input_name: str
    reference_name: str
    river_network_format: str = "pcr_d8"
    coords_names: List[str] = field(default_factory=list)

    @property
    def input_path(self) -> str:
        return os.path.join(INPUT_DIR, self.input_name)

    @property
    def reference_path(self) -> str:
        return os.path.join(REFERENCE_DIR, self.reference_name)


# Register one entry per (input, reference) pair. Extend this list to cover more
# LDD inputs; everything below is data-driven off it.
CASES: List[LddVecCase] = [
    LddVecCase(
        name="pcr_d8_ldd",
        input_name="ldd.nc",
        reference_name="ldd.shp",
        river_network_format="pcr_d8",
    ),
]


def _feature_multiset(gdf: gpd.GeoDataFrame) -> Counter:
    """Reduce a GeoDataFrame to an order-independent multiset of features.

    Each feature is represented by its integer ``ldd`` attribute and its
    geometry normalised to a rounded WKT string, so two shapefiles are equal
    when they contain the same features regardless of row order.
    """
    items = []
    for row in gdf.itertuples(index=False):
        geom = wkt.loads(row.geometry.wkt) if isinstance(row.geometry, str) else row.geometry
        normalised = wkt.dumps(geom, rounding_precision=_COORD_PRECISION)
        items.append((int(row.ldd), normalised))
    return Counter(items)


def _shapefile_components(base_shp: str) -> List[str]:
    """Return the existing sibling components of a ``.shp`` path."""
    base = base_shp[: -len(".shp")] if base_shp.endswith(".shp") else base_shp
    present = []
    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg"):
        candidate = f"{base}{ext}"
        if os.path.exists(candidate):
            present.append(ext)
    return present


@pytest.fixture(params=CASES, ids=lambda c: c.name)
def case(request) -> LddVecCase:
    return request.param


def test_input_and_reference_exist(case: LddVecCase):
    """The declared input netCDF and the core reference components exist."""
    assert os.path.isfile(case.input_path), f"missing input: {case.input_path}"
    for ext in (".shp", ".shx", ".dbf"):
        base = case.reference_path[: -len(".shp")]
        assert os.path.isfile(f"{base}{ext}"), f"missing reference component {ext}"


def test_conversion_matches_reference(case: LddVecCase, tmp_path):
    """Converting the input reproduces the reference shapefile.

    The produced shapefile set must have the same feature count, the same
    (ldd, geometry) multiset, and the same CRS as the reference.
    """
    output_shp = os.path.join(str(tmp_path), "out.shp")

    returned = ldd_to_vector(
        input_file=case.input_path,
        output_file=output_shp,
        river_network_format=case.river_network_format,
        coords_names=list(case.coords_names),
    )

    # The core components must have been written (Req 7.5). A .shp/.shx/.dbf set
    # is always expected.
    produced = _shapefile_components(returned)
    for ext in (".shp", ".shx", ".dbf"):
        assert ext in produced, f"ldd_to_vector did not write {ext}"

    produced_gdf = gpd.read_file(returned)
    reference_gdf = gpd.read_file(case.reference_path)

    # Same number of features.
    assert len(produced_gdf) == len(reference_gdf), (
        f"feature count {len(produced_gdf)} != reference {len(reference_gdf)}"
    )

    # Single integer 'ldd' attribute column plus geometry.
    assert "ldd" in produced_gdf.columns
    assert list(produced_gdf.columns) == list(reference_gdf.columns)

    # Order-independent equality of (ldd, geometry) features.
    assert _feature_multiset(produced_gdf) == _feature_multiset(reference_gdf), (
        "produced features do not match the reference features"
    )

    # The CRS is carried through and matches the reference.
    if reference_gdf.crs is None:
        assert produced_gdf.crs is None
    else:
        assert produced_gdf.crs is not None
        assert produced_gdf.crs.equals(reference_gdf.crs), (
            f"CRS {produced_gdf.crs} != reference {reference_gdf.crs}"
        )


def test_missing_input_raises(case: LddVecCase, tmp_path):
    """A non-existent input path aborts with a ValueError and writes nothing."""
    missing = os.path.join(str(tmp_path), "does_not_exist.nc")
    output_shp = os.path.join(str(tmp_path), "out.shp")

    with pytest.raises(ValueError):
        ldd_to_vector(
            input_file=missing,
            output_file=output_shp,
            river_network_format=case.river_network_format,
            coords_names=list(case.coords_names),
        )

    assert not os.path.exists(output_shp), "no output must be written on failure"


def test_unsupported_format_raises(case: LddVecCase, tmp_path):
    """An unsupported river_network_format aborts with a ValueError."""
    output_shp = os.path.join(str(tmp_path), "out.shp")

    with pytest.raises(ValueError):
        ldd_to_vector(
            input_file=case.input_path,
            output_file=output_shp,
            river_network_format="not_a_format",
            coords_names=list(case.coords_names),
        )
