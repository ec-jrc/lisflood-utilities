"""
Helper to (re)generate the small synthetic input files used by tests/test_twsmaps.py.

It creates, in this folder:
    - locations.nc : outlet locations/IDs map (variable "Band1", lat/lon, GloFAS-like)
    - grand.shp    : GRanD reservoir outlines (GRAND_ID, AREA_SKM)
    - hylak.shp    : HydroLAKES outlines (Hylak_id, Lake_area)
    - glwd.shp     : GLWD outlines (GLWD_ID, AREA_SKM)
    - table.xlsx   : ID match table (sheet "as_reservoirs")
and, in the "efas_double_ids" subfolder, the inputs for the EFAS domain case with a
double GRAND_ID in the table (see make_inputs_efas_double_ids).

Run: python tests/data/twsmaps/make_inputs.py
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import box

HERE = os.path.dirname(os.path.abspath(__file__))


def make_inputs(folder=HERE):
    # 4x4 geographic grid, 1-degree cells centred on integer coordinates
    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([3.0, 2.0, 1.0, 0.0])  # descending, as in typical NetCDF maps

    # locations map: reservoir with LISFLOOD id 100, outlet at row 3 (lat=0), col 0 (lon=0)
    loc = np.zeros((4, 4), dtype=np.int32)
    loc[3, 0] = 100
    ds = xr.Dataset(
        {"Band1": (("lat", "lon"), loc)},
        coords={"lat": lat, "lon": lon},
    )
    ds.to_netcdf(os.path.join(folder, "locations.nc"))

    # GRanD reservoir polygon covering the top-left 2x2 block of cells
    # (cells are centred on integers with 1-degree size)
    grand = gpd.GeoDataFrame(
        {"GRAND_ID": [1], "AREA_SKM": [50000.0]},
        geometry=[box(-0.5, 1.5, 1.5, 3.5)],
        crs="EPSG:4326",
    )
    grand.to_file(os.path.join(folder, "grand.shp"))

    # HydroLAKES and GLWD shapefiles are needed by the GloFAS branch even though
    # this test only uses a GRAND source reservoir
    hylak = gpd.GeoDataFrame(
        {"Hylak_id": [10], "Lake_area": [100.0]},
        geometry=[box(2.5, -0.5, 3.5, 0.5)],
        crs="EPSG:4326",
    )
    hylak.to_file(os.path.join(folder, "hylak.shp"))

    glwd = gpd.GeoDataFrame(
        {"GLWD_ID": [20], "AREA_SKM": [100.0]},
        geometry=[box(2.5, 1.5, 3.5, 2.5)],
        crs="EPSG:4326",
    )
    glwd.to_file(os.path.join(folder, "glwd.shp"))

    # ID matching table: LISFLOOD id 100 -> GRAND id 1
    tab = pd.DataFrame({
        "GDW_ID": [100],
        "GRAND_ID": [1],
        "RES_ID": [np.nan],
        "GLWD_ID": [np.nan],
        "CATCH_SRC": ["GRAND"],
        "AREA_SKM": [50000.0],
    })
    with pd.ExcelWriter(os.path.join(folder, "table.xlsx")) as writer:
        tab.to_excel(writer, sheet_name="as_reservoirs", index=False)


def make_inputs_polygon_not_found(folder=None):
    """Inputs for the "polygon not found" fallback test.

    Same as make_inputs, but the locations map contains an extra reservoir (LISFLOOD
    id 200) that is present in the ID table (source GRAND, GRAND_ID 99) but whose
    polygon does NOT exist in the GRanD shapefile. With the fix in twsmaps, this
    reservoir must still be assigned to its outlet location pixel.
    """
    if folder is None:
        folder = os.path.join(HERE, "polygon_not_found")
    os.makedirs(folder, exist_ok=True)

    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([3.0, 2.0, 1.0, 0.0])

    # reservoir 100 (has a polygon) at outlet (row 3, col 0);
    # reservoir 200 (polygon NOT in shapefile) at outlet (row 3, col 3)
    loc = np.zeros((4, 4), dtype=np.int32)
    loc[3, 0] = 100
    loc[3, 3] = 200
    ds = xr.Dataset(
        {"Band1": (("lat", "lon"), loc)},
        coords={"lat": lat, "lon": lon},
    )
    ds.to_netcdf(os.path.join(folder, "locations.nc"))

    # GRanD shapefile only contains the polygon for reservoir 100 (GRAND_ID 1)
    grand = gpd.GeoDataFrame(
        {"GRAND_ID": [1], "AREA_SKM": [50000.0]},
        geometry=[box(-0.5, 1.5, 1.5, 3.5)],
        crs="EPSG:4326",
    )
    grand.to_file(os.path.join(folder, "grand.shp"))

    hylak = gpd.GeoDataFrame(
        {"Hylak_id": [10], "Lake_area": [100.0]},
        geometry=[box(2.5, -0.5, 3.5, 0.5)],
        crs="EPSG:4326",
    )
    hylak.to_file(os.path.join(folder, "hylak.shp"))

    glwd = gpd.GeoDataFrame(
        {"GLWD_ID": [20], "AREA_SKM": [100.0]},
        geometry=[box(2.5, 1.5, 3.5, 2.5)],
        crs="EPSG:4326",
    )
    glwd.to_file(os.path.join(folder, "glwd.shp"))

    # table: id 100 -> GRAND_ID 1 (exists); id 200 -> GRAND_ID 99 (missing in shapefile)
    tab = pd.DataFrame({
        "GDW_ID": [100, 200],
        "GRAND_ID": [1, 99],
        "RES_ID": [np.nan, np.nan],
        "GLWD_ID": [np.nan, np.nan],
        "CATCH_SRC": ["GRAND", "GRAND"],
        "AREA_SKM": [50000.0, 50000.0],
    })
    with pd.ExcelWriter(os.path.join(folder, "table.xlsx")) as writer:
        tab.to_excel(writer, sheet_name="as_reservoirs", index=False)


def make_inputs_efas_double_ids(folder=None):
    """Inputs for the EFAS domain case with a double source ID in the table.

    Reservoirs 300 and 301 share the same GRAND_ID (5) in the ID table; the one with
    the largest AREA_SKM (301) must be kept and the other one (300) removed from the
    locations map by the delete_double workaround. Reservoir 302 is sourced from GDW
    (EFAS third source). The locations map uses the EFAS variable name "reservoirs".
    """
    if folder is None:
        folder = os.path.join(HERE, "efas_double_ids")
    os.makedirs(folder, exist_ok=True)

    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([3.0, 2.0, 1.0, 0.0])

    # reservoir 300 (double GRAND_ID, small area -> removed) at outlet (row 3, col 0);
    # reservoir 301 (double GRAND_ID, largest area -> kept) at outlet (row 3, col 1);
    # reservoir 302 (GDW source) at outlet (row 0, col 0)
    loc = np.zeros((4, 4), dtype=np.int32)
    loc[3, 0] = 300
    loc[3, 1] = 301
    loc[0, 0] = 302
    ds = xr.Dataset(
        {"reservoirs": (("lat", "lon"), loc)},
        coords={"lat": lat, "lon": lon},
    )
    ds.to_netcdf(os.path.join(folder, "locations.nc"))

    # GRanD polygon (GRAND_ID 5, shared by 300 and 301) covering the cells
    # (rows 2-3) x (cols 1-3); the outlet cell of 300 (row 3, col 0) is not covered
    grand = gpd.GeoDataFrame(
        {"GRAND_ID": [5], "AREA_SKM": [50000.0]},
        geometry=[box(0.5, -0.5, 3.5, 1.5)],
        crs="EPSG:4326",
    )
    grand.to_file(os.path.join(folder, "grand.shp"))

    # HydroLAKES shapefile is loaded by the EFAS branch even though no reservoir
    # in the table uses the HYLAK source
    hylak = gpd.GeoDataFrame(
        {"Hylak_id": [10], "Lake_area": [100.0]},
        geometry=[box(-0.5, -0.5, 0.5, 0.5)],
        crs="EPSG:4326",
    )
    hylak.to_file(os.path.join(folder, "hylak.shp"))

    # GDW polygon for reservoir 302 (in EFAS the LISFLOOD ids are the GDW ids)
    gdw = gpd.GeoDataFrame(
        {"GDW_ID": [302], "AREA_SKM": [50000.0]},
        geometry=[box(-0.5, 2.5, 0.5, 3.5)],
        crs="EPSG:4326",
    )
    gdw.to_file(os.path.join(folder, "gdw.shp"))

    # table: 300 and 301 share GRAND_ID 5 with different areas; 302 is a GDW reservoir
    # (in the EFAS table the LISFLOOD ids are the GDW ids, as in the GFZ script)
    tab = pd.DataFrame({
        "GDW_ID": [300, 301, 302],
        "GRAND_ID": [5, 5, np.nan],
        "HYLAK_ID": [np.nan, np.nan, np.nan],
        "CATCH_SRC": ["GRAND", "GRAND", "GDW"],
        "AREA_SKM": [100.0, 50000.0, 50000.0],
    })
    with pd.ExcelWriter(os.path.join(folder, "table.xlsx")) as writer:
        tab.to_excel(writer, sheet_name="as_reservoirs", index=False)


if __name__ == "__main__":
    make_inputs()
    make_inputs_polygon_not_found()
    make_inputs_efas_double_ids()
    print("Test inputs generated in", HERE)
