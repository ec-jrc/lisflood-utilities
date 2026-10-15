"""
Laura Jensen, 2026/06/01
Contact: laura.jensen@gfz.de

This script may be used under CC BY license (https://creativecommons.org/licenses/by/4.0/), which
enables reusers to distribute, remix, adapt, and build upon the material in any medium or format,
so long as attribution is given to the creator.

Creator:
Dr.-Ing. Laura Jensen
laura.jensen@gfz.de
Section 1.3: Earth System Modelling
GFZ Helmholtz Centre for Geosciences

This script can be used for the creation of lakes and reservoirs extent maps in netcdf format,
which are needed for the output of Total Water Storage maps in OS LISFLOOD.
For further documentation view:
https://github.com/ec-jrc/lisflood-code/blob/feature/docs/docs/4_Static-Maps_reservoirs-lakes/index.md

"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import shapely


def create_extent_map(domain, type, is_geographic, file_out, thresh1, thresh2,
                      file_shp1, file_shp2, file_tab, file_loc, file_shp3=None):

    # ******************************************************
    # *** set input file names and read the ID match table
    # ******************************************************
    # --- file_shp1: shapefile with GRanD reservoir outlines (GRAND_ID, AREA_SKM)
    # --- file_shp2: shapefile with HydroLAKES outlines (Hylak_id, Lake_area)
    # --- file_shp3: shapefile with GLWD outlines (GLWD_ID, AREA_SKM); GloFAS domain only
    # --- file_tab:  Excel table with the ID match between LISFLOOD and the source datasets
    # --- file_loc:  netcdf file with outlet locations and IDs of the lakes/reservoirs
    if domain == 'GloFAS':
        if file_shp3 is None:
            raise ValueError('file_shp3 (GLWD shapefile) is required for the GloFAS domain.')

        if type == 'reservoir':
            # --- table with ID match
            tab = pd.read_excel(file_tab, sheet_name="as_reservoirs")
            id_lisf = tab.GDW_ID.values
            id_oth1 = tab.GRAND_ID.values
            id_oth2 = tab.RES_ID.values
            id_oth3 = tab.GLWD_ID.values
            id_src = tab.CATCH_SRC.values
            varname = "Band1"

        elif type == 'lake':
            # --- table with ID match
            tab = pd.read_excel(file_tab)
            id_lisf = tab.FID.values
            id_oth1 = tab.GRAND_ID.values
            id_oth2 = tab.HYLAK_ID.values
            id_oth3 = tab.GLWD_ID.values
            id_src = tab.CATCH_SRC.values
            varname = "Band1"

        else:
            raise ValueError('Invalid type.')

    elif domain == 'ETRS89':
        if type == 'reservoir':
            # --- table with ID match
            tab = pd.read_excel(file_tab, sheet_name="as_reservoirs")
            id_lisf = tab.RES_ID.values
            id_oth1 = tab.GRAND_ID.values
            id_oth2 = tab.HYLAK_ID.values
            id_src = tab.CATCH_SRC.values
            varname = "resnew"

        elif type == 'lake':
            # --- table with ID match
            tab = pd.read_excel(file_tab)
            id_lisf = tab.FID.values
            id_oth1 = tab.GRAND_ID.values
            id_oth2 = tab.HYLAK_ID.values
            id_src = tab.CATCH_SRC.values
            varname = "lakesnew"

        else:
            raise ValueError('Invalid type.')
    else:
        raise ValueError('Invalid type.')

    # ******************************************************
    # *** load res/lakes IDs loc map
    # ******************************************************
    # --- read raster file with locs of lakes/reservoirs ---
    ds = xr.open_dataset(file_loc)
    locM = ds[varname].values.astype(np.int_)

    # --- compute pixel_area
    if is_geographic:
        xname, yname = "lon", "lat"
    else:
        xname, yname = "x", "y"

    x = ds[xname].values
    y = ds[yname].values
    cols, rows = len(x), len(y)

    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])

    if not is_geographic:
        pixel_area = abs(dx) * abs(dy) * np.ones(y.shape)
    else:
        dx_rad = np.deg2rad(dx)
        dy_rad = np.deg2rad(dy)
        lat_rad = np.deg2rad(y)
        pixel_area = 6371000.0**2 * dx_rad * dy_rad * np.cos(lat_rad)

    # ******************************************************
    # *** load shapefiles
    # ******************************************************
    print('Load shapefiles...')

    gdf = gpd.read_file(file_shp1)
    polygons1 = list(gdf.geometry)
    ids1 = list(gdf["GRAND_ID"])
    area1 = list(gdf["AREA_SKM"]*1000*1000)

    gdf = gpd.read_file(file_shp2)
    polygons2 = list(gdf.geometry)
    ids2 = list(gdf["Hylak_id"])
    area2 = list(gdf["Lake_area"]*1000*1000)

    if domain == 'GloFAS':
        gdf = gpd.read_file(file_shp3)
        polygons3 = list(gdf.geometry)
        ids3 = list(gdf["GLWD_ID"])
        area3 = list(gdf["AREA_SKM"]*1000*1000)

    # ******************************************************
    # *** create grid polygons from raster file
    # ******************************************************
    print('Create grid polygons...')

    xx, yy = np.meshgrid(x, y)

    xmin = (xx - dx / 2).ravel()
    ymin = (yy - dy / 2).ravel()
    xmax = (xx + dx / 2).ravel()
    ymax = (yy + dy / 2).ravel()

    cell_polygons = shapely.box(xmin, ymin, xmax, ymax)

    rows_idx, cols_idx = np.indices((rows, cols))
    rows_idx = rows_idx.ravel()
    cols_idx = cols_idx.ravel()
    cell_indices = list(zip(rows_idx, cols_idx))

    grid = gpd.GeoDataFrame(
        {"index": cell_indices},
        geometry=cell_polygons,
        crs=gdf.crs
    )

    # ******************************************************
    # *** compute lake/res extent
    # ******************************************************
    print('Start loop over lakes/reservoirs...')
    out = np.full((rows, cols), 0, dtype=np.int32)  # 0 = NoData

    # loop over each lake/reservoir defined in nc file with outlet locations
    ids_locs = np.unique(locM[locM > 0])
    for i, id_l in enumerate(ids_locs):
        loc = np.where(locM == id_l)
        idx_ref = np.where(id_lisf == id_l)
        if len(idx_ref[0]) == 0:
            print(i, id_l, 'not found in table, lake/res assigned to outlet loc pixel.')
            out[loc[0][0], loc[1][0]] = id_l
        else:
            id_ref = 0
            ids = []
            # if GRAND is the source:
            if id_src[idx_ref] == 'GRAND':
                id_ref = id_oth1[idx_ref]
                ids = ids1
                area = area1
                polygons = polygons1
            # if HYLAK is the source:
            elif id_src[idx_ref] == 'HYLAK':
                id_ref = id_oth2[idx_ref]
                ids = ids2
                area = area2
                polygons = polygons2
            # if GLWD is the source:
            elif (domain == 'GloFAS') and (id_src[idx_ref] == 'GLWD'):
                id_ref = id_oth3[idx_ref]
                ids = ids3
                area = area3
                polygons = polygons3
            else:
                print(i, id_l, "source not available, lake/res assigned to outlet loc pixel.")
                out[loc[0][0], loc[1][0]] = id_l

            if (id_ref > 0) and ~np.isnan(id_ref) and np.any(ids == id_ref):
                idx_shp = np.where(ids == id_ref)[0][0]
                # res is too small
                if (area[idx_shp] / pixel_area[loc[0][0]]) < thresh1:
                    print(i, id_l, 'is too small, lake/res assigned to outlet loc pixel.')
                    out[loc[0][0], loc[1][0]] = id_l
                else:
                    poly = polygons[idx_shp]
                    # use the spatial index to only test the cells that actually
                    # intersect the polygon instead of scanning the whole grid
                    hit_pos = grid.sindex.query(poly, predicate="intersects")
                    count = 0
                    for pos in hit_pos:
                        row, col = cell_indices[pos]
                        cell_geom = cell_polygons[pos]
                        inter = cell_geom.intersection(poly)
                        if not inter.is_empty:
                            frac = inter.area / cell_geom.area
                            # print(id_l, inter.area, cell_geom.area, frac)
                            if frac > thresh2:
                                if out[row, col] == 0:
                                    count += 1
                                    out[row, col] = id_l
                                else:
                                    print(i, id_l, 'overlaps with another lake/res, pixel was not assigned.')

                    if count == 0:
                        out[loc[0][0], loc[1][0]] = id_l
                        print(i, id_l, len(hit_pos), "area fractions too small; lake/res assigned to outlet loc pixel")
                    else:
                        print(i, id_l, count,'of',len(hit_pos), "pixels assigned")
            else:
                print(i, id_l, "POLYGON NOT FOUND, lake/res assigned to outlet loc pixel.")
                out[loc[0][0], loc[1][0]] = id_l

    # ******************************************************
    # *** write out file as copy from input
    # ******************************************************
    out_ds = ds.copy()
    out_ds = out_ds.drop_vars(varname)
    if is_geographic:
        out_ds["polygon_id"] = (("lat", "lon"), out)
    else:
        out_ds["polygon_id"] = (("y", "x"), out)

    out_ds["polygon_id"] = out_ds["polygon_id"].astype("int32")

    encoding = {
        "polygon_id": {
            "zlib": True,
            "complevel": 4,  # 1-9
            "_FillValue": 0,
        }
    }
    out_ds["polygon_id"].attrs.update({
        "standard_name": "polygon_id",
        "long_name": "polygon ids",
        "units": "-",
    })
    if domain == 'ETRS89':
        out_ds["polygon_id"].attrs.update({
            "grid_mapping": "lambert_azimuthal_equal_area"
        })

    out_ds.to_netcdf(file_out, encoding=encoding)

    # close the datasets to release the netCDF file handles
    out_ds.close()
    ds.close()


def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to create lakes and reservoirs extent maps in netcdf format,
        needed for the output of Total Water Storage maps in OS LISFLOOD.
        """,
        prog=prog,
    )
    # --- domain: can be 'GloFAS' or 'ETRS89'
    parser.add_argument("-d", "--domain", required=True, choices=['GloFAS', 'ETRS89'],
                        help="domain: can be 'GloFAS' or 'ETRS89'")
    # --- type: can be 'lake' or 'reservoir'
    parser.add_argument("-t", "--type", required=True, choices=['lake', 'reservoir'],
                        help="type: can be 'lake' or 'reservoir'")
    # --- is_geographic: True if input file has geographic coordinates (in degree), e.g. GloFAS
    #                   False if input file has projected coordinates (in meter), e.g. ETRS89 Use Case
    parser.add_argument("-g", "--is-geographic", action="store_true", default=False,
                        help="set if the input file has geographic coordinates (in degree), e.g. GloFAS; "
                             "otherwise projected coordinates (in meter) are assumed, e.g. ETRS89 Use Case")
    parser.add_argument("-o", "--output", required=True,
                        help="output netcdf file path and name")
    # --- input file paths
    parser.add_argument("--file-shp1", required=True,
                        help="shapefile with GRanD reservoir outlines (GRAND_ID, AREA_SKM)")
    parser.add_argument("--file-shp2", required=True,
                        help="shapefile with HydroLAKES outlines (Hylak_id, Lake_area)")
    parser.add_argument("--file-shp3", required=False, default=None,
                        help="shapefile with GLWD outlines (GLWD_ID, AREA_SKM); required for the GloFAS domain")
    parser.add_argument("--file-tab", required=True,
                        help="Excel table with the ID match between LISFLOOD and the source datasets")
    parser.add_argument("--file-loc", required=True,
                        help="netcdf file with outlet locations and IDs of the lakes/reservoirs")
    # --- thresh1: lakes/res where (total_lake_area / cell_area) < thresh1 are assigned to
    # ---          the pixel of their outlet location as defined in ec_res.nc / ec_lakes.nc;
    # ---          other lakes/res may potentially cover more than one pixel
    # ---          default: thresh1 = 0.05
    parser.add_argument("--thresh1", required=False, type=float, default=0.05,
                        help="lakes/res where (total_lake_area / cell_area) < thresh1 are assigned to "
                             "the pixel of their outlet location; default: 0.05")
    # --- thresh2: percentage of grid cell covered by lake; if exceeded: grid cell is
    # ---          assigned with lake/res ID
    # ---          if thresh2 is exceeded in no pixel: lake/res ID is assigned at outlet
    # ---          location as defined in ec_res.nc / ec_lakes.nc
    # ---          default: thresh2 = 0.07
    parser.add_argument("--thresh2", required=False, type=float, default=0.07,
                        help="percentage of grid cell covered by lake; if exceeded the grid cell is "
                             "assigned with lake/res ID; default: 0.07")

    args = parser.parse_args()

    try:
        create_extent_map(args.domain, args.type, args.is_geographic,
                          args.output, args.thresh1, args.thresh2,
                          file_shp1=args.file_shp1, file_shp2=args.file_shp2,
                          file_tab=args.file_tab, file_loc=args.file_loc,
                          file_shp3=args.file_shp3)
    except Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1)


def main_script():
    # Run the tool and then terminate the process with os._exit to avoid a
    # segmentation fault that the GDAL (geopandas) and HDF5/netCDF libraries can
    # raise while cleaning up at interpreter shutdown. The actual work is already
    # finished at this point, so exiting immediately is safe.
    exit_code = main()
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0 if exit_code is None else exit_code)


if __name__ == '__main__':
    main_script()
