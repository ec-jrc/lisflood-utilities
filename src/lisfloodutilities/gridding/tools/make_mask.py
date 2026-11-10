#!/usr/bin/env python
"""Create a NetCDF boolean mask map from a GeoTIFF.

Valid (non-nodata) cells in the input GeoTIFF become mask value 1;
all other cells become 0. Georeferencing (lat/lon, CRS) is preserved.
"""
import argparse

import numpy as np
import rasterio
from netCDF4 import Dataset


def make_mask(src_tif: str, dst_nc: str, crop: bool = False) -> None:
    with rasterio.open(src_tif) as src:
        data = src.read(1)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs

    # Valid cells = cells that are not nodata (and finite)
    if nodata is not None:
        valid = data != nodata
    else:
        valid = np.ones(data.shape, dtype=bool)
    if np.issubdtype(data.dtype, np.floating):
        valid &= np.isfinite(data)

    rows, cols = data.shape

    # Optionally crop to the bounding box of valid cells
    if crop and valid.any():
        ridx = np.where(valid.any(axis=1))[0]
        cidx = np.where(valid.any(axis=0))[0]
        r0, r1 = ridx[0], ridx[-1] + 1
        c0, c1 = cidx[0], cidx[-1] + 1
        valid = valid[r0:r1, c0:c1]
        # shift the affine origin to the cropped window
        transform = transform * rasterio.Affine.translation(c0, r0)
        rows, cols = valid.shape

    mask = valid.astype(np.int8)

    # Build coordinate axes from the affine transform (cell centres)
    a, b, c, d, e, f = (transform.a, transform.b, transform.c,
                        transform.d, transform.e, transform.f)
    lon = c + a * (np.arange(cols) + 0.5)
    lat = f + e * (np.arange(rows) + 0.5)

    with Dataset(dst_nc, "w", format="NETCDF4") as nc:
        nc.createDimension("lat", rows)
        nc.createDimension("lon", cols)

        vlat = nc.createVariable("lat", "f8", ("lat",))
        vlat.standard_name = "latitude"
        vlat.long_name = "latitude"
        vlat.units = "degrees_north"
        vlat[:] = lat

        vlon = nc.createVariable("lon", "f8", ("lon",))
        vlon.standard_name = "longitude"
        vlon.long_name = "longitude"
        vlon.units = "degrees_east"
        vlon[:] = lon

        vmask = nc.createVariable(
            "mask", "i1", ("lat", "lon"),
            zlib=True, complevel=4, fill_value=0,
        )
        vmask.long_name = "boolean mask (1 = valid cell)"
        vmask.grid_mapping = "crs"
        vmask[:] = mask

        vcrs = nc.createVariable("crs", "i4")
        if crs is not None:
            try:
                vcrs.spatial_ref = crs.to_wkt()
                epsg = crs.to_epsg()
                if epsg is not None:
                    vcrs.epsg_code = f"EPSG:{epsg}"
            except Exception:
                pass

        nc.Conventions = "CF-1.8"
        nc.source = src_tif

    print(f"Wrote {dst_nc}: {rows} x {cols}, {int(mask.sum())} valid cells")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("src", help="input GeoTIFF")
    p.add_argument("dst", help="output NetCDF")
    p.add_argument("--crop", action="store_true",
                help="crop to bounding box of valid cells")
    args = p.parse_args()
    make_mask(args.src, args.dst, crop=args.crop)
