"""
Copyright 2019-2026 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

Adjust gridded annual maximum discharge using a linear sqrt(R2)-weighted
trend correction.

Original script: adjust_annual_maxima_trends_linear_r2_oper.py
Author: Ervin Zsoter (ECMWF)
Adapted as a lisflood-utilities tool by the EC-JRC.

For each grid cell with a complete valid annual-maxima series, fit an
ordinary least-squares linear trend and shift each value towards the fitted
trend level at the chosen reference year. Weight the correction by sqrt(R2)
and cap adjusted values at zero:

    x = year - mean(years)
    T(year) = slope*x + intercept
    R2 = 1 - SSE/SST
    weight = sqrt(clip(R2, 0, 1))
    Qadj(year) = max(0, Q(year) + weight*(T(reference_year) - T(year)))

Constant series receive zero weight. Input values are otherwise unchanged
apart from the weighted correction and non-negative cap.

Input requirements:
  - A three-dimensional variable ordered (time, latitude, longitude), with
    a coordinate variable for each dimension.
  - At least two time steps, each belonging to a distinct calendar year.
  - A CF time coordinate with units; calendar defaults to "standard".
  - The reference year must occur in the input time coordinate.

Years are taken directly from decoded timestamps. No adjustment is made
for timestamps representing the end of a preceding annual period.

A grid cell is processed only if every annual value is finite and has
absolute magnitude below 1e19. Masked input values are treated as missing.
If any annual value is invalid, the entire output series for that cell is
missing. Missing years are not fitted or reconstructed.

Processing reads and writes latitude chunks. The default chunk size is
25 rows; the input variable defaults to "var7" and the reference year to
2025. POINT station inputs are not supported.

Output:
  - Compressed float32 variable "annual_maximum_adjusted", with fill value
    1e20 and units copied from the input variable when available.
  - The three input coordinate variables and their attributes.
  - Global metadata describing the source, years, reference year, method,
    valid fitted pixel count and number of values capped at zero.

Other input variables and original global attributes are not copied.

With --diagnostics, also write:
  linear_trend                        fitted trend at each input year
  linear_slope                        slope per year
  linear_intercept                    trend value at mean(years)
  linear_r2                           R2 before clipping for the weight
  linear_sqrt_r2_weight               sqrt(clip(R2, 0, 1))
  linear_trend_reference_value        fitted trend at the reference year
  linear_first_year_weighted_adjustment
                                      correction at the first input time step,
                                      before applying the non-negative cap

Usage:
  adjust_trends INPUT OUTPUT --variable var7 --reference-year 2025 \
      --latitude-chunk 25

Add --diagnostics to include diagnostic fields.
Add --overwrite to replace an existing output file. An existing output is
deleted before input validation; input and output paths must be different.
"""

import argparse
import datetime as dt
import os
import sys
import time

import numpy as np
from netCDF4 import Dataset, num2date

FILL = np.float32(1.0e20)
MISSING_LIMIT = 1.0e19


def parse_args(argv=None):
    """Parse the command-line arguments of the adjust_trends tool."""
    parser = argparse.ArgumentParser(
        description="Linear sqrt(R2)-weighted adjustment of annual maxima.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=os.path.basename(argv[0]) if argv else 'adjust_trends',
    )
    parser.add_argument("input_file", help="Input NetCDF file with annual maxima")
    parser.add_argument("output_file", help="Output NetCDF file with adjusted annual maxima")
    parser.add_argument("--variable", default="var7", help="Annual-maximum discharge variable")
    parser.add_argument("--reference-year", type=int, default=2025,
                        help="Trend level to which all annual maxima are adjusted")
    parser.add_argument("--latitude-chunk", type=int, default=25,
                        help="Number of latitude rows processed at once")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Also write fitted linear trend, slope, intercept, R2, sqrt(R2) weight "
                             "and adjustment diagnostics")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite the output file if it exists")
    return parser.parse_args(argv[1:] if argv is not None else None)


def read_years(tvar):
    """Return the calendar years of a NetCDF time coordinate variable."""
    if not hasattr(tvar, "units"):
        raise ValueError(f"Time variable {tvar.name!r} has no units attribute")
    dates = num2date(tvar[:], units=tvar.units, calendar=getattr(tvar, "calendar", "standard"))
    return np.asarray([d.year for d in dates], dtype=np.int32)


def copy_variable(src, dst, name):
    """Copy a coordinate variable and its attributes from the source to the destination dataset."""
    v = src.variables[name]
    fill = getattr(v, "_FillValue", None)
    out = dst.createVariable(name, v.dtype, v.dimensions,
                             **({"fill_value": fill} if fill is not None else {}))
    out[:] = v[:]
    for att in v.ncattrs():
        if att != "_FillValue":
            out.setncattr(att, v.getncattr(att))


def create_3d(nc, name, dims, lat_chunk, ny, nx):
    """Create a compressed 3-D float32 output variable."""
    return nc.createVariable(name, "f4", dims, fill_value=FILL, zlib=True, complevel=4, shuffle=True,
                             chunksizes=(1, min(lat_chunk, ny), min(720, nx)))


def create_2d(nc, name, dims):
    """Create a compressed 2-D float32 output variable."""
    return nc.createVariable(name, "f4", dims, fill_value=FILL, zlib=True, complevel=4, shuffle=True)


def valid_mask(data):
    """Return the mask of cells whose whole time series is finite and below the missing limit."""
    return np.all(np.isfinite(data), axis=0) & np.all(np.abs(data) < MISSING_LIMIT, axis=0)


def fit_r2(y, fitted):
    """Compute R2 and the sqrt(clip(R2,0,1)) weight for each column of y."""
    residual = y - fitted
    sse = np.sum(residual * residual, axis=0, dtype=np.float64)
    mean = np.mean(y, axis=0, dtype=np.float64)
    sst = np.sum((y - mean[None, :]) ** 2, axis=0, dtype=np.float64)
    r2 = np.zeros(y.shape[1], dtype=np.float64)
    ok = sst > 0.0
    r2[ok] = 1.0 - sse[ok] / sst[ok]
    return r2, np.sqrt(np.clip(r2, 0.0, 1.0))


def process(a):
    """
    Adjust gridded annual maximum discharge using a linear sqrt(R2)-weighted trend correction.

    For each grid cell with a complete valid annual-maxima series, fit an ordinary least-squares
    linear trend and shift each value towards the fitted trend level at the chosen reference year,
    weighted by sqrt(R2) and capped at zero:

        x = year - mean(years)
        T(year) = slope*x + intercept
        R2 = 1 - SSE/SST
        weight = sqrt(clip(R2, 0, 1))
        Qadj(year) = max(0, Q(year) + weight*(T(reference_year) - T(year)))

    :param a: parsed command-line arguments (see parse_args)
    """
    if a.latitude_chunk < 1:
        raise ValueError("--latitude-chunk must be at least 1")
    if os.path.exists(a.output_file):
        if not a.overwrite:
            raise FileExistsError(f"Output exists: {a.output_file!r}; use --overwrite")
        os.remove(a.output_file)

    start = time.perf_counter()
    with Dataset(a.input_file) as src:
        if a.variable not in src.variables:
            raise KeyError(f"Variable {a.variable!r} not found")
        qvar = src.variables[a.variable]
        if qvar.ndim != 3:
            raise ValueError(f"{a.variable!r} must be 3-D; found {qvar.dimensions}")

        tname, yname, xname = qvar.dimensions
        for name in (tname, yname, xname):
            if name not in src.variables:
                raise KeyError(f"Coordinate variable {name!r} not found")

        years = read_years(src.variables[tname])
        nt, ny, nx = len(years), len(src.dimensions[yname]), len(src.dimensions[xname])
        if nt < 2:
            raise ValueError("At least two annual maxima are required")
        if np.unique(years).size != nt:
            raise ValueError("More than one time step found for at least one year")
        if a.reference_year not in years:
            raise ValueError(f"Reference year {a.reference_year} not present in input")

        centre_year = float(np.mean(years))
        x = years.astype(np.float64) - centre_year
        xref = float(a.reference_year) - centre_year
        design = np.column_stack((x, np.ones(nt)))
        inverse = np.linalg.pinv(design)
        reference = np.array([xref, 1.0])
        nchunk = (ny + a.latitude_chunk - 1) // a.latitude_chunk

        print("=" * 72)
        print("LINEAR sqrt(R2) ANNUAL-MAXIMUM ADJUSTMENT")
        print("=" * 72)
        print(f"Input              : {a.input_file}")
        print(f"Variable           : {a.variable}")
        print(f"Output             : {a.output_file}")
        print(f"Grid               : {ny} x {nx}")
        print(f"Years              : {years.min()}-{years.max()} ({nt})")
        print(f"Reference year     : {a.reference_year}")
        print("Weight             : sqrt(clip(R2,0,1))")
        print(f"Diagnostics        : {a.diagnostics}")
        print(f"Latitude chunk     : {a.latitude_chunk}")
        print("=" * 72)

        with Dataset(a.output_file, "w", format="NETCDF4") as dst:
            dst.createDimension(tname, nt)
            dst.createDimension(yname, ny)
            dst.createDimension(xname, nx)
            for name in (tname, yname, xname):
                copy_variable(src, dst, name)

            dims3 = (tname, yname, xname)
            dims2 = (yname, xname)
            adjusted = create_3d(dst, "annual_maximum_adjusted", dims3, a.latitude_chunk, ny, nx)
            adjusted.long_name = "annual maximum discharge adjusted using sqrt(R2)-weighted pixel-wise linear trend"
            adjusted.reference_year = a.reference_year
            adjusted.adjustment_weight = "sqrt(clip(R2,0,1))"
            if hasattr(qvar, "units"):
                adjusted.units = qvar.units

            if a.diagnostics:
                trend_var = create_3d(dst, "linear_trend", dims3, a.latitude_chunk, ny, nx)
                slope_var = create_2d(dst, "linear_slope", dims2)
                intercept_var = create_2d(dst, "linear_intercept", dims2)
                r2_var = create_2d(dst, "linear_r2", dims2)
                weight_var = create_2d(dst, "linear_sqrt_r2_weight", dims2)
                ref_var = create_2d(dst, "linear_trend_reference_value", dims2)
                first_adj_var = create_2d(dst, "linear_first_year_weighted_adjustment", dims2)
                r2_var.units = weight_var.units = "1"
                if hasattr(qvar, "units"):
                    trend_var.units = intercept_var.units = ref_var.units = first_adj_var.units = qvar.units
                    slope_var.units = f"{qvar.units} year-1"

            dst.title = "Linear sqrt(R2)-weighted trend-adjusted annual maxima"
            dst.source_file = os.path.abspath(a.input_file)
            dst.source_variable = a.variable
            dst.year_range = f"{years.min()}-{years.max()}"
            dst.number_of_years = nt
            dst.reference_year = a.reference_year
            dst.trend_centre_year = centre_year
            dst.trend_model = "linear OLS"
            dst.adjustment_type = "additive"
            dst.adjustment_formula = "Q_adjusted=max(0,Q+sqrt(clip(R2,0,1))*(T(reference_year)-T(year)))"
            dst.history = f"{dt.datetime.now(dt.timezone.utc).isoformat()} created by {os.path.basename(__file__)}"

            valid_total = 0
            zero_caps = 0

            for ichunk, j0 in enumerate(range(0, ny, a.latitude_chunk), start=1):
                j1 = min(j0 + a.latitude_chunk, ny)
                clock = time.perf_counter()
                chunk_ny = j1 - j0
                npixel = chunk_ny * nx

                raw = qvar[:, j0:j1, :]
                data = raw.filled(np.nan) if np.ma.isMaskedArray(raw) else np.asarray(raw)
                data = np.asarray(data, dtype=np.float64).reshape(nt, npixel)
                valid = valid_mask(data)
                nvalid = np.count_nonzero(valid)
                valid_total += nvalid

                out = np.full((nt, npixel), FILL, dtype=np.float32)

                if a.diagnostics:
                    trend_out = np.full((nt, npixel), FILL, dtype=np.float32)
                    diag_out = [np.full(npixel, FILL, dtype=np.float32) for _ in range(6)]

                if nvalid:
                    y = data[:, valid]
                    coef = inverse @ y
                    fitted = design @ coef
                    ref_value = reference @ coef
                    r2, weight = fit_r2(y, fitted)
                    yadj = y + weight[None, :] * (ref_value[None, :] - fitted)
                    zero_caps += int(np.count_nonzero(yadj < 0.0))
                    yadj = np.maximum(yadj, 0.0)
                    out[:, valid] = yadj.astype(np.float32)

                    if a.diagnostics:
                        trend_out[:, valid] = fitted.astype(np.float32)
                        diag_out[0][valid] = coef[0].astype(np.float32)
                        diag_out[1][valid] = coef[1].astype(np.float32)
                        diag_out[2][valid] = r2.astype(np.float32)
                        diag_out[3][valid] = weight.astype(np.float32)
                        diag_out[4][valid] = ref_value.astype(np.float32)
                        diag_out[5][valid] = (weight * (ref_value - fitted[0])).astype(np.float32)

                adjusted[:, j0:j1, :] = out.reshape(nt, chunk_ny, nx)

                if a.diagnostics:
                    trend_var[:, j0:j1, :] = trend_out.reshape(nt, chunk_ny, nx)
                    for var, values in zip((slope_var, intercept_var, r2_var, weight_var, ref_var, first_adj_var), diag_out):
                        var[j0:j1, :] = values.reshape(chunk_ny, nx)

                print(f"Chunk {ichunk:4d}/{nchunk:4d} lat {j0:4d}:{j1:4d}  valid={nvalid:,}  time={time.perf_counter()-clock:.1f}s")

            dst.valid_pixel_count = int(valid_total)
            dst.values_capped_at_zero = int(zero_caps)

    print("=" * 72)
    print(f"Valid fitted pixels : {valid_total:,}")
    print(f"Values capped at 0  : {zero_caps:,}")
    print(f"Output              : {a.output_file}")
    print(f"Elapsed time        : {time.perf_counter()-start:.1f}s")
    print("=" * 72)


def main(argv=sys.argv):
    """Command-line entry point of the adjust_trends tool."""
    args = parse_args(argv)
    try:
        process(args)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


def main_script():
    sys.exit(main())


if __name__ == "__main__":
    main_script()
