#   ldd2vec converts a local drain direction (LDD) raster into a vector line
#   network (ESRI Shapefile). Each flow cell becomes a line segment pointing to
#   its downstream neighbour, and each pit/sink cell becomes a diagonal cross.
#
#   Usage:
#   ldd2vec -i ldd.nc -o ldd.shp -f pcr_d8 -E y x

import sys
import os
import tempfile
from typing import Dict, List, Optional, Tuple
import dataclasses
import numpy as np
import xarray as xr
import geopandas
from shapely.geometry import LineString
import pyproj
from lisfloodutilities.mctrivers.helpers import (
    X_COORD_NAMES,
    Y_COORD_NAMES,
)

# River network format direction tables.
#
# Each supported river network format defines the set of valid flow-direction
# LDD values (mapped to a compass direction), the pit/sink value, and the
# no-data/fill value. The compass directions are converted to array-index
# offsets later using the sign of each coordinate axis step (see
# _direction_offset).
#
# WARNING: These mappings MUST be verified against the earthkit.hydro,
# PCRaster, and ESRI documentation before the tool is considered correct.
@dataclasses.dataclass(frozen=True)
class FormatSpec:
    """Encoding of a river network format.

    Attributes:
        flow: mapping of flow LDD value -> compass direction (e.g. "NE").
        pit: the LDD value denoting a pit/sink cell (no downstream neighbour).
        fill: the LDD value denoting a no-data / fill cell.
    """
    flow: Dict[int, str]
    pit: int
    fill: int


FORMAT_SPECS: Dict[str, FormatSpec] = {
    "pcr_d8": FormatSpec(
        flow={1: "SW", 2: "S", 3: "SE", 4: "W", 6: "E", 7: "NW", 8: "N", 9: "NE"},
        pit=5,
        fill=255,
    ),
    "esri_d8": FormatSpec(
        flow={1: "E", 2: "SE", 4: "S", 8: "SW", 16: "W", 32: "NW", 64: "N", 128: "NE"},
        pit=0,
        fill=255,
    ),
    "merit_d8": FormatSpec(
        flow={1: "E", 2: "SE", 4: "S", 8: "SW", 16: "W", 32: "NW", 64: "N", 128: "NE"},
        pit=0,
        fill=247,
    ),
}


def _validate_river_network_format(river_network_format: str) -> None:
    """Validate a candidate river_network_format against FORMAT_SPECS.

    Defence in depth behind the argparse ``choices`` restriction: leaves the
    input unchanged and raises ``ValueError`` naming the unsupported value when
    it is not a key of ``FORMAT_SPECS`` (Req 2.3).

    :param river_network_format: candidate river network format name.
    :raises ValueError: when ``river_network_format`` is not a supported format,
        naming the unsupported value and listing the supported formats.
    """
    if river_network_format not in FORMAT_SPECS:
        supported = ", ".join(sorted(FORMAT_SPECS.keys()))
        raise ValueError(
            f"Unsupported river_network_format '{river_network_format}'. "
            f"Supported formats: {supported}."
        )


def extract_coords(ds: xr.Dataset, coords_names: List[str] = []) -> Tuple[str, str]:
    """Resolve the x/longitude and y/latitude coordinate variable names.

    Adapted from ``lddfix.extract_coords``. When ``coords_names`` is empty the
    names are auto-detected by intersecting the dataset coordinates with the
    package naming conventions ``X_COORD_NAMES``/``Y_COORD_NAMES``, requiring
    exactly one match per axis (Req 6.2). When a two-element ``[y, x]`` list is
    supplied those names are used directly, after verifying both are present in
    the file (Req 6.3). Any other number of names is rejected (Req 6.4).

    :param ds: input dataset whose coordinate names are resolved.
    :param coords_names: optional explicit ``[y, x]`` coordinate names.
    :return: ``(x_proj, y_proj)`` coordinate variable names.
    :raises ValueError: when auto-detection is ambiguous (zero or multiple
        matches on an axis, Req 6.6), when a supplied name is absent from the
        file (Req 6.5), or when ``coords_names`` does not consist of exactly two
        names (Req 6.4). The message requests explicit ``-E`` names.
    """
    x_proj, y_proj = '', ''
    if len(coords_names) == 0:
        # Auto-detect coordinates by matching against known patterns (Req 6.2).
        available_coords = set(ds.coords)
        x_matches = available_coords.intersection(X_COORD_NAMES)
        y_matches = available_coords.intersection(Y_COORD_NAMES)

        # A unique match is required on each axis; zero or multiple candidates
        # mean the names cannot be determined, so request explicit -E (Req 6.6).
        if len(x_matches) != 1 or len(y_matches) != 1:
            raise ValueError(
                f"Unable to uniquely identify coordinates. "
                f"Found {len(x_matches)} x-coords: {x_matches}, "
                f"{len(y_matches)} y-coords: {y_matches}. "
                f"Available: {list(available_coords)}. "
                f"Please specify coordinates explicitly using the -E argument "
                f"(format: 'lat lon')."
            )

        x_proj = str(next(iter(x_matches)))
        y_proj = str(next(iter(y_matches)))
    elif len(coords_names) == 2:
        # Use provided coordinates (y, x order as per documentation) (Req 6.3).
        y_proj, x_proj = coords_names
        # A supplied name that is absent from the file is rejected (Req 6.5).
        missing = [name for name in (y_proj, x_proj) if name not in ds.variables]
        if missing:
            raise ValueError(
                f"Supplied coordinate name(s) not found in the input file: "
                f"{missing}. Available variables: {list(ds.variables)}. "
                f"Please specify coordinates explicitly using the -E argument "
                f"(format: 'lat lon')."
            )
    else:
        # Exactly two names are required when an explicit list is given (Req 6.4).
        raise ValueError(
            f"coords_names must have 0 or 2 elements, got {len(coords_names)}: "
            f"{coords_names}. Please specify coordinates explicitly using the "
            f"-E argument (format: 'lat lon')."
        )

    return x_proj, y_proj


def _resolve_ldd_var_name(ds: xr.Dataset, x_proj: str, y_proj: str) -> str:
    """Resolve the name of the LDD data variable in the input dataset.

    Adapted from ``lddfix._resolve_ldd_var_name``. The input is a 2-D LDD raster
    carrying a single direction value per cell, so the LDD variable is the lone
    data variable spanning both the resolved x and y coordinate axes. Coordinate
    variables (``x_proj``/``y_proj``) are excluded.

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :return: name of the LDD data variable.
    :raises ValueError: when zero or more than one candidate data variable spans
        both coordinate axes, identifying the file as not containing a usable
        LDD raster (Req 1.5).
    """
    candidates = [
        str(name)
        for name, var in ds.data_vars.items()
        if name not in (x_proj, y_proj)
        and x_proj in var.dims
        and y_proj in var.dims
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"The input file does not contain a usable LDD raster. Found "
            f"{len(candidates)} candidate data variable(s) spanning the "
            f"'{y_proj}'/'{x_proj}' axes: {candidates}. A usable LDD raster "
            f"requires exactly one 2-D data variable on those axes. "
            f"Available data variables: {list(ds.data_vars)}."
        )
    return candidates[0]


def _axis_steps(
    ds: xr.Dataset, x_proj: str, y_proj: str
) -> Tuple[float, float, float, float]:
    """Signed per-axis step and half-cell sizes of the input grid.

    The *signed* step of each axis is the spacing between the first two
    coordinate values (``x[1] - x[0]`` and ``y[1] - y[0]``). The sign is kept
    deliberately: it encodes the axis orientation (ascending vs descending) that
    ``_direction_offset`` uses to turn a compass direction into an array-index
    offset, so that a descending y axis (the common netCDF layout) maps "north"
    to a decreasing row index (Req 3.3). This differs from
    ``lddfix._axis_half_resolution``, which discards the sign because it only
    needs a selection tolerance.

    The half-cell sizes are the absolute half-spacings ``half_w = |x_step| / 2``
    and ``half_h = |y_step| / 2``; they locate the four cell corners half a cell
    width in x and half a cell height in y from the cell centre (Req 3.5) and are
    therefore orientation-independent.

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :return: ``(x_step, y_step, half_w, half_h)`` where ``x_step``/``y_step``
        are signed and ``half_w``/``half_h`` are non-negative.
    :raises ValueError: when either axis has fewer than two coordinate values,
        so a step cannot be determined.
    """
    x_coords = np.asarray(ds[x_proj].values)
    y_coords = np.asarray(ds[y_proj].values)

    if x_coords.size < 2 or y_coords.size < 2:
        raise ValueError(
            "Cannot determine axis steps: each spatial axis requires at least "
            f"two coordinate values, but the x axis '{x_proj}' has "
            f"{x_coords.size} and the y axis '{y_proj}' has {y_coords.size}."
        )

    x_step = float(x_coords[1]) - float(x_coords[0])
    y_step = float(y_coords[1]) - float(y_coords[0])
    half_w = abs(x_step) / 2.0
    half_h = abs(y_step) / 2.0

    return x_step, y_step, half_w, half_h


def _resolve_crs(
    ds: xr.Dataset, ldd_var_name: str
) -> Tuple[Optional[pyproj.CRS], bool]:
    """Resolve the input CRS from the LDD data variable's attributes.

    The CRS attributes are read from ``ds[ldd_var_name].attrs`` in priority
    order -- ``spatial_ref``, ``crs_wkt``, ``esri_pe_string``, ``grid_mapping``
    -- following the attribute set carried by
    ``lisfloodutilities.mctrivers.helpers.copy_coordinates_and_attributes``. For
    ``spatial_ref``/``crs_wkt``/``esri_pe_string`` the attribute value is a
    WKT/PE string interpreted directly via ``pyproj.CRS.from_user_input``. For
    ``grid_mapping`` the attribute names another variable whose own attributes
    carry the CRS (``crs_wkt``/``spatial_ref``/``esri_pe_string``), so that
    variable's attributes are looked up and interpreted.

    :param ds: input dataset carrying the LDD raster and any CRS metadata.
    :param ldd_var_name: name of the LDD data variable whose attributes hold the
        CRS metadata.
    :return: ``(crs, warn)`` where

        * ``(pyproj.CRS, False)`` when a CRS attribute is present and
          interpretable (Req 5.1);
        * ``(None, False)`` when no CRS attribute is present at all -- the
          output is written with an unset CRS (Req 5.2);
        * ``(None, True)`` when a CRS attribute is present but none of the
          candidates can be interpreted; the warn flag signals the caller to
          emit a warning while still writing a complete output with unset CRS
          and exit zero (Req 5.3).
    """
    attrs = ds[ldd_var_name].attrs

    # Collect candidate CRS source strings in priority order. A candidate is
    # recorded whenever its attribute is present (so we can distinguish "no CRS
    # attribute present" from "present but uninterpretable"), regardless of
    # whether it ultimately interprets.
    candidates: List[str] = []
    crs_present = False

    # spatial_ref / crs_wkt / esri_pe_string hold a WKT/PE string directly.
    for attr_name in ("spatial_ref", "crs_wkt", "esri_pe_string"):
        if attr_name in attrs:
            crs_present = True
            value = attrs[attr_name]
            if value is not None and str(value).strip():
                candidates.append(str(value))

    # grid_mapping names another variable whose attributes carry the CRS.
    if "grid_mapping" in attrs:
        crs_present = True
        grid_mapping_name = attrs["grid_mapping"]
        if grid_mapping_name is not None and str(grid_mapping_name) in ds.variables:
            gm_attrs = ds[str(grid_mapping_name)].attrs
            for attr_name in ("crs_wkt", "spatial_ref", "esri_pe_string"):
                if attr_name in gm_attrs:
                    value = gm_attrs[attr_name]
                    if value is not None and str(value).strip():
                        candidates.append(str(value))

    # No CRS attribute present at all: unset CRS, no warning (Req 5.2).
    if not crs_present:
        return None, False

    # Try each candidate in priority order; the first interpretable one wins.
    for candidate in candidates:
        try:
            return pyproj.CRS.from_user_input(candidate), False
        except Exception:
            continue

    # A CRS attribute was present but none could be interpreted (Req 5.3).
    return None, True


# Compass direction -> spatial unit offset (Delta-lon, Delta-lat).
#
# Defined in spatial terms first (east is increasing longitude, north is
# increasing latitude), independent of array-index orientation. The per-axis
# sign of the coordinate step converts these to array-index offsets in
# _direction_offset. See design Data Models (Compass direction -> array-index
# offset).
COMPASS_OFFSETS: Dict[str, Tuple[int, int]] = {
    "E": (1, 0),
    "W": (-1, 0),
    "N": (0, 1),
    "S": (0, -1),
    "NE": (1, 1),
    "NW": (-1, 1),
    "SE": (1, -1),
    "SW": (-1, -1),
}


def _direction_offset(
    value: int, fmt: FormatSpec, x_step: float, y_step: float
) -> Optional[Tuple[int, int]]:
    """Map a flow LDD value to an array-index offset ``(drow, dcol)``.

    The value is looked up in the format's compass table to obtain a spatial
    unit offset ``(dlon, dlat)`` (east/north positive), then converted to an
    array-index offset using the sign of each coordinate axis step (Req 2.4,
    3.1)::

        dcol = dlon * (1 if x_step > 0 else -1)
        drow = dlat * (1 if y_step > 0 else -1)

    A descending y axis (``y_step < 0``, the common netCDF layout where latitude
    decreases as the row index increases) therefore flips "north" to
    ``drow = -1``, so the offset always points to the correct *physical*
    downstream neighbour regardless of axis orientation (Req 3.1, Property 8).

    Only values present in ``fmt.flow`` produce an offset. The pit/sink value,
    the fill/no-data value, and any unrecognised value all return ``None`` -- the
    caller distinguishes these cases.

    :param value: the LDD cell value to map.
    :param fmt: the resolved :class:`FormatSpec` for the river network format.
    :param x_step: signed x-axis step (``x[1] - x[0]``); its sign orients the
        column offset.
    :param y_step: signed y-axis step (``y[1] - y[0]``); its sign orients the
        row offset.
    :return: ``(drow, dcol)`` array-index offset for a flow value, or ``None``
        for a pit, fill, or unrecognised value.
    """
    compass = fmt.flow.get(value)
    if compass is None:
        return None

    dlon, dlat = COMPASS_OFFSETS[compass]
    dcol = dlon * (1 if x_step > 0 else -1)
    drow = dlat * (1 if y_step > 0 else -1)
    return drow, dcol


def _flow_segment(
    row: int,
    col: int,
    drow: int,
    dcol: int,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
) -> Optional[LineString]:
    """Build the flow line segment from a source cell to its downstream cell.

    The segment is a two-vertex :class:`~shapely.geometry.LineString` running
    from the source cell centre ``(x_coords[col], y_coords[row])`` to the
    downstream cell centre ``(x_coords[col + dcol], y_coords[row + drow])``,
    where ``(drow, dcol)`` is the array-index offset produced by
    :func:`_direction_offset` (Req 3.1, 3.2, 3.3). Cell centres are taken
    directly from the 1-D coordinate arrays, so the segment orientation follows
    the physical grid regardless of axis ascending/descending layout.

    When the downstream neighbour falls outside the grid -- the source cell sits
    on the boundary and its flow direction points off the raster -- no segment
    can be drawn and ``None`` is returned so the caller can omit the cell
    (Req 3.7).

    :param row: source cell row index.
    :param col: source cell column index.
    :param drow: row-index offset to the downstream neighbour.
    :param dcol: column-index offset to the downstream neighbour.
    :param x_coords: 1-D array of x/longitude coordinate centres.
    :param y_coords: 1-D array of y/latitude coordinate centres.
    :return: a two-vertex ``LineString`` from the source centre to the
        downstream centre, or ``None`` when the downstream index is out of the
        grid bounds.
    """
    down_row = row + drow
    down_col = col + dcol

    # Omit cells whose downstream neighbour lies off the raster (Req 3.7).
    if not (0 <= down_row < len(y_coords) and 0 <= down_col < len(x_coords)):
        return None

    source = (float(x_coords[col]), float(y_coords[row]))
    downstream = (float(x_coords[down_col]), float(y_coords[down_row]))
    return LineString([source, downstream])


def _pit_cross(
    row: int,
    col: int,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    half_w: float,
    half_h: float,
) -> List[LineString]:
    """Build the two diagonal line segments of a pit cell's X-shaped cross.

    A pit/sink cell has no downstream neighbour, so instead of a flow segment it
    is drawn as an X through its four corners so it stays visible in the output
    (Req 3.4). The cell is centred at ``(cx, cy) = (x_coords[col], y_coords[row])``
    and its four corners sit half a cell width in x and half a cell height in y
    from that centre (Req 3.5). Corner "top"/"bottom" are defined in *spatial*
    terms -- "top" is the larger y value -- so the cross is geometrically correct
    whether the y axis is ascending or descending:

        * top-left     ``(cx - half_w, cy + half_h)``
        * top-right    ``(cx + half_w, cy + half_h)``
        * bottom-left  ``(cx - half_w, cy - half_h)``
        * bottom-right ``(cx + half_w, cy - half_h)``

    Two diagonals are returned (Req 3.4): diagonal A joins top-left to
    bottom-right, diagonal B joins top-right to bottom-left. Both pass through
    the cell centre ``(cx, cy)`` (Req 3.5, 3.6).

    :param row: pit cell row index.
    :param col: pit cell column index.
    :param x_coords: 1-D array of x/longitude coordinate centres.
    :param y_coords: 1-D array of y/latitude coordinate centres.
    :param half_w: half the cell width in x (non-negative).
    :param half_h: half the cell height in y (non-negative).
    :return: a list ``[diagonal_A, diagonal_B]`` of two two-vertex
        ``LineString``s forming the pit cross.
    """
    cx = float(x_coords[col])
    cy = float(y_coords[row])

    top_left = (cx - half_w, cy + half_h)
    top_right = (cx + half_w, cy + half_h)
    bottom_left = (cx - half_w, cy - half_h)
    bottom_right = (cx + half_w, cy - half_h)

    diagonal_a = LineString([top_left, bottom_right])
    diagonal_b = LineString([top_right, bottom_left])
    return [diagonal_a, diagonal_b]


def _build_features(
    ldd: np.ndarray,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    river_network_format: str,
    x_step: float,
    y_step: float,
    half_w: float,
    half_h: float,
) -> Tuple[List[LineString], List[int]]:
    """Classify every LDD cell and assemble parallel geometry/attribute lists.

    Iterates every ``(row, col)`` of the 2-D LDD array (indexed ``[row, col]``
    on the ``(y, x)`` axes) and classifies each cell's value ``v`` against the
    resolved :class:`FormatSpec` for ``river_network_format`` (Req 3.1, 3.2):

    * ``v == fmt.fill`` -- a no-data/fill cell -- is omitted entirely (Req 3.8).
    * ``v == fmt.pit`` -- a pit/sink cell -- contributes the two diagonal
      ``LineString``s of its :func:`_pit_cross`, each carrying attribute ``v``
      (Req 3.4, 3.6, 4.1).
    * ``v`` in ``fmt.flow`` -- a flow cell -- yields an array-index offset via
      :func:`_direction_offset` and a segment via :func:`_flow_segment`. When the
      segment is in-grid it is appended with attribute ``v``; when the downstream
      neighbour falls off the raster the segment is ``None`` and the cell is
      omitted (Req 3.2, 3.7).
    * any other value is unrecognised: the cell is omitted, a warning naming the
      value and its ``(row, col)`` cell is collected, and iteration continues so
      the remaining cells are still processed (Req 3.9).

    After the loop, any collected unrecognised-value warnings are written to
    ``sys.stderr`` (non-fatal -- they do not stop processing or change the exit
    status, Req 3.9).

    The returned lists are parallel: ``attributes[i]`` is the source LDD value of
    ``geometries[i]`` (Req 4.1, 4.2). Each attribute is the integer LDD value of
    the source cell. Returning them in a single deterministic iteration order
    keeps feature count and attribute order reproducible.

    :param ldd: 2-D LDD raster indexed ``[row, col]`` on the ``(y, x)`` axes.
    :param x_coords: 1-D array of x/longitude coordinate centres.
    :param y_coords: 1-D array of y/latitude coordinate centres.
    :param river_network_format: river network format name; a key of
        ``FORMAT_SPECS``.
    :param x_step: signed x-axis step; its sign orients flow column offsets.
    :param y_step: signed y-axis step; its sign orients flow row offsets.
    :param half_w: half the cell width in x, for pit-cross corners.
    :param half_h: half the cell height in y, for pit-cross corners.
    :return: ``(geometries, attributes)`` -- parallel lists of ``LineString``
        geometries and their integer source LDD values.
    """
    fmt = FORMAT_SPECS[river_network_format]

    geometries: List[LineString] = []
    attributes: List[int] = []
    warnings: List[str] = []

    n_rows = ldd.shape[0]
    n_cols = ldd.shape[1]

    for row in range(n_rows):
        for col in range(n_cols):
            v = int(ldd[row, col])

            # No-data / fill cell: omit entirely (Req 3.8).
            if v == fmt.fill:
                continue

            # Pit/sink cell: draw an X cross, each diagonal carrying v (Req 3.4,
            # 3.6, 4.1).
            if v == fmt.pit:
                for diagonal in _pit_cross(
                    row, col, x_coords, y_coords, half_w, half_h
                ):
                    geometries.append(diagonal)
                    attributes.append(v)
                continue

            # Flow cell: map to a downstream offset and build the segment
            # (Req 3.1, 3.2). An out-of-grid downstream neighbour yields a None
            # segment, so the cell is omitted (Req 3.7).
            if v in fmt.flow:
                offset = _direction_offset(v, fmt, x_step, y_step)
                drow, dcol = offset
                segment = _flow_segment(
                    row, col, drow, dcol, x_coords, y_coords
                )
                if segment is not None:
                    geometries.append(segment)
                    attributes.append(v)
                continue

            # Unrecognised value: omit the cell, collect a warning naming the
            # value and cell, and keep processing the rest (Req 3.9).
            warnings.append(
                f"Unrecognised LDD value {v} at cell (row={row}, col={col}) "
                f"for river network format '{river_network_format}'; "
                f"cell omitted."
            )

    # Emit any collected unrecognised-value warnings (non-fatal, Req 3.9).
    for message in warnings:
        print(f"ldd2vec: warning: {message}", file=sys.stderr)

    return geometries, attributes


# Shapefile component extensions written by an "ESRI Shapefile" driver write.
# ``.shp``/``.shx``/``.dbf`` are always produced; ``.prj`` is produced only when
# the GeoDataFrame has a CRS (Req 5.2, 5.3); ``.cpg`` (character encoding) may
# also be produced. Every component of a given base name is moved/cleaned so no
# stale component of a prior set survives an overwrite (Req 7.6).
_SHAPEFILE_EXTENSIONS: Tuple[str, ...] = (".shp", ".shx", ".dbf", ".prj", ".cpg")


def _write_shapefile(
    gdf: "geopandas.GeoDataFrame", output_file: str, crs: Optional[pyproj.CRS]
) -> str:
    """Write the vector network to an ESRI Shapefile set atomically.

    The common base name shared by every output component is derived from
    ``output_file`` by stripping a single trailing ``.shp`` extension, if
    present; all components (``.shp``/``.shx``/``.dbf``/``.prj``/``.cpg``) share
    that base name (Req 7.4). The final ``.shp`` path is ``{base}.shp``.

    The write is performed atomically so a failure never leaves a partial set at
    the target (Req 7.7). The GeoDataFrame is first written into a private
    temporary directory created as a sibling of the target (so it shares the
    target's filesystem and ``os.replace`` moves stay atomic). Only once that
    write succeeds are the produced components moved onto the target base name.
    Before each move the matching pre-existing target component is removed, and
    every target component extension is cleaned first, so no component of a prior
    set with the same base name survives the overwrite (Req 7.6).

    When ``crs`` is ``None`` the GeoDataFrame has an unset CRS and ``to_file``
    naturally omits the ``.prj`` projection definition, so the set is written
    without one (Req 5.2, 5.3). When ``crs`` is supplied it is assigned to the
    GeoDataFrame (overriding any existing value) before writing so the ``.prj``
    records it (Req 5.1).

    On any failure during the write or move, every component already placed at
    the target base name and every temporary component is removed before the
    original exception is re-raised, so the target location is left with no
    partial component (Req 7.7).

    :param gdf: the vector network to write; one ``ldd`` attribute column plus a
        geometry column.
    :param output_file: the requested output shapefile set name; a trailing
        ``.shp`` is stripped to form the shared base name.
    :param crs: the CRS to assign to the output, or ``None`` for an unset CRS
        (no ``.prj`` projection definition).
    :return: the written ``{base}.shp`` path.
    :raises Exception: re-raises any error from deriving paths, writing, or
        moving components, after removing every partial temporary and target
        component so no partial set is left behind (Req 7.7).
    """
    # Derive the shared base name by stripping a single trailing .shp (Req 7.4).
    if output_file.endswith(".shp"):
        base = output_file[: -len(".shp")]
    else:
        base = output_file

    target_dir = os.path.dirname(base) or os.curdir
    target_shp = f"{base}.shp"

    # Assign the CRS so .prj records it; None leaves the CRS unset so to_file
    # omits .prj (Req 5.1, 5.2, 5.3).
    if crs is not None:
        gdf = gdf.set_crs(crs, allow_override=True)

    tmp_dir: Optional[str] = None
    target_touched = False
    try:
        # Write into a private temp directory that is a sibling of the target so
        # it shares the target's filesystem and os.replace stays atomic. The
        # whole write happens here, before any target component is touched, so a
        # failure at this stage cannot disturb a pre-existing set at the target.
        tmp_dir = tempfile.mkdtemp(dir=target_dir)
        tmp_base = os.path.join(tmp_dir, os.path.basename(base))
        tmp_shp = f"{tmp_base}.shp"

        gdf.to_file(tmp_shp, driver="ESRI Shapefile")

        # The temp write succeeded; from here on the target base name is being
        # modified. A failure during the removal/move phase below could leave a
        # partial set, so these steps mark the target as touched and the except
        # handler cleans every target component (Req 7.7).
        target_touched = True

        # Remove every pre-existing target component first so no stale component
        # of a prior set with this base name survives (Req 7.6).
        for ext in _SHAPEFILE_EXTENSIONS:
            existing = f"{base}{ext}"
            if os.path.exists(existing):
                os.remove(existing)

        # Move each component the temp write actually produced onto the target
        # base name (Req 7.5, 7.6).
        for ext in _SHAPEFILE_EXTENSIONS:
            produced = f"{tmp_base}{ext}"
            if os.path.exists(produced):
                os.replace(produced, f"{base}{ext}")
    except Exception:
        # Leave no partial set behind. If the failure happened before the target
        # was touched (e.g. the temp write itself failed) the pre-existing set is
        # left intact; only the temporary components are cleaned (in finally). If
        # the target was being modified, every target component is removed so no
        # partial set survives (Req 7.7).
        if target_touched:
            for ext in _SHAPEFILE_EXTENSIONS:
                target_component = f"{base}{ext}"
                try:
                    if os.path.exists(target_component):
                        os.remove(target_component)
                except OSError:
                    pass
        raise
    finally:
        # Clean up the temporary directory and any leftover temp components.
        if tmp_dir is not None and os.path.isdir(tmp_dir):
            for name in os.listdir(tmp_dir):
                try:
                    os.remove(os.path.join(tmp_dir, name))
                except OSError:
                    pass
            try:
                os.rmdir(tmp_dir)
            except OSError:
                pass

    return target_shp


def ldd_to_vector(
    input_file: str,
    output_file: str = "ldd.shp",
    river_network_format: str = "pcr_d8",
    coords_names: List[str] = [],
) -> str:
    """Convert an LDD raster netCDF into a vector line-network shapefile set.

    Orchestrates the full read -> resolve -> build -> write pipeline (design
    Architecture, high-level flow):

    1. Validate ``river_network_format`` against ``FORMAT_SPECS`` as defence in
       depth behind the argparse ``choices`` restriction (Req 2.3).
    2. Check the input file exists, raising ``ValueError`` naming the path when
       it does not (Req 1.3).
    3. Open the dataset with ``xarray``; any failure to open propagates to the
       caller (``main``) which reports an unreadable-file error (Req 1.4).
    4. Resolve the ``(x_proj, y_proj)`` coordinate names via
       :func:`extract_coords` (Req 6).
    5. Resolve the LDD data variable via :func:`_resolve_ldd_var_name` (Req 1.5).
    6. Resolve the input CRS via :func:`_resolve_crs`; when the warn flag is set
       (a CRS attribute was present but uninterpretable) a non-fatal warning is
       emitted to stderr and processing continues with an unset CRS (Req 5.3).
    7. Read the 2-D LDD array (transposed to ``(y, x)`` so it is indexed
       ``[row, col]``) and the 1-D x/y coordinate arrays, and compute the signed
       axis steps and half-cell sizes via :func:`_axis_steps` (Req 3.3, 3.5).
    8. Build the parallel geometry/attribute lists via :func:`_build_features`
       (Req 3, 4).
    9. Assemble a :class:`geopandas.GeoDataFrame` carrying the geometries, a
       single ``ldd`` integer attribute column, and the resolved CRS (Req 4, 5).
    10. Write the output shapefile set via :func:`_write_shapefile` (Req 7).

    :param input_file: path to the input LDD netCDF file.
    :param output_file: requested output shapefile set name; a trailing ``.shp``
        is stripped to form the shared base name. Defaults to ``ldd.shp``.
    :param river_network_format: river network format name; a key of
        ``FORMAT_SPECS``. Defaults to ``pcr_d8``.
    :param coords_names: optional explicit ``[y, x]`` coordinate variable names;
        empty to auto-detect.
    :return: the written ``{base}.shp`` path.
    :raises ValueError: when ``river_network_format`` is unsupported (Req 2.3),
        the input path does not exist (Req 1.3), the coordinate names cannot be
        resolved (Req 6), or no usable LDD raster is found (Req 1.5).
    """
    # Step 1: validate the format (defence in depth behind argparse, Req 2.3).
    _validate_river_network_format(river_network_format)

    # Step 2: the input file must exist (Req 1.3).
    if not os.path.isfile(input_file):
        raise ValueError(f"Input file not found: '{input_file}'.")

    # Step 3: open the dataset; a failure to open propagates (Req 1.4).
    with xr.open_dataset(input_file) as ds:
        # Step 4: resolve coordinate names (Req 6).
        x_proj, y_proj = extract_coords(ds, coords_names)

        # Step 5: resolve the LDD data variable (Req 1.5).
        ldd_var_name = _resolve_ldd_var_name(ds, x_proj, y_proj)

        # Step 6: resolve the CRS; a present-but-uninterpretable CRS warns but
        # does not abort (Req 5.1, 5.2, 5.3).
        crs, crs_warning = _resolve_crs(ds, ldd_var_name)
        if crs_warning:
            print(
                "ldd2vec: warning: the input CRS could not be interpreted; "
                "the output shapefile set is written with an unset CRS.",
                file=sys.stderr,
            )

        # Step 7: read the LDD array and coordinate arrays, compute axis steps.
        # The LDD variable's dims may be ordered (y, x); transpose so the array
        # is indexed [row=y, col=x] (Req 3.3).
        ldd_da = ds[ldd_var_name].transpose(y_proj, x_proj)
        ldd = np.asarray(ldd_da.values)
        x_coords = np.asarray(ds[x_proj].values)
        y_coords = np.asarray(ds[y_proj].values)
        x_step, y_step, half_w, half_h = _axis_steps(ds, x_proj, y_proj)

        # Step 8: classify every cell and build the geometry/attribute lists
        # (Req 3, 4).
        geometries, attributes = _build_features(
            ldd,
            x_coords,
            y_coords,
            river_network_format,
            x_step,
            y_step,
            half_w,
            half_h,
        )

    # Step 9: assemble the output GeoDataFrame carrying the single ldd attribute
    # column and the resolved CRS (Req 4, 5).
    gdf = geopandas.GeoDataFrame(
        {"ldd": attributes}, geometry=geometries, crs=crs
    )

    # Step 10: write the output shapefile set and return its base .shp path
    # (Req 7).
    return _write_shapefile(gdf, output_file, crs)


def getarg():
    """ Get program arguments.

    Builds the ``ldd2vec`` argparse CLI, mirroring ``lddfix.getarg()``, with
    both a short and a long flag for every argument (Req 8.2). The returned
    namespace carries exactly one attribute per argument, named after the long
    flag's ``dest``: ``input``, ``output``, ``format``, ``coordsnames``.

    Argument-level validation happens here so an invalid invocation is rejected
    before any processing, with a non-zero exit and no output written (Req 8.6):

    - ``-i/--input`` is ``required=True`` so a missing input exits non-zero
      naming the argument before any input is read (Req 1.2, 8.6).
    - ``-f/--format`` is restricted to the keys of ``FORMAT_SPECS`` via argparse
      ``choices`` so an unsupported format exits non-zero before reading the
      input (Req 2.1, 2.2, 2.3).
    - ``-o/--output`` length must be 1-255 characters; otherwise
      ``parser.error(...)`` exits non-zero before any output is written
      (Req 7.1, 7.2, 7.3).
    - ``-E/--coordsnames`` takes one or more coordinate names (``y x`` order);
      the element count is validated downstream in :func:`extract_coords`
      (Req 6.1, 6.4).

    :return: args:  namespace of program arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input LDD netCDF file (ldd.nc)')
    parser.add_argument('-o', '--output', type=str, required=False,
                        default='ldd.shp',
                        help='Output shapefile set name (default: ldd.shp)')
    parser.add_argument('-f', '--format', type=str, required=False,
                        default='pcr_d8', choices=sorted(FORMAT_SPECS.keys()),
                        help='River network format (default: pcr_d8)')
    parser.add_argument('-E', '--coordsnames', type=str, nargs='+', required=False,
                        default=[],
                        help='Coords names for lat, lon (in this order with '
                             'space!) from the netcdf file used')
    args = parser.parse_args()  # assign namespace to args

    # Validate the output file name length is 1-255 chars (Req 7.1, 7.2).
    # parser.error() prints the message to stderr and exits with a non-zero
    # status, so an invalid name is rejected before any output is written.
    if not (1 <= len(args.output) <= 255):
        parser.error(
            f"output file name must be between 1 and 255 characters, got "
            f"{len(args.output)} characters."
        )

    return args

def main():
    """Run ldd2vec from the command line.

    Parses and validates arguments via :func:`getarg` (argparse exits non-zero
    on a missing/invalid argument before any processing runs, so argument
    failures never reach :func:`ldd_to_vector` and no output is written,
    Req 8.6), maps the parsed namespace onto an :func:`ldd_to_vector` call, and
    translates the outcome into a process exit status: ``0`` on success,
    non-zero on any failure (missing/unreadable input, unresolved coordinates,
    no usable LDD raster, write failure, ...). Any such failure is caught here,
    reported to stderr as ``ldd2vec: error: {exc}``, and converted into a
    non-zero return code rather than propagating as an uncaught traceback.
    Returning an int keeps ``main`` testable and lets the ``bin/ldd2vec``
    launcher terminate with ``sys.exit(main())`` (Req 8.6).

    :return: ``0`` on success, ``1`` on any failure.
    """
    args = getarg()
    try:
        ldd_to_vector(
            input_file=args.input,
            output_file=args.output,
            river_network_format=args.format,
            coords_names=args.coordsnames,
        )
    except Exception as exc:
        print(f"ldd2vec: error: {exc}", file=sys.stderr)
        return 1
    return 0
