#   lddfix repairs a local drain direction (LDD) netCDF so the result is
#   hydrologically sound. Every run ends with earthkit.hydro's river network
#   repair as a mandatory final step. Before that final repair the user may
#   apply manual LDD value substitutions, either a single cell coordinate plus a
#   value on the command line, a tab-separated (TSV) file listing many entries,
#   or both.
#
#   Usage:
#   lddfix -i ldd.nc -o ldd_fixed.nc -f pcr_d8 -c <lat> <lon> <value> -t subs.tsv -E y x

import sys
import tempfile
import os
from typing import List, Optional, Tuple
import dataclasses
import xarray as xr
import numpy as np
from earthkit.hydro.river_network import repair as river_network_repair
from lisfloodutilities.mctrivers.helpers import (
    X_COORD_NAMES,
    Y_COORD_NAMES,
    copy_coordinates_and_attributes,
)


# D8 raster formats from earthkit.hydro that encode one direction value per
# cell and therefore support lddfix's per-cell value substitution. The table
# is validated against the --format argument keys (Req 2.1, 2.3), and each
# substitution value is validated against the set for the configured format
# (Req 2.5, 3.5, 4.6).
_POWERS_OF_TWO_D8 = {0, 1, 2, 4, 8, 16, 32, 64, 128}
VALID_LDD_VALUES = {
    "pcr_d8": set(range(1, 10)),    # 1..9 (5 = pit/sink); fill 255
    "esri_d8": _POWERS_OF_TWO_D8,   # 0 + 2**0..2**7 (0 = sink); fill 255
    "merit_d8": _POWERS_OF_TWO_D8,  # same encoding as esri_d8; fill 247
}

# Per-format missing/fill values, relevant when writing output and handling
# no-data cells (Req 7.4).
FILL_VALUES = {
    "pcr_d8": 255,
    "esri_d8": 255,
    "merit_d8": 247,
}


def _validate_river_network_format(river_network_format):
    """Validate a candidate river_network_format against VALID_LDD_VALUES.

    Leaves the input unchanged and raises ValueError naming the unsupported
    value when it is not a supported format. (Req 2.1, 2.3)
    """
    if river_network_format not in VALID_LDD_VALUES:
        supported = ", ".join(sorted(VALID_LDD_VALUES.keys()))
        raise ValueError(
            f"Unsupported river_network_format '{river_network_format}'. "
            f"Supported formats: {supported}."
        )


def extract_coords(ds: xr.Dataset, coords_names: List[str] = []) -> Tuple[str, str]:
    """Resolve the x/longitude and y/latitude coordinate variable names.

    Adapted from ``mctrivers.extract_coords``. When ``coords_names`` is empty
    the names are auto-detected by intersecting the dataset coordinates with
    the package naming conventions ``X_COORD_NAMES``/``Y_COORD_NAMES``. When a
    two-element ``[y, x]`` list is supplied those names are used directly,
    after verifying both are present in the file.

    :param ds: input dataset whose coordinate names are resolved.
    :param coords_names: optional explicit ``[y, x]`` coordinate names.
    :return: ``(x_proj, y_proj)`` coordinate variable names.
    :raises ValueError: when auto-detection is ambiguous (zero or multiple
        matches on an axis), when a supplied name is absent from the file, or
        when ``coords_names`` has an unexpected number of elements. The message
        requests explicit ``-E`` names. (Req 6.2, 6.3, 6.4, 6.5)
    """
    x_proj, y_proj = '', ''
    if len(coords_names) == 0:
        # Auto-detect coordinates by matching against known patterns (Req 6.2)
        available_coords = set(ds.coords)
        x_matches = available_coords.intersection(X_COORD_NAMES)
        y_matches = available_coords.intersection(Y_COORD_NAMES)

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
        # Use provided coordinates (y, x order as per documentation) (Req 6.3)
        y_proj, x_proj = coords_names
        # A supplied name that is absent from the file is rejected (Req 6.4)
        missing = [name for name in (y_proj, x_proj) if name not in ds.variables]
        if missing:
            raise ValueError(
                f"Supplied coordinate name(s) not found in the input file: "
                f"{missing}. Available variables: {list(ds.variables)}. "
                f"Please specify coordinates explicitly using the -E argument "
                f"(format: 'lat lon')."
            )
    else:
        raise ValueError(
            f"coords_names must have 0 or 2 elements, got {len(coords_names)}: "
            f"{coords_names}. Please specify coordinates explicitly using the "
            f"-E argument (format: 'lat lon')."
        )

    return x_proj, y_proj


def _axis_half_resolution(coords: np.ndarray) -> float:
    """Half the grid resolution of a coordinate axis.

    The grid resolution is the absolute spacing between adjacent coordinate
    values; half of it is the selection tolerance used to decide whether a
    requested coordinate falls within a cell (Req 3.2). The sign of the spacing
    (ascending vs descending axis) is irrelevant, so the absolute value is used.

    :param coords: 1-D array of coordinate centre values for a single axis.
    :return: half the absolute spacing between adjacent values.
    :raises ValueError: when the axis has fewer than two coordinate values and
        the resolution cannot be determined.
    """
    if coords.size < 2:
        raise ValueError(
            "Cannot determine grid resolution from an axis with fewer than two "
            "coordinate values."
        )
    return abs(float(coords[1]) - float(coords[0])) / 2.0


def select_cell_index(
    ds: xr.Dataset,
    x_proj: str,
    y_proj: str,
    coord_x: float,
    coord_y: float,
) -> Tuple[int, int]:
    """Map a coordinate pair to a single ``(row, col)`` grid cell index.

    The grid resolution per axis is the absolute spacing between adjacent
    coordinate values. Selection uses ``xarray`` ``.sel(method='nearest')`` with
    ``tolerance`` set to half that resolution on each axis: a coordinate within
    half a cell width of a centre resolves to that cell (Req 3.2), and a pair
    lying exactly on a boundary shared by two cells resolves to the cell whose
    centre is nearest, which is ``method='nearest'`` semantics (Req 3.3). When
    the nearest centre is farther than half a cell width the pair lies outside
    the grid envelope; ``.sel`` raises ``KeyError`` which is translated into an
    out-of-grid ``ValueError`` identifying the coordinate (Req 3.6, 4.7, 5.4).

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :param coord_x: x/longitude component of the coordinate pair.
    :param coord_y: y/latitude component of the coordinate pair.
    :return: ``(row, col)`` integer cell index, where ``row`` is the position
        along the y axis and ``col`` the position along the x axis.
    :raises ValueError: when the coordinate pair lies outside the grid envelope
        (farther than half a cell width from the nearest centre on either axis).
    """
    x_coords = np.asarray(ds[x_proj].values)
    y_coords = np.asarray(ds[y_proj].values)

    x_tol = _axis_half_resolution(x_coords)
    y_tol = _axis_half_resolution(y_coords)

    # A coordinate exactly on the boundary between two cells sits exactly
    # half a cell width from the nearest centre. Floating-point rounding can
    # make that computed distance marginally exceed the raw half-resolution,
    # so .sel would wrongly reject an exact-boundary coordinate as out-of-grid.
    # Widen the sel tolerance by a tiny relative epsilon to admit distances of
    # exactly half a cell (inclusive boundary, Req 3.3) while still rejecting
    # coordinates meaningfully farther than half a cell (genuinely out-of-grid).
    _TOL_EPS = 1e-9
    x_sel_tol = x_tol * (1.0 + _TOL_EPS)
    y_sel_tol = y_tol * (1.0 + _TOL_EPS)

    # Resolve each axis independently so an out-of-grid component can be named.
    try:
        x_sel = ds[x_proj].sel({x_proj: coord_x}, method="nearest", tolerance=x_sel_tol)
    except KeyError:
        raise ValueError(
            f"Coordinate x={coord_x} is outside the input raster grid "
            f"(nearest {x_proj} centre is farther than half a cell width)."
        )
    try:
        y_sel = ds[y_proj].sel({y_proj: coord_y}, method="nearest", tolerance=y_sel_tol)
    except KeyError:
        raise ValueError(
            f"Coordinate y={coord_y} is outside the input raster grid "
            f"(nearest {y_proj} centre is farther than half a cell width)."
        )

    selected_x = x_sel[x_proj].values
    selected_y = y_sel[y_proj].values

    # Translate the selected centre values back to positional indices.
    col = int(np.argwhere(x_coords == selected_x).ravel()[0])
    row = int(np.argwhere(y_coords == selected_y).ravel()[0])

    return row, col


@dataclasses.dataclass
class Substitution:
    """One manual LDD value substitution, resolved to a coordinate pair.

    Instances form an ordered, in-memory list that preserves application order:
    TSV rows first in file order, then the single-coordinate substitution last.
    That ordering is what gives the single-coordinate value its override
    precedence over a TSV row targeting the same cell (Req 5.2).

    :ivar coord_y: coordinate pair component on the y/latitude axis.
    :ivar coord_x: coordinate pair component on the x/longitude axis.
    :ivar value: LDD value to assign at the selected cell.
    :ivar source: provenance tag, ``'tsv:row=<n>'`` or ``'cli'``, used in error
        messages.
    """

    coord_y: float   # coordinate pair component in y/lat
    coord_x: float   # coordinate pair component in x/lon
    value: int       # LDD value to assign
    source: str      # 'tsv:row=<n>' or 'cli'


def read_substitution_tsv(tsv_file: str, valid_values) -> List[Substitution]:
    """Parse a TSV substitution file into an ordered list of ``Substitution``.

    Modelled on ``cutmaps.helpers.read_column_file``: the file is read line by
    line and each line is split on tabs. Every non-blank line must yield exactly
    three tab-separated fields in ``(y, x, value)`` order (i.e. lat, lon, ldd
    value), producing a ``Substitution(coord_y, coord_x, value, source='tsv:row=<n>')``
    where ``<n>`` is the 1-based line number of that row in the file. Blank or
    whitespace-only lines are skipped and do not become substitutions (Req 4.10).

    Parsing is **atomic**: the whole file is scanned first and, if any row is
    malformed, a ``ValueError`` naming the offending 1-based row number and the
    reason is raised and **no** ``Substitution`` objects are returned (Req 7.1,
    7.2, 7.3, 7.5). A row is malformed when it does not have exactly three
    fields (Req 4.5), when any of its fields is non-numeric (Req 4.5), or when
    its integer LDD value is not in ``valid_values`` (Req 2.5, 4.6). Because the
    scan completes before any value is accepted, a later bad row prevents the
    earlier good rows from being applied — nothing is partially applied.

    :param tsv_file: path to the tab-separated substitution file.
    :param valid_values: set/collection of acceptable integer LDD values for the
        configured river network format (``VALID_LDD_VALUES[format]``).
    :return: substitutions in file order, one per non-blank row.
    :raises ValueError: when the file does not exist (naming the path, Req 4.4),
        or when any row is malformed (naming the 1-based row and the reason). No
        substitutions are produced on any such failure.
    """
    if not os.path.isfile(tsv_file):
        raise ValueError(
            f"TSV substitution file not found: '{tsv_file}'."
        )

    try:
        with open(tsv_file) as handle:
            raw_lines = handle.readlines()
    except OSError as exc:
        raise ValueError(
            f"Unable to read TSV substitution file '{tsv_file}': {exc}."
        )

    substitutions: List[Substitution] = []
    # Scan every row first; any malformed row aborts before returning so the
    # parse is atomic (no partial list of substitutions is produced).
    for row_number, raw_line in enumerate(raw_lines, start=1):
        # A blank or whitespace-only line is skipped and is not a substitution.
        if raw_line.strip() == "":
            continue

        # Split on tabs; strip the trailing newline and surrounding whitespace
        # on each field so a well-formed row is parsed regardless of a trailing
        # EOL. Exactly three fields (y, x, value) are required.
        fields = [field.strip() for field in raw_line.rstrip("\n").split("\t")]
        if len(fields) != 3:
            raise ValueError(
                f"Malformed TSV substitution row {row_number}: expected 3 "
                f"tab-separated fields (y, x, value), got {len(fields)}."
            )

        y_field, x_field, value_field = fields

        try:
            coord_y = float(y_field)
            coord_x = float(x_field)
        except ValueError:
            raise ValueError(
                f"Malformed TSV substitution row {row_number}: coordinate "
                f"fields must be numeric, got y='{y_field}', x='{x_field}'."
            )

        # The LDD value must be an integer; reject non-numeric or non-integer
        # (e.g. '3.5') values naming the row.
        try:
            value_float = float(value_field)
        except ValueError:
            raise ValueError(
                f"Malformed TSV substitution row {row_number}: value field "
                f"must be a numeric LDD value, got '{value_field}'."
            )
        if not value_float.is_integer():
            raise ValueError(
                f"Malformed TSV substitution row {row_number}: LDD value must "
                f"be an integer, got '{value_field}'."
            )
        value = int(value_float)

        if value not in valid_values:
            raise ValueError(
                f"Malformed TSV substitution row {row_number}: LDD value "
                f"{value} is not a valid value for the configured format "
                f"(valid values: {sorted(valid_values)})."
            )

        substitutions.append(
            Substitution(
                coord_y=coord_y,
                coord_x=coord_x,
                value=value,
                source=f"tsv:row={row_number}",
            )
        )

    return substitutions


def apply_substitutions(
    ds: xr.Dataset,
    ldd_var_name: str,
    substitutions: List[Substitution],
    x_proj: str,
    y_proj: str,
) -> np.ndarray:
    """Apply an ordered list of validated substitutions to an LDD array.

    Copies the LDD variable's values so the input dataset is never mutated, then
    assigns each substitution's value at the ``(row, col)`` cell returned by
    :func:`select_cell_index`, applying the substitutions in list order. Only
    the targeted cells change; every other cell keeps its original value. When
    two substitutions resolve to the same cell the later one wins, which is how
    a single-coordinate substitution (appended last) overrides a TSV row on the
    same cell (Req 5.2). (Req 3.4, 4.3)

    The substitutions are assumed to be already validated (value in the valid
    set, coordinate in-grid); callers perform those checks up front so a failure
    aborts before any repair runs. Passing an empty list returns an unchanged
    copy of the LDD values.

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param ldd_var_name: name of the LDD data variable within ``ds``.
    :param substitutions: ordered substitutions to apply, in application order.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :return: a new array, a copy of the LDD variable's values, with the
        substitution values assigned at their selected cells.
    """
    # Copy so the input dataset's array is never mutated in place.
    ldd_values = np.array(ds[ldd_var_name].values, copy=True)

    # The LDD variable is a 2-D (y, x) raster; map (row, col) onto its axes.
    ldd_dims = ds[ldd_var_name].dims
    y_axis = ldd_dims.index(y_proj)
    x_axis = ldd_dims.index(x_proj)

    for sub in substitutions:
        row, col = select_cell_index(ds, x_proj, y_proj, sub.coord_x, sub.coord_y)
        # Build an index tuple honouring the variable's actual dimension order
        # so (row, col) target the y and x axes regardless of their layout.
        index: List[int] = [0, 0]
        index[y_axis] = row
        index[x_axis] = col
        ldd_values[tuple(index)] = sub.value

    return ldd_values


def _resolve_ldd_var_name(ds: xr.Dataset, x_proj: str, y_proj: str) -> str:
    """Resolve the name of the LDD data variable in the input dataset.

    The input is a 2-D LDD raster carrying a single direction value per cell, so
    the LDD variable is the lone data variable spanning both the resolved x and y
    coordinate axes. Coordinate variables (``x_proj``/``y_proj``) are excluded.

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :return: name of the LDD data variable.
    :raises ValueError: when zero or more than one candidate data variable spans
        both coordinate axes and the LDD variable cannot be uniquely determined.
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
            f"Unable to uniquely identify the LDD data variable in the input "
            f"file. Found {len(candidates)} candidate variable(s) spanning the "
            f"'{y_proj}'/'{x_proj}' axes: {candidates}. "
            f"Available data variables: {list(ds.data_vars)}."
        )
    return candidates[0]


def _collect_substitutions(
    ds: xr.Dataset,
    x_proj: str,
    y_proj: str,
    river_network_format: str,
    coordinate: Optional[List[float]],
    tsv_file: Optional[str],
) -> List[Substitution]:
    """Assemble and validate the ordered list of manual substitutions.

    The ordered list is TSV rows first (in file order) followed by the single
    single-coordinate CLI substitution last, so the CLI value overrides a TSV
    row targeting the same cell (Req 5.2). Every substitution's value is checked
    against ``VALID_LDD_VALUES[river_network_format]`` and every coordinate is
    resolved via :func:`select_cell_index`, which raises an out-of-grid
    ``ValueError`` for a coordinate pair outside the input grid envelope
    (Req 3.6, 4.7, 5.4). Two TSV rows resolving to the same grid cell are a
    conflict and abort the run (Req 4.9); a CLI/TSV collision is resolved by
    last-wins ordering, not treated as a conflict (Req 5.2). Any failure raises
    before any repair/write, leaving the input unchanged and producing no output
    (Req 3.5, 4.6, 4.8, 5.4).

    :param ds: input dataset carrying the LDD raster and its coordinate axes.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :param river_network_format: the configured, already-validated format.
    :param coordinate: optional CLI single-coordinate substitution as
        ``[coord_y, coord_x, value]`` in ``(y, x, value)`` order (Req 3.1).
    :param tsv_file: optional path to the TSV substitution file.
    :return: ordered, validated substitutions ready to apply (TSV first, CLI
        last). An empty list when no substitution was requested.
    :raises ValueError: on an invalid value, an out-of-grid coordinate, or a
        duplicate TSV cell. No substitution is applied on any such failure.
    """
    valid_values = VALID_LDD_VALUES[river_network_format]

    substitutions: List[Substitution] = []

    # TSV rows first, in file order. read_substitution_tsv already validates
    # each row's structure and value against valid_values (atomic parse).
    if tsv_file is not None:
        substitutions.extend(read_substitution_tsv(tsv_file, valid_values))

    # Then the single-coordinate CLI substitution, appended last so it overrides
    # a TSV row on the same cell (Req 5.2). The CLI value has not been validated
    # yet (the TSV path validates in read_substitution_tsv), so validate it here.
    if coordinate is not None:
        if len(coordinate) != 3:
            raise ValueError(
                f"Single-coordinate substitution must be 3 values "
                f"(y, x, value), got {len(coordinate)}: {coordinate}."
            )
        coord_y, coord_x, raw_value = coordinate
        value_float = float(raw_value)
        if not value_float.is_integer():
            raise ValueError(
                f"Single-coordinate substitution LDD value must be an integer, "
                f"got '{raw_value}'."
            )
        value = int(value_float)
        if value not in valid_values:
            raise ValueError(
                f"Single-coordinate substitution LDD value {value} is not a "
                f"valid value for the configured format '{river_network_format}' "
                f"(valid values: {sorted(valid_values)})."
            )
        substitutions.append(
            Substitution(
                coord_y=float(coord_y),
                coord_x=float(coord_x),
                value=value,
                source="cli",
            )
        )

    # Resolve every coordinate to a grid cell so an out-of-grid coordinate (from
    # either source) aborts before any repair (Req 3.6, 4.7, 5.4), and detect
    # two TSV rows resolving to the same cell as a conflict (Req 4.9). The
    # CLI/TSV collision is intentional override, not a conflict (Req 5.2).
    tsv_cells = {}
    for sub in substitutions:
        try:
            cell = select_cell_index(ds, x_proj, y_proj, sub.coord_x, sub.coord_y)
        except ValueError as exc:
            raise ValueError(
                f"Substitution ({sub.source}) rejected: {exc}"
            )
        if sub.source.startswith("tsv"):
            if cell in tsv_cells:
                raise ValueError(
                    f"Duplicate TSV substitution cell: rows '{tsv_cells[cell]}' "
                    f"and '{sub.source}' both resolve to grid cell "
                    f"(row={cell[0]}, col={cell[1]})."
                )
            tsv_cells[cell] = sub.source

    return substitutions


def getarg():
    """ Get program arguments.

    Builds the ``lddfix`` argparse CLI, mirroring the ``mctrivers`` getarg
    style, with both a short and a long flag for every argument (Req 8.2). The
    returned namespace carries one attribute per argument, named after the long
    flag's ``dest``: ``input``, ``output``, ``format``, ``coordinate``, ``tsv``,
    ``coordsnames``.

    Argument-level validation happens here so an invalid invocation is rejected
    before any repair, with a non-zero exit and no output written (Req 8.5):

    - ``-f/--format`` is restricted to the keys of ``VALID_LDD_VALUES`` via
      argparse ``choices`` so an unsupported format exits non-zero (Req 2.1,
      2.3).
    - ``-o/--output`` length must be 1-255 characters; otherwise
      ``parser.error(...)`` exits non-zero (Req 7.1, 7.2).
    - ``-c/--coordinate`` takes exactly three floats ``lat lon value`` in
      ``(y, x, value)`` order (Req 3.1); argparse ``nargs=3`` rejects a wrong
      arity at parse time.

    :return: args:  namespace of program arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input LDD netCDF file (ldd.nc)')
    parser.add_argument('-o', '--output', type=str, required=False,
                        default='ldd_fixed.nc',
                        help='Output file name (default: ldd_fixed.nc)')
    parser.add_argument('-f', '--format', type=str, required=False,
                        default='pcr_d8', choices=sorted(VALID_LDD_VALUES.keys()),
                        help='River network format (default: pcr_d8)')
    parser.add_argument('-c', '--coordinate', type=float, nargs=3, required=False,
                        default=None, metavar=('LAT', 'LON', 'VALUE'),
                        help='Single-coordinate substitution: lat lon value '
                             '(in y x value order)')
    parser.add_argument('-t', '--tsv', type=str, required=False, default=None,
                        help='TSV substitution file with lat lon value rows')
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


def fix_ldd(
    input_file: str,
    output_file: str = "ldd_fixed.nc",
    river_network_format: str = "pcr_d8",
    coordinate: Optional[List[float]] = None,
    tsv_file: Optional[str] = None,
    coords_names: List[str] = [],
) -> str:
    """Repair an LDD netCDF, optionally applying manual substitutions first.

    :return: path of the written output file
    """
    _validate_river_network_format(river_network_format)

    # Open the input dataset and resolve coordinate variable names (Req 6). The
    # substituted dataset (when any substitutions apply) is built inside this
    # context so its data is captured while the input is open, then staged to a
    # temporary netCDF the final repair reads from.
    with xr.open_dataset(input_file) as ds:
        x_proj, y_proj = extract_coords(ds, coords_names)

        # Resolve the LDD data variable (single value-per-cell raster).
        ldd_var_name = _resolve_ldd_var_name(ds, x_proj, y_proj)

        # Assemble the ordered, validated substitution list (TSV first, CLI
        # last). Any invalid value, out-of-grid coordinate, or duplicate TSV
        # cell raises here, before any repair runs, leaving the input unchanged
        # and writing no output (Req 3.5, 3.6, 4.6, 4.7, 4.8, 4.9, 5.2, 5.4).
        substitutions = _collect_substitutions(
            ds,
            x_proj,
            y_proj,
            river_network_format,
            coordinate,
            tsv_file,
        )

        if substitutions:
            # Apply the validated substitutions to an in-memory copy of the LDD
            # array (the input dataset is never mutated), then build a modified
            # dataset preserving dims, coords and attrs with the LDD variable
            # replaced by the substituted values (Req 3.4, 4.3).
            substituted_values = apply_substitutions(
                ds, ldd_var_name, substitutions, x_proj, y_proj
            )
            modified_ds = ds.copy()
            modified_ds[ldd_var_name] = xr.DataArray(
                substituted_values,
                dims=ds[ldd_var_name].dims,
                coords=ds[ldd_var_name].coords,
                attrs=ds[ldd_var_name].attrs,
            )
            modified_ds.attrs = dict(ds.attrs)
            # Materialise so the data survives after the input dataset closes.
            modified_ds = modified_ds.load()
        else:
            modified_ds = None

    try:
        if modified_ds is not None:
            # Stage the substituted raster to a temporary netCDF so the
            # file-to-file river_network_repair can read it. Mirror the
            # mctrivers tempfile staging approach: create the temp file, write
            # and close it fully before the repair reads it, and remove it in a
            # finally block (Req 2.4, 3.4, 4.3).
            temp_fd, temp_path = tempfile.mkstemp(suffix=".nc")
            os.close(temp_fd)
            try:
                modified_ds.to_netcdf(temp_path)
                modified_ds.close()
                _repair_to_output(
                    temp_path, output_file, river_network_format
                )
            finally:
                if os.path.exists(temp_path):
                    os.remove(temp_path)
        else:
            # No substitutions: repair the input file directly (Req 5.3).
            _repair_to_output(input_file, output_file, river_network_format)
    finally:
        if modified_ds is not None:
            modified_ds.close()

    # river_network_repair writes the output netCDF itself, and its handling of
    # coordinate variables and CRS/grid_mapping metadata is not guaranteed: it
    # renames the LDD variable to 'ldd' and drops the input's coordinate
    # attributes (standard_name/units/...) and any CRS attrs. Reassert the
    # input's coordinate variables and CRS/grid_mapping attributes on the output
    # so the written file stays consistent with the input (Req 1.3, 7.4, 7.5).
    _reassert_output_metadata(input_file, output_file, x_proj, y_proj)

    return output_file


def _reassert_output_metadata(
    input_file: str, output_file: str, x_proj: str, y_proj: str
) -> None:
    """Reassert input coordinate variables and CRS attrs on the repaired output.

    ``river_network_repair`` writes the final netCDF itself and does not
    guarantee preservation of the input's coordinate variables or
    CRS/``grid_mapping`` metadata. This step re-opens the input (read-only, never
    mutated) and the output and, following the
    :func:`mctrivers.helpers.copy_coordinates_and_attributes` pattern, restores
    the input's ``x_proj``/``y_proj`` coordinates (values and attributes) and the
    CRS attributes (``esri_pe_string``, ``spatial_ref``, ``crs_wkt``,
    ``grid_mapping``) plus any ``grid_mapping`` coordinate onto the output's LDD
    data variable. The corrected dataset is then written back to ``output_file``,
    overwriting the repair's output (Req 1.3, 7.4, 7.5).

    Guarding Property 14: if the repair altered the grid envelope or reordered
    axes — the output LDD's 2-D shape differs from the input's, or the
    ``x_proj``/``y_proj`` coordinate axis values are not equal — this function
    **fails loudly** by raising ``RuntimeError`` and removing the now
    inconsistent output file rather than silently writing a file whose grid does
    not match the input (Req 7.4).

    :param input_file: path of the (untouched) input LDD netCDF.
    :param output_file: path of the repaired output netCDF to correct in place.
    :param x_proj: name of the x/longitude coordinate variable.
    :param y_proj: name of the y/latitude coordinate variable.
    :raises RuntimeError: when the output grid shape or coordinate axis values do
        not match the input (grid envelope changed or axes reordered).
    """
    # Open the input read-only; it must never be modified by this step.
    with xr.open_dataset(input_file) as input_ds:
        in_x = np.asarray(input_ds[x_proj].values)
        in_y = np.asarray(input_ds[y_proj].values)

        # Load the full output into memory then close the handle before the
        # rewrite so writing back to the same path never hits an open-file
        # error (relevant on Windows, harmless on Linux).
        with xr.open_dataset(output_file) as out_handle:
            output_ds = out_handle.load()

        # Resolve the LDD data variable in the repaired output. river_network_repair
        # renames it to 'ldd'; fall back to the lone data variable otherwise.
        if "ldd" in output_ds.data_vars:
            out_ldd_name = "ldd"
        elif len(output_ds.data_vars) == 1:
            out_ldd_name = str(next(iter(output_ds.data_vars)))
        else:
            raise RuntimeError(
                f"Unable to identify the repaired LDD variable in the output "
                f"'{output_file}'. Data variables: {list(output_ds.data_vars)}."
            )

        out_ldd = output_ds[out_ldd_name]

        # --- Property 14 guard: grid shape ---------------------------------
        # The input LDD is a (y, x) raster; compare the output LDD's 2-D shape
        # against the input's expected (len(y), len(x)). Fail loudly on any
        # mismatch (grid envelope changed or axes reordered) and remove the
        # inconsistent output rather than keep it.
        expected_shape = (in_y.size, in_x.size)
        if tuple(out_ldd.shape) != expected_shape:
            if os.path.exists(output_file):
                os.remove(output_file)
            raise RuntimeError(
                f"Repaired output grid shape {tuple(out_ldd.shape)} does not "
                f"match the input grid shape {expected_shape} "
                f"(y='{y_proj}', x='{x_proj}'); the repair altered the grid "
                f"envelope or reordered axes. Refusing to write an inconsistent "
                f"output file."
            )

        # --- Property 14 guard: coordinate axis values --------------------
        # The output must carry both coordinate axes with values equal to the
        # input's (exact match, or np.allclose for floating-point coords).
        for axis_name, in_vals in ((x_proj, in_x), (y_proj, in_y)):
            if axis_name not in output_ds.coords and axis_name not in output_ds.variables:
                # Missing axis entirely: cannot verify equality -> inconsistent.
                if os.path.exists(output_file):
                    os.remove(output_file)
                raise RuntimeError(
                    f"Repaired output is missing the '{axis_name}' coordinate "
                    f"axis present in the input; refusing to write an "
                    f"inconsistent output file."
                )
            out_vals = np.asarray(output_ds[axis_name].values)
            if out_vals.shape != in_vals.shape or not _coords_equal(out_vals, in_vals):
                if os.path.exists(output_file):
                    os.remove(output_file)
                raise RuntimeError(
                    f"Repaired output '{axis_name}' coordinate values do not "
                    f"match the input (grid envelope changed or axes "
                    f"reordered). Refusing to write an inconsistent output file."
                )

        # --- Reassert coordinates and CRS/grid_mapping metadata ------------
        # Restore the input's coordinate variables (values + attributes) and CRS
        # attrs onto the output LDD variable, following the shared mctrivers
        # helper pattern (esri_pe_string/spatial_ref/crs_wkt/grid_mapping plus a
        # grid_mapping coordinate when present).
        out_ldd = copy_coordinates_and_attributes(
            source=input_ds, target=out_ldd, y_proj=y_proj, x_proj=x_proj
        )

        # copy_coordinates_and_attributes only copies the CRS attrs it finds on
        # the input's first data variable; also restore the coordinate axes'
        # own attributes (standard_name/long_name/units) that the repair drops.
        for axis_name in (x_proj, y_proj):
            out_ldd[axis_name].attrs = dict(input_ds[axis_name].attrs)

        # Reassert the input's LDD variable name. river_network_repair renames
        # the LDD data variable to 'ldd'; the output must preserve the input's
        # own variable name (e.g. 'ldd_repaired') so the written file stays
        # consistent with the input (Req 1.3, 7.4, 7.5). Resolve the input's
        # LDD variable name and, when the repaired output's variable differs,
        # rename it back and restore the input variable's attributes.
        input_ldd_name = _resolve_ldd_var_name(input_ds, x_proj, y_proj)
        output_ds[out_ldd_name] = out_ldd
        if out_ldd_name != input_ldd_name:
            output_ds = output_ds.rename({out_ldd_name: input_ldd_name})
        # Preserve the input LDD variable's attributes where reasonable (the
        # repair drops them); the CRS/grid_mapping attrs reasserted above by
        # copy_coordinates_and_attributes take precedence, so merge rather than
        # overwrite them.
        merged_attrs = dict(input_ds[input_ldd_name].attrs)
        merged_attrs.update(output_ds[input_ldd_name].attrs)
        output_ds[input_ldd_name].attrs = merged_attrs

    # Rewrite the corrected dataset to the same path, overwriting the repair's
    # output. The output dataset was loaded into memory and its source handle
    # closed above, so the same-path rewrite is safe.
    try:
        output_ds.to_netcdf(output_file)
    finally:
        output_ds.close()


def _coords_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """Compare two coordinate axes for equality.

    Floating-point coordinate axes are compared with ``np.allclose`` to tolerate
    representation round-trips through netCDF; non-floating axes (or any array
    that cannot be compared with ``np.allclose``) fall back to exact ``==``.

    :param a: first coordinate axis values.
    :param b: second coordinate axis values.
    :return: ``True`` when the axes are equal within tolerance.
    """
    if np.issubdtype(a.dtype, np.floating) or np.issubdtype(b.dtype, np.floating):
        return bool(np.allclose(a, b, rtol=0.0, atol=1e-9, equal_nan=True))
    return bool(np.array_equal(a, b))


def _repair_to_output(
    source_file: str, output_file: str, river_network_format: str
) -> None:
    """Run river_network_repair file-to-file, cleaning up on failure.

    ``river_network_repair`` reads ``source_file`` and writes ``output_file``.
    On any repair failure, a partially written output file is removed and the
    original error is re-raised so the run aborts without leaving a partial
    Output File behind (Req 1.7, 7.6).

    :param source_file: path of the netCDF the repair reads from (the input
        file for a repair-only run, or the staged temp file when substitutions
        were applied).
    :param output_file: path of the netCDF the repair writes.
    :param river_network_format: the configured, already-validated format.
    """
    try:
        river_network_repair(
            source_file,
            output_file,
            river_network_format=river_network_format,
            input_source="file",
        )
    except Exception:
        # Remove any partial output so a failed repair never leaves a partially
        # written Output File (Req 1.7, 7.6), then re-raise the original error.
        if os.path.exists(output_file):
            os.remove(output_file)
        raise


def main():
    """Run lddfix from the command line.

    Parses and validates arguments via :func:`getarg` (argparse exits non-zero
    on a missing/invalid argument before any repair runs, so argument failures
    never reach :func:`fix_ldd` and no output is written, Req 1.4, 8.5), maps the
    parsed namespace onto a :func:`fix_ldd` call, and translates the outcome into
    a process exit status: ``0`` on success, non-zero on any validation or repair
    failure (missing/invalid input, invalid substitution, out-of-grid coordinate,
    repair failure, ...). Any such failure is caught here, reported to stderr, and
    converted into a non-zero return code rather than propagating as an uncaught
    traceback (Req 1.7). Returning an int keeps ``main`` testable and lets the
    ``bin/lddfix`` launcher terminate with ``sys.exit(main())`` (Req 8.5).

    :return: ``0`` on success, non-zero on any validation or repair failure.
    """
    args = getarg()
    try:
        fix_ldd(
            input_file=args.input,
            output_file=args.output,
            river_network_format=args.format,
            coordinate=args.coordinate,
            tsv_file=args.tsv,
            coords_names=args.coordsnames,
        )
    except Exception as exc:
        print(f"lddfix: error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
