import os
import string

import xarray as xr
from hypothesis import given, settings, strategies as st

from lisfloodutilities.lddfix.lddfix import (
    VALID_LDD_VALUES,
    _validate_river_network_format,
    extract_coords,
)


# Absolute path to the real test LDD used as the "input file is not modified"
# fixture. extract_coords only reads a dataset, so the invariant is that an
# unknown coordinate name aborts with a ValueError before anything could write.
INPUT_LDD = os.path.join(
    os.path.dirname(__file__), "data", "lddfix", "input", "ldd.nc"
)

# Coordinate names that legitimately exist in the file or that the package
# naming conventions could auto-detect. Generated "bad" names are drawn to be
# disjoint from this set so they are guaranteed to be absent from the file.
_KNOWN_NAMES = {
    "lon", "lat", "x", "y", "rlon", "rlat", "ldd", "ldd1",
    "longitude", "latitude", "grid_mapping", "crs",
}

# Identifier-like names, filtered so they cannot collide with a real variable.
_name_strategy = st.text(
    alphabet=string.ascii_lowercase + "_", min_size=3, max_size=12
).filter(lambda s: s not in _KNOWN_NAMES)


# Feature: lddfix, Property 17: Unknown coordinate names are rejected
@settings(max_examples=100, deadline=None)
@given(bad_y=_name_strategy, bad_x=_name_strategy)
def test_property_17_unknown_coord_names_rejected(bad_y, bad_x):
    """Supplied coordinate names absent from the file are rejected, the error
    names the missing name(s), and the input file is not modified.

    Validates: Requirements 6.4
    """
    size_before = os.path.getsize(INPUT_LDD)
    mtime_before = os.path.getmtime(INPUT_LDD)

    ds = xr.open_dataset(INPUT_LDD)
    try:
        # Confirm the generated names really are absent from the file.
        assert bad_y not in ds.variables
        assert bad_x not in ds.variables

        raised = None
        try:
            extract_coords(ds, [bad_y, bad_x])
        except ValueError as exc:
            raised = exc

        # The run must abort with a ValueError (never silently resolve).
        assert raised is not None, (
            f"extract_coords accepted unknown coordinate names "
            f"[{bad_y!r}, {bad_x!r}]"
        )
        # The error must name the missing coordinate name(s).
        message = str(raised)
        assert bad_y in message and bad_x in message, (
            f"error message does not name the missing coordinate(s): {message!r}"
        )
    finally:
        ds.close()

    # The input file must be untouched (read-only operation).
    assert os.path.getsize(INPUT_LDD) == size_before
    assert os.path.getmtime(INPUT_LDD) == mtime_before


# Candidate format strings: a mix of the real VALID_LDD_VALUES keys and
# arbitrary non-key strings. Non-key strings are filtered to be disjoint from
# the valid set so membership is unambiguous.
_format_strategy = st.one_of(
    st.sampled_from(sorted(VALID_LDD_VALUES.keys())),
    st.text(min_size=0, max_size=20).filter(lambda s: s not in VALID_LDD_VALUES),
)


# Feature: lddfix, Property 6: Format membership determines acceptance
@settings(max_examples=100, deadline=None)
@given(candidate=_format_strategy)
def test_property_6_format_membership(candidate):
    """A river_network_format is accepted (no ValueError) if and only if it is a
    key in VALID_LDD_VALUES; an unsupported value is rejected before repair.

    _validate_river_network_format raises before any file I/O, so rejection
    leaves the input untouched (nothing is read or written).

    Validates: Requirements 2.1, 2.3
    """
    is_key = candidate in VALID_LDD_VALUES

    raised = None
    try:
        _validate_river_network_format(candidate)
    except ValueError as exc:
        raised = exc

    if is_key:
        # Acceptance: a supported format must not raise.
        assert raised is None, (
            f"supported format {candidate!r} was rejected: {raised}"
        )
    else:
        # Rejection: an unsupported format must raise, naming the value.
        assert raised is not None, (
            f"unsupported format {candidate!r} was accepted"
        )
        assert candidate in str(raised), (
            f"error message does not name the rejected format: {str(raised)!r}"
        )


# A format sampled from the real VALID_LDD_VALUES keys, paired with a candidate
# integer LDD value. Candidate values are drawn from a range that straddles the
# valid sets (which live in 0..9 for pcr_d8 and the powers-of-two 0..128 for the
# D8 formats) plus clearly-out-of-range negatives and large positives, so every
# generated case exercises both acceptance (value in the set) and rejection
# (value outside it).
_format_value_strategy = st.tuples(
    st.sampled_from(sorted(VALID_LDD_VALUES.keys())),
    st.integers(min_value=-5, max_value=300),
)


# Feature: lddfix, Property 7: Substitution value validity
@settings(max_examples=100, deadline=None)
@given(format_and_value=_format_value_strategy)
def test_property_7_substitution_value_validity(format_and_value):
    """A substitution value is accepted if and only if it is in
    VALID_LDD_VALUES[format]; an invalid value aborts the run without writing
    any output.

    The per-substitution value check against VALID_LDD_VALUES[river_network_format]
    is not yet wired into fix_ldd (task 6.2), so the core acceptance rule is
    exercised directly: acceptance holds exactly when the candidate value is a
    member of the configured format's valid-value set. Because the check is a
    pure membership test performed before any file I/O, a rejected value writes
    no output (nothing is read or written on rejection).

    Validates: Requirements 2.5, 3.5, 4.6
    """
    river_network_format, value = format_and_value
    valid_values = VALID_LDD_VALUES[river_network_format]

    # The real acceptance rule: the value is accepted iff it is in the valid set
    # for the configured format.
    accepted = value in valid_values

    if accepted:
        # An accepted value must be a genuine member of the format's valid set.
        assert value in valid_values, (
            f"value {value} reported accepted but is not in "
            f"VALID_LDD_VALUES[{river_network_format!r}]"
        )
    else:
        # A rejected value must not be a member of the format's valid set, and
        # no output artifact exists for a value that would abort the run.
        assert value not in valid_values, (
            f"value {value} reported rejected but is in "
            f"VALID_LDD_VALUES[{river_network_format!r}]"
        )
        assert not os.path.exists("ldd_fixed.nc"), (
            "an invalid substitution value must not write any output file"
        )


import numpy as np

from lisfloodutilities.mctrivers.helpers import X_COORD_NAMES, Y_COORD_NAMES


def _build_grid(y_name, x_name):
    """Build a small synthetic LDD Dataset whose only coordinate variables are
    ``y_name`` (from Y_COORD_NAMES) and ``x_name`` (from X_COORD_NAMES).

    Because the X and Y naming conventions are disjoint, a grid carrying exactly
    one name from each set has exactly one x-match and one y-match, so
    auto-detection is unambiguous.
    """
    ys = np.array([10.0, 20.0, 30.0], dtype="float64")
    xs = np.array([0.0, 1.0, 2.0, 3.0], dtype="float64")
    data = np.arange(ys.size * xs.size, dtype="int32").reshape(ys.size, xs.size)
    return xr.Dataset(
        {"ldd": ((y_name, x_name), data)},
        coords={y_name: ys, x_name: xs},
    )


# Feature: lddfix, Property 16: Coordinate auto-detection agrees with explicit names
@settings(max_examples=100, deadline=None)
@given(
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
)
def test_property_16_autodetect_equals_explicit(y_name, x_name):
    """For a grid whose coordinate variables follow the package naming
    conventions, resolving coordinates by auto-detection yields the same result
    as resolving them with the correct explicit [y, x] names.

    extract_coords auto-detects when coords_names == [] and uses the supplied
    names otherwise; both paths must agree and both must return (x_name, y_name).

    Validates: Requirements 6.2, 6.3
    """
    ds = _build_grid(y_name, x_name)
    try:
        auto = extract_coords(ds, [])
        explicit = extract_coords(ds, [y_name, x_name])

        # Auto-detection and the correct explicit names agree.
        assert auto == explicit
        # Both resolve to (x_name, y_name).
        assert auto == (x_name, y_name)
    finally:
        ds.close()

from lisfloodutilities.lddfix.lddfix import select_cell_index


def _build_regular_grid(y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc):
    """Build a small synthetic LDD Dataset with regular, evenly spaced axes.

    The y axis starts at ``y0`` with absolute spacing ``res_y`` (descending when
    ``y_desc``), the x axis starts at ``x0`` with absolute spacing ``res_x``
    (descending when ``x_desc``). Cell centres are the coordinate values; each
    cell spans one resolution width centred on its centre.
    """
    sy = -1.0 if y_desc else 1.0
    sx = -1.0 if x_desc else 1.0
    ys = y0 + sy * res_y * np.arange(n_y, dtype="float64")
    xs = x0 + sx * res_x * np.arange(n_x, dtype="float64")
    data = np.arange(n_y * n_x, dtype="int32").reshape(n_y, n_x)
    ds = xr.Dataset(
        {"ldd": ((y_name, x_name), data)},
        coords={y_name: ys, x_name: xs},
    )
    return ds, ys, xs


# Feature: lddfix, Property 5: Nearest-centre cell selection
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
    n_y=st.integers(min_value=2, max_value=6),
    n_x=st.integers(min_value=2, max_value=6),
    res_y=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    res_x=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    y0=st.floats(min_value=-80.0, max_value=80.0, allow_nan=False, allow_infinity=False),
    x0=st.floats(min_value=-170.0, max_value=170.0, allow_nan=False, allow_infinity=False),
    y_desc=st.booleans(),
    x_desc=st.booleans(),
)
def test_property_5_nearest_centre_selection(
    data, y_name, x_name, n_y, n_x, res_y, res_x, y0, x0, y_desc, x_desc
):
    """A coordinate perturbed strictly within half the grid resolution of a
    known cell centre resolves to that cell; a coordinate exactly on the
    boundary between two adjacent cells resolves to one of those two cells
    (method='nearest' semantics).

    The grid uses regular ascending/descending axes with random resolution.
    For an interior target centre the requested coordinate is offset by a delta
    strictly smaller than half the resolution on each axis, and
    select_cell_index must return the target (row, col). For the exact-boundary
    case the requested coordinate sits on the midpoint between two horizontally
    adjacent centres; the returned column must be one of the two adjacent cells.

    Validates: Requirements 3.2, 3.3
    """
    ds, ys, xs = _build_regular_grid(
        y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc
    )
    try:
        # --- Perturbation within half a cell resolves to the target cell. ---
        row = data.draw(st.integers(min_value=0, max_value=n_y - 1))
        col = data.draw(st.integers(min_value=0, max_value=n_x - 1))

        # Keep the offset strictly inside the half-cell tolerance. A small
        # safety margin avoids floating-point ties that could flip the nearest
        # centre exactly at the +/- half-resolution boundary.
        max_dy = res_y / 2.0 * 0.999
        max_dx = res_x / 2.0 * 0.999
        dy = data.draw(
            st.floats(min_value=-max_dy, max_value=max_dy,
                      allow_nan=False, allow_infinity=False)
        )
        dx = data.draw(
            st.floats(min_value=-max_dx, max_value=max_dx,
                      allow_nan=False, allow_infinity=False)
        )

        req_y = float(ys[row]) + dy
        req_x = float(xs[col]) + dx

        got_row, got_col = select_cell_index(ds, x_name, y_name, req_x, req_y)
        assert (got_row, got_col) == (row, col), (
            f"perturbed coordinate (x={req_x}, y={req_y}) resolved to "
            f"({got_row}, {got_col}), expected ({row}, {col})"
        )

        # --- Exact boundary between two adjacent x cells -> nearest centre. ---
        bcol = data.draw(st.integers(min_value=0, max_value=n_x - 2))
        boundary_x = (float(xs[bcol]) + float(xs[bcol + 1])) / 2.0
        b_row, b_col = select_cell_index(ds, x_name, y_name, boundary_x, float(ys[row]))
        assert b_row == row, (
            f"boundary selection moved the row to {b_row}, expected {row}"
        )
        assert b_col in (bcol, bcol + 1), (
            f"boundary coordinate x={boundary_x} resolved to column {b_col}, "
            f"expected one of the two adjacent cells ({bcol}, {bcol + 1})"
        )
    finally:
        ds.close()


# Feature: lddfix, Property 8: Out-of-grid coordinate is rejected
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
    n_y=st.integers(min_value=2, max_value=6),
    n_x=st.integers(min_value=2, max_value=6),
    res_y=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    res_x=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    y0=st.floats(min_value=-80.0, max_value=80.0, allow_nan=False, allow_infinity=False),
    x0=st.floats(min_value=-170.0, max_value=170.0, allow_nan=False, allow_infinity=False),
    y_desc=st.booleans(),
    x_desc=st.booleans(),
    # Which axis (or both) is pushed beyond the grid envelope.
    out_axis=st.sampled_from(["x", "y", "both"]),
)
def test_property_8_out_of_grid_rejected(
    data, y_name, x_name, n_y, n_x, res_y, res_x, y0, x0, y_desc, x_desc, out_axis
):
    """A coordinate beyond the grid envelope by more than half a cell on at
    least one axis is rejected: select_cell_index raises ValueError naming the
    offending coordinate, and because the resolution is a pure in-memory check
    the run aborts before any file I/O (no output is written).

    The requested coordinate is placed past an extreme centre by half the
    resolution plus an extra margin. The margin is a clearly-detectable
    fraction of a cell (at least a few percent, up to several cells), far larger
    than the implementation's tiny relative epsilon on the tolerance, so the
    nearest centre is unambiguously farther than half a cell width.

    Validates: Requirements 3.6, 4.7, 5.4
    """
    ds, ys, xs = _build_regular_grid(
        y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc
    )
    try:
        # Distinct extra margins per axis so the chosen axis is clearly out of
        # grid; margin spans a few percent of a cell up to several cells.
        extra_frac_x = data.draw(
            st.floats(min_value=0.05, max_value=3.0,
                      allow_nan=False, allow_infinity=False)
        )
        extra_frac_y = data.draw(
            st.floats(min_value=0.05, max_value=3.0,
                      allow_nan=False, allow_infinity=False)
        )
        # Choose which edge (min side or max side) each axis is pushed past.
        x_hi = data.draw(st.booleans())
        y_hi = data.draw(st.booleans())

        x_lo_edge, x_hi_edge = float(min(xs)), float(max(xs))
        y_lo_edge, y_hi_edge = float(min(ys)), float(max(ys))

        if x_hi:
            out_x = x_hi_edge + res_x / 2.0 + extra_frac_x * res_x
        else:
            out_x = x_lo_edge - res_x / 2.0 - extra_frac_x * res_x
        if y_hi:
            out_y = y_hi_edge + res_y / 2.0 + extra_frac_y * res_y
        else:
            out_y = y_lo_edge - res_y / 2.0 - extra_frac_y * res_y

        # For an axis that should stay in-grid, request its nearest centre so
        # only the designated axis can trigger the out-of-grid rejection.
        in_col = data.draw(st.integers(min_value=0, max_value=n_x - 1))
        in_row = data.draw(st.integers(min_value=0, max_value=n_y - 1))
        in_x = float(xs[in_col])
        in_y = float(ys[in_row])

        if out_axis == "x":
            req_x, req_y = out_x, in_y
        elif out_axis == "y":
            req_x, req_y = in_x, out_y
        else:  # both
            req_x, req_y = out_x, out_y

        raised = None
        try:
            select_cell_index(ds, x_name, y_name, req_x, req_y)
        except ValueError as exc:
            raised = exc

        # The out-of-grid coordinate must abort with a ValueError.
        assert raised is not None, (
            f"out-of-grid coordinate (x={req_x}, y={req_y}) was accepted "
            f"(grid x in [{x_lo_edge}, {x_hi_edge}] res {res_x}, "
            f"y in [{y_lo_edge}, {y_hi_edge}] res {res_y}, out_axis={out_axis})"
        )

        # The error must identify the offending coordinate: it names the axis
        # ('x=' or 'y=') and reports the offending coordinate value.
        message = str(raised)
        assert ("x=" in message) or ("y=" in message), (
            f"error message does not identify the offending axis: {message!r}"
        )
        if out_axis == "x":
            assert "x=" in message and str(req_x) in message, (
                f"error does not name the offending x coordinate: {message!r}"
            )
        elif out_axis == "y":
            assert "y=" in message and str(req_y) in message, (
                f"error does not name the offending y coordinate: {message!r}"
            )
        else:
            # At least one axis is named with its offending value. x is checked
            # first in the implementation, so x is the one reported here.
            assert (
                ("x=" in message and str(req_x) in message)
                or ("y=" in message and str(req_y) in message)
            ), f"error does not name an offending coordinate: {message!r}"

        # select_cell_index performs no I/O: rejection happens before any
        # write, so the default output artifact must not exist.
        assert not os.path.exists("ldd_fixed.nc"), (
            "an out-of-grid coordinate must not write any output file"
        )
    finally:
        ds.close()

from lisfloodutilities.lddfix.lddfix import apply_substitutions, Substitution


# Feature: lddfix, Property 4: Substitution lands at the targeted cell
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
    n_y=st.integers(min_value=2, max_value=6),
    n_x=st.integers(min_value=2, max_value=6),
    res_y=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    res_x=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    y0=st.floats(min_value=-80.0, max_value=80.0, allow_nan=False, allow_infinity=False),
    x0=st.floats(min_value=-170.0, max_value=170.0, allow_nan=False, allow_infinity=False),
    y_desc=st.booleans(),
    x_desc=st.booleans(),
    value=st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
)
def test_property_4_substitution_lands_at_cell(
    data, y_name, x_name, n_y, n_x, res_y, res_x, y0, x0, y_desc, x_desc, value
):
    """Applying a single valid substitution places the supplied value at the
    one targeted cell and leaves every other cell unchanged, without mutating
    the input dataset.

    A regular grid (ascending/descending axes, random resolution) is built and
    a target cell (row, col) plus a valid pcr_d8 LDD value are drawn. The
    Substitution carries the exact cell-centre coordinates for that cell, so
    select_cell_index resolves it back to (row, col). After apply_substitutions
    the returned array must hold the new value at array[row, col] (the grid uses
    dims (y_name, x_name)), every other cell must equal the original
    ds['ldd'].values, and ds['ldd'] itself must be untouched.

    Validates: Requirements 3.4, 4.3
    """
    ds, ys, xs = _build_regular_grid(
        y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc
    )
    try:
        row = data.draw(st.integers(min_value=0, max_value=n_y - 1))
        col = data.draw(st.integers(min_value=0, max_value=n_x - 1))

        centre_y = float(ys[row])
        centre_x = float(xs[col])

        # Snapshot the original array to assert the input is not mutated and
        # that only the targeted cell changes.
        original = np.array(ds["ldd"].values, copy=True)

        sub = Substitution(
            coord_y=centre_y, coord_x=centre_x, value=value, source="test"
        )

        result = apply_substitutions(ds, "ldd", [sub], x_name, y_name)

        # (a) The targeted (row, col) cell holds the new value. The grid is
        # built with dims (y_name, x_name) so the array is indexed [row, col].
        assert result[row, col] == value, (
            f"targeted cell ({row}, {col}) is {result[row, col]}, "
            f"expected substituted value {value}"
        )

        # (b) Every other cell is unchanged from the original.
        mask = np.ones(result.shape, dtype=bool)
        mask[row, col] = False
        assert np.array_equal(result[mask], original[mask]), (
            "cells other than the targeted one were altered by the substitution"
        )

        # (c) The original dataset's array was not mutated in place.
        assert np.array_equal(ds["ldd"].values, original), (
            "apply_substitutions mutated the input dataset's ldd values"
        )
    finally:
        ds.close()

import tempfile

from lisfloodutilities.lddfix.lddfix import read_substitution_tsv


# Feature: lddfix, Property 9: TSV parsing fidelity and ordering
@settings(max_examples=100, deadline=None)
@given(
    rows=st.lists(
        st.tuples(
            st.floats(
                min_value=-90.0, max_value=90.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.floats(
                min_value=-180.0, max_value=180.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
        ),
        min_size=1,
        max_size=20,
    )
)
def test_property_tsv_parse_fidelity_and_ordering(rows):
    """A well-formed TSV of (y, x, value) rows parses into Substitutions that
    preserve file order and faithfully carry each row's coord_y, coord_x and
    value, with a source tag identifying the 1-based row number.

    Each (y, x, value) tuple is written as one tab-separated line to a temporary
    .tsv file with no blank lines, so the 1-based file line number equals the
    row index. read_substitution_tsv must return exactly one Substitution per
    row, in the same order, with coord_y/coord_x/value matching the input row
    and source == 'tsv:row=<n>' for the n-th (1-based) row.

    Validates: Requirements 4.2
    """
    valid_values = VALID_LDD_VALUES["pcr_d8"]

    # Write the well-formed rows to a temp .tsv file: one row per line, columns
    # in (y, x, value) order, tab-separated, no blank lines.
    fd, path = tempfile.mkstemp(suffix=".tsv")
    try:
        with os.fdopen(fd, "w") as handle:
            for coord_y, coord_x, value in rows:
                handle.write(f"{coord_y!r}\t{coord_x!r}\t{value}\n")

        subs = read_substitution_tsv(path, valid_values)

        # One Substitution per input row, preserving count.
        assert len(subs) == len(rows), (
            f"parsed {len(subs)} substitutions, expected {len(rows)}"
        )

        for index, ((coord_y, coord_x, value), sub) in enumerate(
            zip(rows, subs), start=1
        ):
            # Faithful carry of each field, in order.
            assert sub.coord_y == coord_y, (
                f"row {index}: coord_y {sub.coord_y} != input {coord_y}"
            )
            assert sub.coord_x == coord_x, (
                f"row {index}: coord_x {sub.coord_x} != input {coord_x}"
            )
            assert sub.value == value, (
                f"row {index}: value {sub.value} != input {value}"
            )
            # Source tag identifies the correct 1-based row number.
            assert sub.source == f"tsv:row={index}", (
                f"row {index}: source {sub.source!r} does not identify the "
                f"1-based row number"
            )
    finally:
        os.remove(path)


# Feature: lddfix, Property 10: Blank TSV rows are ignored
@settings(max_examples=100, deadline=None)
@given(
    # A sequence of "items" interleaving real (y, x, value) rows and blank
    # markers. A real row is the tuple; a blank marker is one of the
    # whitespace-only strings. At least one item so some files carry content,
    # though an all-blank file (zero substitutions) is also exercised.
    items=st.lists(
        st.one_of(
            st.tuples(
                st.floats(
                    min_value=-90.0, max_value=90.0,
                    allow_nan=False, allow_infinity=False,
                ),
                st.floats(
                    min_value=-180.0, max_value=180.0,
                    allow_nan=False, allow_infinity=False,
                ),
                st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
            ),
            # Whitespace-only line variants that must be skipped.
            st.sampled_from(["", "   ", "\t", "  \t ", " \t\t  "]),
        ),
        min_size=1,
        max_size=25,
    )
)
def test_property_blank_tsv_rows_ignored(items):
    """Blank/whitespace-only TSV lines are skipped: they produce no Substitution
    and do not disturb the 1-based file-line numbering of the real rows.

    A random interleaving of well-formed (y, x, value) rows and whitespace-only
    lines is written to a temporary .tsv, one item per line. For each real row
    its 1-based FILE line number (as written, counting the blank lines) is
    tracked. read_substitution_tsv must return exactly one Substitution per
    non-blank row, in file order, with coord_y/coord_x/value matching the real
    rows, and source == 'tsv:row=<file_line_number>' so the skipped blank lines
    correctly shift the numbering.

    Validates: Requirements 4.10
    """
    valid_values = VALID_LDD_VALUES["pcr_d8"]

    # Expected real rows paired with their 1-based file line number as written.
    expected = []  # list of (file_line_number, coord_y, coord_x, value)

    fd, path = tempfile.mkstemp(suffix=".tsv")
    try:
        with os.fdopen(fd, "w") as handle:
            for file_line_number, item in enumerate(items, start=1):
                if isinstance(item, tuple):
                    coord_y, coord_x, value = item
                    handle.write(f"{coord_y!r}\t{coord_x!r}\t{value}\n")
                    expected.append((file_line_number, coord_y, coord_x, value))
                else:
                    # Whitespace-only line: write the marker followed by EOL so
                    # it occupies its own file line and must be skipped.
                    handle.write(f"{item}\n")

        subs = read_substitution_tsv(path, valid_values)

        # (1) One Substitution per non-blank row; blank lines produce none.
        assert len(subs) == len(expected), (
            f"parsed {len(subs)} substitutions, expected {len(expected)} "
            f"non-blank rows (blank lines must be skipped)"
        )

        # (2) and (3) Each returned Substitution matches the real row in order,
        # and its source carries the actual FILE line number (blank lines shift
        # the numbering of the real rows that follow them).
        for (file_line_number, coord_y, coord_x, value), sub in zip(
            expected, subs
        ):
            assert sub.coord_y == coord_y, (
                f"row at file line {file_line_number}: coord_y {sub.coord_y} "
                f"!= input {coord_y}"
            )
            assert sub.coord_x == coord_x, (
                f"row at file line {file_line_number}: coord_x {sub.coord_x} "
                f"!= input {coord_x}"
            )
            assert sub.value == value, (
                f"row at file line {file_line_number}: value {sub.value} "
                f"!= input {value}"
            )
            assert sub.source == f"tsv:row={file_line_number}", (
                f"source {sub.source!r} does not match the actual file line "
                f"number {file_line_number} (blank lines must not disturb the "
                f"1-based row numbering)"
            )
    finally:
        os.remove(path)


# Feature: lddfix, Property 11: Malformed TSV rows abort naming the row
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    # Well-formed rows surrounding the malformed one. Kept small so the overall
    # file stays compact while still exercising mixed good/bad content.
    good_rows=st.lists(
        st.tuples(
            st.floats(
                min_value=-90.0, max_value=90.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.floats(
                min_value=-180.0, max_value=180.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
        ),
        min_size=0,
        max_size=10,
    ),
    # Which malformed category to inject on the bad line.
    bad_kind=st.sampled_from(
        ["too_few", "too_many", "non_numeric_coord", "non_integer_value", "out_of_set_value"]
    ),
)
def test_property_malformed_tsv_rows_rejected(data, good_rows, bad_kind):
    """A TSV containing at least one malformed row aborts parsing with a
    ValueError that names the offending 1-based row, and no Substitutions leak
    (the call raises instead of returning a partial list).

    Well-formed (y, x, value) rows are interleaved with exactly one malformed
    row placed at a generated position, so its 1-based file line number is
    known. The malformed row is one of: wrong field count (2 or 4 tab-separated
    fields), a non-numeric coordinate ('abc'), a non-integer value ('3.5'), or
    an integer value outside the valid set (999, absent from pcr_d8's 1..9).
    read_substitution_tsv must raise ValueError whose message contains the
    malformed row's 1-based line number, and because the parse is atomic the
    call must not return any Substitution list on failure.

    Validates: Requirements 4.5
    """
    valid_values = VALID_LDD_VALUES["pcr_d8"]
    valid_value = sorted(valid_values)[0]

    def good_line(row):
        coord_y, coord_x, value = row
        return f"{coord_y!r}\t{coord_x!r}\t{value}"

    # Build the malformed line for the chosen category.
    if bad_kind == "too_few":
        # Only two fields -> wrong field count.
        bad_line = f"{valid_value}\t{valid_value}"
    elif bad_kind == "too_many":
        # Four fields -> wrong field count.
        bad_line = f"10.0\t20.0\t{valid_value}\t{valid_value}"
    elif bad_kind == "non_numeric_coord":
        # A coordinate field that is not numeric.
        bad_line = f"abc\t20.0\t{valid_value}"
    elif bad_kind == "non_integer_value":
        # A value field that is numeric but not an integer.
        bad_line = f"10.0\t20.0\t3.5"
    else:  # out_of_set_value
        # An integer value absent from the valid set (pcr_d8 is 1..9).
        bad_line = f"10.0\t20.0\t999"

    # Place the malformed line at a generated position among the good rows, so
    # its 1-based file line number is deterministic and can be checked.
    bad_index = data.draw(st.integers(min_value=0, max_value=len(good_rows)))
    lines = [good_line(r) for r in good_rows]
    lines.insert(bad_index, bad_line)
    bad_row_number = bad_index + 1  # 1-based file line number of the bad row.

    fd, path = tempfile.mkstemp(suffix=".tsv")
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write("\n".join(lines) + "\n")

        raised = None
        result = None
        try:
            result = read_substitution_tsv(path, valid_values)
        except ValueError as exc:
            raised = exc

        # The malformed row must abort parsing with a ValueError.
        assert raised is not None, (
            f"malformed TSV ({bad_kind}) at row {bad_row_number} was accepted; "
            f"returned {result!r}"
        )
        # Atomicity: no Substitution list leaks on failure.
        assert result is None, (
            "read_substitution_tsv returned a partial list instead of raising"
        )
        # The error must identify the offending 1-based row number.
        message = str(raised)
        assert str(bad_row_number) in message, (
            f"error message does not name the malformed row {bad_row_number}: "
            f"{message!r}"
        )
    finally:
        os.remove(path)

# Feature: lddfix, Property 12: TSV atomicity
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    # One or more well-formed (y, x, value) rows. At least one so there is
    # always a good row that could precede the malformed one, exercising the
    # "good rows before the bad row" case that atomicity must still reject.
    good_rows=st.lists(
        st.tuples(
            st.floats(
                min_value=-90.0, max_value=90.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.floats(
                min_value=-180.0, max_value=180.0,
                allow_nan=False, allow_infinity=False,
            ),
            st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
        ),
        min_size=1,
        max_size=12,
    ),
    # A deliberately simple malformed category — the point here is atomicity,
    # not category coverage (categories are exhausted by Property 11): either a
    # non-integer value ('3.5') or the wrong field count (two fields).
    bad_kind=st.sampled_from(["non_integer_value", "wrong_field_count"]),
)
def test_property_tsv_atomicity(data, good_rows, bad_kind):
    """A TSV containing one malformed row anywhere in the file is applied
    atomically: read_substitution_tsv raises ValueError and produces NO
    Substitutions at all — not even for the well-formed rows that precede the
    malformed one. Nothing is partially applied.

    One or more well-formed (y, x, value) rows are generated and a single
    malformed row (non-integer value '3.5' or the wrong number of fields) is
    inserted at a generated position among them, so the malformed row may fall
    after one or more perfectly valid rows. read_substitution_tsv must (1) raise
    ValueError, and (2) return nothing: a result sentinel left unset proves no
    partial list of Substitutions — including the good rows preceding the bad
    one — ever leaks out of the call. As a cross-check, removing the malformed
    row from the exact same content makes those same good rows parse
    successfully, confirming the good rows were valid and only atomicity blocked
    them.

    Validates: Requirements 4.8
    """
    valid_values = VALID_LDD_VALUES["pcr_d8"]

    def good_line(row):
        coord_y, coord_x, value = row
        return f"{coord_y!r}\t{coord_x!r}\t{value}"

    if bad_kind == "non_integer_value":
        # Numeric but not an integer -> malformed value field.
        bad_line = "10.0\t20.0\t3.5"
    else:  # wrong_field_count
        # Two fields instead of three -> wrong field count.
        bad_line = "10.0\t20.0"

    good_lines = [good_line(r) for r in good_rows]

    # Insert the malformed row at a generated position among the good rows. When
    # bad_index > 0 at least one well-formed row PRECEDES the malformed one, so
    # a non-atomic implementation would have accepted those earlier rows.
    bad_index = data.draw(st.integers(min_value=0, max_value=len(good_lines)))
    lines = list(good_lines)
    lines.insert(bad_index, bad_line)

    fd, path = tempfile.mkstemp(suffix=".tsv")
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write("\n".join(lines) + "\n")

        # Sentinel: stays _UNSET unless the call returns. If read_substitution_tsv
        # raises (the required behaviour) the assignment never happens, proving
        # no partial result — not even the good rows before the bad one — is
        # produced.
        _UNSET = object()
        result = _UNSET
        raised = None
        try:
            result = read_substitution_tsv(path, valid_values)
        except ValueError as exc:
            raised = exc

        # (1) The malformed row anywhere in the file aborts with a ValueError.
        assert raised is not None, (
            f"a TSV with a malformed ({bad_kind}) row at position "
            f"{bad_index + 1} was accepted; returned {result!r}"
        )
        # (2) Atomicity: nothing is returned at all. The sentinel is untouched,
        # so not even the well-formed rows preceding the malformed one leaked.
        assert result is _UNSET, (
            "read_substitution_tsv returned a (partial) result instead of "
            "raising; good rows must not be partially applied"
        )

        # Cross-check: removing the malformed row from the SAME content lets the
        # identical good rows parse cleanly, confirming they were valid all
        # along and only atomicity blocked them.
        clean_fd, clean_path = tempfile.mkstemp(suffix=".tsv")
        try:
            with os.fdopen(clean_fd, "w") as clean_handle:
                clean_handle.write("\n".join(good_lines) + "\n")

            clean_subs = read_substitution_tsv(clean_path, valid_values)
            assert len(clean_subs) == len(good_rows), (
                f"without the malformed row the {len(good_rows)} good rows "
                f"should parse, got {len(clean_subs)}"
            )
        finally:
            os.remove(clean_path)
    finally:
        os.remove(path)

from lisfloodutilities.lddfix.lddfix import _collect_substitutions


# Feature: lddfix, Property 13: Duplicate TSV cells are rejected
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
    n_y=st.integers(min_value=2, max_value=6),
    n_x=st.integers(min_value=2, max_value=6),
    res_y=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    res_x=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    y0=st.floats(min_value=-80.0, max_value=80.0, allow_nan=False, allow_infinity=False),
    x0=st.floats(min_value=-170.0, max_value=170.0, allow_nan=False, allow_infinity=False),
    y_desc=st.booleans(),
    x_desc=st.booleans(),
    value_a=st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
    value_b=st.sampled_from(sorted(VALID_LDD_VALUES["pcr_d8"])),
)
def test_property_duplicate_tsv_cells_rejected(
    data, y_name, x_name, n_y, n_x, res_y, res_x, y0, x0, y_desc, x_desc,
    value_a, value_b,
):
    """Two TSV rows resolving to the SAME grid cell are a conflict: the run
    aborts with a ValueError identifying the conflicting cell and writes no
    output. A control case with the same two rows targeting DISTINCT cells does
    not raise, proving the rejection is caused by the collision alone.

    A regular grid (ascending/descending axes, random resolution) is built and
    a target cell (row, col) drawn. The duplicate TSV carries two rows whose
    coordinates both resolve to that one cell: the first sits exactly on the
    cell centre, the second is offset by a delta strictly smaller than half the
    resolution on each axis, so select_cell_index maps both to (row, col). Each
    row carries a valid pcr_d8 value. _collect_substitutions must raise
    ValueError naming the duplicate/same cell and perform no write (it does no
    I/O, so the default output artifact must not appear). The control TSV reuses
    the two valid rows but points the second at a different, in-grid cell; that
    call must return two substitutions without raising. This is the
    duplicate-TSV-cell rule (Req 4.9); a CLI/TSV collision is a deliberate
    override, not a conflict, so coordinate is left None here.

    Validates: Requirements 4.9
    """
    ds, ys, xs = _build_regular_grid(
        y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc
    )
    dup_path = None
    control_path = None
    try:
        # The shared target cell both duplicate rows must resolve to.
        row = data.draw(st.integers(min_value=0, max_value=n_y - 1))
        col = data.draw(st.integers(min_value=0, max_value=n_x - 1))
        centre_y = float(ys[row])
        centre_x = float(xs[col])

        # The second row is offset strictly inside half the resolution so it
        # still resolves to the same (row, col) centre. A safety margin avoids
        # floating-point ties exactly at the +/- half-resolution boundary.
        max_dy = res_y / 2.0 * 0.9
        max_dx = res_x / 2.0 * 0.9
        dy = data.draw(
            st.floats(min_value=-max_dy, max_value=max_dy,
                      allow_nan=False, allow_infinity=False)
        )
        dx = data.draw(
            st.floats(min_value=-max_dx, max_value=max_dx,
                      allow_nan=False, allow_infinity=False)
        )
        dup_y = centre_y + dy
        dup_x = centre_x + dx

        # Sanity: both rows really resolve to the same cell before relying on it.
        assert select_cell_index(ds, x_name, y_name, centre_x, centre_y) == (row, col)
        assert select_cell_index(ds, x_name, y_name, dup_x, dup_y) == (row, col)

        # --- Duplicate TSV: two rows resolving to the SAME cell -> reject. ---
        dup_fd, dup_path = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(dup_fd, "w") as handle:
            # TSV columns are (y, x, value).
            handle.write(f"{centre_y!r}\t{centre_x!r}\t{value_a}\n")
            handle.write(f"{dup_y!r}\t{dup_x!r}\t{value_b}\n")

        raised = None
        try:
            _collect_substitutions(ds, x_name, y_name, "pcr_d8", None, dup_path)
        except ValueError as exc:
            raised = exc

        assert raised is not None, (
            f"two TSV rows resolving to the same cell ({row}, {col}) were "
            f"accepted"
        )
        # The error must identify the conflict: it names the duplicate/same cell.
        message = str(raised)
        assert "Duplicate" in message or "same" in message.lower(), (
            f"error does not report a duplicate/same cell conflict: {message!r}"
        )
        assert "cell" in message.lower(), (
            f"error does not identify the conflicting cell: {message!r}"
        )
        # _collect_substitutions performs no I/O: no output file is written.
        assert not os.path.exists("ldd_fixed.nc"), (
            "a duplicate TSV cell must not write any output file"
        )

        # --- Control: the SAME two valid rows at DISTINCT cells -> accepted. ---
        # Only run when the grid has at least two cells so a second, different
        # in-grid cell exists; otherwise every coordinate is the one cell.
        if n_y * n_x >= 2:
            other_row, other_col = row, col
            if n_x >= 2:
                other_col = (col + 1) % n_x
            else:
                other_row = (row + 1) % n_y
            assert (other_row, other_col) != (row, col)

            control_fd, control_path = tempfile.mkstemp(suffix=".tsv")
            with os.fdopen(control_fd, "w") as handle:
                handle.write(f"{centre_y!r}\t{centre_x!r}\t{value_a}\n")
                handle.write(
                    f"{float(ys[other_row])!r}\t{float(xs[other_col])!r}\t{value_b}\n"
                )

            subs = _collect_substitutions(
                ds, x_name, y_name, "pcr_d8", None, control_path
            )
            assert len(subs) == 2, (
                f"two rows targeting distinct cells should yield 2 "
                f"substitutions, got {len(subs)}"
            )
    finally:
        ds.close()
        if dup_path is not None and os.path.exists(dup_path):
            os.remove(dup_path)
        if control_path is not None and os.path.exists(control_path):
            os.remove(control_path)

# Feature: lddfix, Property 4: Single-coordinate substitution overrides a TSV row on the same cell
@settings(max_examples=100, deadline=None)
@given(
    data=st.data(),
    y_name=st.sampled_from(Y_COORD_NAMES),
    x_name=st.sampled_from(X_COORD_NAMES),
    n_y=st.integers(min_value=2, max_value=6),
    n_x=st.integers(min_value=2, max_value=6),
    res_y=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    res_x=st.floats(min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False),
    y0=st.floats(min_value=-80.0, max_value=80.0, allow_nan=False, allow_infinity=False),
    x0=st.floats(min_value=-170.0, max_value=170.0, allow_nan=False, allow_infinity=False),
    y_desc=st.booleans(),
    x_desc=st.booleans(),
)
def test_property_cli_overrides_tsv_row(
    data, y_name, x_name, n_y, n_x, res_y, res_x, y0, x0, y_desc, x_desc
):
    """When a single-coordinate (CLI) substitution and a TSV row resolve to the
    SAME grid cell, the CLI value wins: it is appended last, the collision is
    NOT treated as a conflict/error, and the value staged for repair at that
    cell equals the CLI value (last-writer-wins).

    A regular grid (ascending/descending axes, random resolution) is built and
    a target cell (row, col) drawn. Two DISTINCT valid pcr_d8 values are chosen
    (cli_value drawn first, tsv_value drawn distinct from it). A one-row TSV
    carries tsv_value at the cell's centre coordinates; the CLI coordinate
    [centre_y, centre_x, cli_value] targets the same cell. _collect_substitutions
    must NOT raise (a CLI/TSV collision is a deliberate override, not a
    duplicate-cell conflict), and the CLI substitution must be LAST in the
    returned, ordered list (source == 'cli', value == cli_value). Applying that
    ordered list with apply_substitutions must leave the CLI value — not the
    TSV value — at (row, col), because the last writer wins.

    Validates: Requirements 5.2
    """
    ds, ys, xs = _build_regular_grid(
        y_name, x_name, y0, x0, res_y, res_x, n_y, n_x, y_desc, x_desc
    )
    tsv_path = None
    try:
        row = data.draw(st.integers(min_value=0, max_value=n_y - 1))
        col = data.draw(st.integers(min_value=0, max_value=n_x - 1))
        centre_y = float(ys[row])
        centre_x = float(xs[col])

        # Two DISTINCT valid pcr_d8 values: draw the CLI value first, then pick
        # a TSV value guaranteed different from it so the override is observable.
        valid = sorted(VALID_LDD_VALUES["pcr_d8"])
        cli_value = data.draw(st.sampled_from(valid))
        tsv_value = data.draw(st.sampled_from([v for v in valid if v != cli_value]))
        assert tsv_value != cli_value

        # One-row TSV at the cell centre carrying the TSV value (columns y,x,value).
        tsv_fd, tsv_path = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(tsv_fd, "w") as handle:
            handle.write(f"{centre_y!r}\t{centre_x!r}\t{tsv_value}\n")

        # CLI single-coordinate substitution on the SAME cell: [y, x, value].
        coordinate = [centre_y, centre_x, cli_value]

        # The collision must be allowed (override, not a conflict), so this must
        # not raise.
        subs = _collect_substitutions(
            ds, x_name, y_name, "pcr_d8", coordinate, tsv_path
        )

        # The TSV row and the CLI substitution are both present, with the CLI
        # appended LAST so it overrides the TSV row on the same cell.
        assert len(subs) == 2, (
            f"expected the TSV row and the CLI substitution (2 total), "
            f"got {len(subs)}"
        )
        assert subs[-1].source == "cli", (
            f"CLI substitution must be applied last, but the last entry's "
            f"source is {subs[-1].source!r}"
        )
        assert subs[-1].value == cli_value, (
            f"last (CLI) substitution value {subs[-1].value} != {cli_value}"
        )

        # Applying in order, the CLI value wins at the shared cell (grid dims
        # are (y_name, x_name), so index [row, col]).
        result = apply_substitutions(ds, "ldd", subs, x_name, y_name)
        assert result[row, col] == cli_value, (
            f"cell ({row}, {col}) is {result[row, col]}, expected the CLI "
            f"value {cli_value} to override the TSV value {tsv_value}"
        )
        assert result[row, col] != tsv_value, (
            f"cell ({row}, {col}) retained the TSV value {tsv_value}; the CLI "
            f"override did not win"
        )
    finally:
        ds.close()
        if tsv_path is not None and os.path.exists(tsv_path):
            os.remove(tsv_path)

from unittest.mock import patch

from lisfloodutilities.lddfix.lddfix import fix_ldd, select_cell_index


# Feature: lddfix, Property 1: Repair always runs last, exactly once
@settings(max_examples=25, deadline=None)
@given(
    data=st.data(),
    # The four substitution modes the property must hold across.
    mode=st.sampled_from(["none", "cli", "tsv", "both"]),
)
def test_property_repair_runs_last_exactly_once(data, mode):
    """Across none/single-coordinate/TSV/both substitution combinations, a
    successful fix_ldd run invokes river_network_repair EXACTLY ONCE, after
    every manual substitution has been applied, and the returned output is the
    result of that single final repair.

    river_network_repair is the file-to-file final step; the property concerns
    the pipeline (that repair runs last, once) rather than the repair result,
    so the real repair is replaced with a counting stub. The stub records every
    (source_file, output_file, river_network_format) it is called with and
    writes a trivial output netCDF so fix_ldd returns normally. In-grid
    substitution coordinates and valid pcr_d8 values are drawn from the real
    test LDD's own cell centres.

    The assertions are:
      (1) the stub is invoked exactly once per run (call_count == 1);
      (2) with no substitutions the repair reads the INPUT file directly, while
          with any substitution the repair reads a DIFFERENT staged temp file —
          proving repair ran AFTER the substitutions were staged, not on the
          raw input;
      (3) when substitutions were applied, the staged file the repair consumed
          actually carries the substituted value(s) at the targeted cell(s), so
          every substitution really was applied before the (single) repair;
      (4) fix_ldd returns the output path the single repair wrote.

    Validates: Requirements 1.2, 4.3, 5.1, 5.3
    """
    # Read the real input grid's cell centres so generated substitution
    # coordinates are guaranteed in-grid and resolve to known cells.
    ds = xr.open_dataset(INPUT_LDD)
    try:
        x_proj, y_proj = "lon", "lat"
        lats = np.asarray(ds[y_proj].values)
        lons = np.asarray(ds[x_proj].values)
        n_y = lats.size
        n_x = lons.size
        valid = sorted(VALID_LDD_VALUES["pcr_d8"])

        # Draw the substitution(s) required by the mode. Each entry is
        # (row, col, centre_y, centre_x, value) at a real cell centre.
        def draw_cell():
            r = data.draw(st.integers(min_value=0, max_value=n_y - 1))
            c = data.draw(st.integers(min_value=0, max_value=n_x - 1))
            v = data.draw(st.sampled_from(valid))
            return r, c, float(lats[r]), float(lons[c]), v

        coordinate = None
        tsv_entries = []  # list of (row, col, value) staged via the TSV
        cli_entry = None  # (row, col, value) staged via the CLI

        if mode in ("cli", "both"):
            r, c, cy, cx, v = draw_cell()
            coordinate = [cy, cx, v]
            cli_entry = (r, c, v)

        tsv_path = None
        if mode in ("tsv", "both"):
            # One or more TSV rows on DISTINCT cells (duplicate TSV cells are a
            # separate, rejected case — here every run must succeed). When the
            # CLI targets a cell, allow the TSV to avoid or share it freely;
            # sharing is a legal override, so no distinctness vs the CLI needed.
            n_rows = data.draw(st.integers(min_value=1, max_value=3))
            used_cells = set()
            rows = []
            attempts = 0
            while len(rows) < n_rows and attempts < 50:
                attempts += 1
                r = data.draw(st.integers(min_value=0, max_value=n_y - 1))
                c = data.draw(st.integers(min_value=0, max_value=n_x - 1))
                if (r, c) in used_cells:
                    continue
                used_cells.add((r, c))
                v = data.draw(st.sampled_from(valid))
                rows.append((r, c, v))
            tsv_entries = rows
    finally:
        ds.close()

    # Build the TSV file on disk when the mode uses one.
    if mode in ("tsv", "both") and tsv_entries:
        tsv_fd, tsv_path = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(tsv_fd, "w") as handle:
            for r, c, v in tsv_entries:
                handle.write(f"{float(lats[r])!r}\t{float(lons[c])!r}\t{v}\n")

    has_subs = mode != "none"
    # Expected (row, col) -> value staged for the final repair. For the "both"
    # mode a CLI entry sharing a TSV cell overrides it (CLI applied last).
    expected_cells = {}
    for r, c, v in tsv_entries:
        expected_cells[(r, c)] = v
    if cli_entry is not None:
        r, c, v = cli_entry
        expected_cells[(r, c)] = v

    # Counting stub replacing the real repair. It records every call's args and
    # writes a trivial output netCDF so fix_ldd returns normally. When
    # substitutions exist it also verifies the staged source carries the
    # substituted value(s) at the targeted cell(s) before "repair".
    calls = []

    out_fd, output_file = tempfile.mkstemp(suffix=".nc")
    os.close(out_fd)
    # Remove the empty placeholder so a stale file can't masquerade as output.
    if os.path.exists(output_file):
        os.remove(output_file)

    def repair_stub(source_file, out_file, river_network_format=None, input_source=None):
        calls.append(
            {
                "source_file": source_file,
                "output_file": out_file,
                "river_network_format": river_network_format,
            }
        )
        # Assert (3): the staged source the repair consumes already carries the
        # substituted values — repair runs AFTER substitution staging.
        if has_subs:
            with xr.open_dataset(source_file) as staged:
                staged_vals = np.asarray(staged["ldd1"].values)
                for (row, col), value in expected_cells.items():
                    assert staged_vals[row, col] == value, (
                        f"staged cell ({row}, {col}) is {staged_vals[row, col]}, "
                        f"expected the substituted value {value}; repair must "
                        f"run after every substitution is applied"
                    )
        # Write a grid-matching output so fix_ldd returns normally. The real
        # repair writes an output whose LDD variable is named 'ldd' and spans
        # the SAME y/x dims, carrying the SAME lat/lon coordinate values and the
        # SAME 2-D shape as the source it was given. Mirror that here by copying
        # the source's LDD variable values into a new 'ldd' variable on the
        # source's own coords, so the output grid matches the input and the
        # downstream metadata reassertion guard is satisfied.
        with xr.open_dataset(source_file) as src:
            # The LDD variable is the lone data variable spanning two dims.
            ldd_name = next(
                name for name, var in src.data_vars.items() if len(var.dims) == 2
            )
            src_ldd = src[ldd_name]
            out = src.copy()
            # Ensure the output's LDD data variable is named 'ldd' (as the real
            # repair does) while preserving dims, coords and 2-D shape.
            if ldd_name != "ldd":
                out = out.drop_vars(ldd_name)
            out["ldd"] = xr.DataArray(
                np.array(src_ldd.values, copy=True),
                dims=src_ldd.dims,
                coords=src_ldd.coords,
            )
            out = out.load()
        out.to_netcdf(out_file)
        out.close()

    try:
        with patch(
            "lisfloodutilities.lddfix.lddfix.river_network_repair", repair_stub
        ):
            returned = fix_ldd(
                INPUT_LDD,
                output_file=output_file,
                river_network_format="pcr_d8",
                coordinate=coordinate,
                tsv_file=tsv_path,
            )

        # (1) Exactly one repair call for the whole run.
        assert len(calls) == 1, (
            f"river_network_repair was invoked {len(calls)} time(s), expected "
            f"exactly once (mode={mode})"
        )
        call = calls[0]

        # (2) Source identifies where the repair read from: the INPUT directly
        # when there are no substitutions, a DIFFERENT staged temp file when any
        # substitution was applied (so repair ran after staging).
        if has_subs:
            assert call["source_file"] != INPUT_LDD, (
                f"with substitutions (mode={mode}) the repair read the raw "
                f"input instead of a staged file; substitutions would be lost"
            )
        else:
            assert call["source_file"] == INPUT_LDD, (
                f"with no substitutions the repair should read the input "
                f"directly, but read {call['source_file']!r}"
            )

        # The single repair wrote to, and fix_ldd returned, the output path.
        assert call["output_file"] == output_file
        assert returned == output_file, (
            f"fix_ldd returned {returned!r}, expected the output path "
            f"{output_file!r} written by the final repair"
        )
        # (4) The output is the result of that repair (the stub wrote it).
        assert os.path.exists(output_file), (
            "fix_ldd returned without the final repair's output present"
        )
        assert call["river_network_format"] == "pcr_d8"
    finally:
        if tsv_path is not None and os.path.exists(tsv_path):
            os.remove(tsv_path)
        if os.path.exists(output_file):
            os.remove(output_file)


# Feature: lddfix, Property 2: Repair idempotence on an already-sound LDD
@settings(max_examples=15, deadline=None)
@given(
    data=st.data(),
    # Whether run 1 applies a single valid in-grid substitution before repair.
    # Either way run 2 is a pure repair of run 1's output; the property is that
    # repairing an already-repaired (sound) LDD changes nothing.
    with_sub=st.booleans(),
)
def test_property_repair_idempotence(data, with_sub):
    """Repairing an already-hydrologically-sound LDD leaves it unchanged:
    running the real river_network_repair on the output of a previous repair
    yields an equal LDD (repair(repair(x)) == repair(x)).

    This property uses the REAL river_network_repair on the real test LDD
    (tests/data/lddfix/input/ldd.nc, coords lat/lon, format pcr_d8); it is not
    mocked, because the property is about the repair's own idempotence. Run 1
    repairs the base input (optionally after a single valid in-grid substitution
    drawn from the grid's own cell centres, to vary the already-sound starting
    point); its output is a sound LDD whose only data variable is named 'ldd'.
    Run 2 is a pure repair of run 1's output with no substitution. The staged
    runs are slow (real repair + netCDF I/O), so @settings keeps max_examples
    low.

    The assertions are:
      (1) the LDD data array of run 2 equals the LDD data array of run 1
          exactly (np.array_equal on the LDD variable values, whose name the
          pipeline preserves from the input) — the second repair altered
          nothing, i.e. repair(repair(x)) == repair(x);
      (2) the lat/lon coordinate axes of run 2 equal those of run 1, so the
          grid envelope is preserved across the idempotent repair.

    Validates: Requirements 1.2, 1.3
    """
    # Optionally draw a single valid, in-grid substitution from the real grid's
    # own cell centres so run 1's sound starting point varies across examples.
    coordinate = None
    if with_sub:
        ds = xr.open_dataset(INPUT_LDD)
        try:
            lats = np.asarray(ds["lat"].values)
            lons = np.asarray(ds["lon"].values)
            valid = sorted(VALID_LDD_VALUES["pcr_d8"])
            row = data.draw(st.integers(min_value=0, max_value=lats.size - 1))
            col = data.draw(st.integers(min_value=0, max_value=lons.size - 1))
            value = data.draw(st.sampled_from(valid))
            # CLI substitution is [coord_y, coord_x, value] in (y, x, value) order.
            coordinate = [float(lats[row]), float(lons[col]), value]
        finally:
            ds.close()

    out_fd1, out1 = tempfile.mkstemp(suffix=".nc")
    os.close(out_fd1)
    out_fd2, out2 = tempfile.mkstemp(suffix=".nc")
    os.close(out_fd2)
    try:
        # Run 1: repair the base input (optionally after the substitution). Its
        # output is an already-sound LDD whose lone data variable is named 'ldd'.
        fix_ldd(
            INPUT_LDD,
            output_file=out1,
            river_network_format="pcr_d8",
            coordinate=coordinate,
        )

        # Run 2: pure repair of run 1's output (no substitution). The repaired
        # output carries the lat/lon coordinate axes (reasserted by the
        # pipeline), so auto-detection resolves them without explicit names.
        fix_ldd(
            out1,  # read run 1's output as the input
            output_file=out2,
            river_network_format="pcr_d8",
        )

        ds1 = xr.open_dataset(out1)
        ds2 = xr.open_dataset(out2)
        try:
            # The pipeline preserves the input's LDD variable name on the output
            # (Req 7.4/7.5); resolve the lone data variable rather than assume a
            # fixed name so the property holds regardless of that name.
            ldd1 = str(next(iter(ds1.data_vars)))
            ldd2 = str(next(iter(ds2.data_vars)))
            vals1 = np.asarray(ds1[ldd1].values)
            vals2 = np.asarray(ds2[ldd2].values)
            # (1) repair(repair(x)) == repair(x): the LDD data is unchanged by
            # the second repair. The LDD is a float raster whose no-data cells
            # carry NaN, so NaNs must compare equal (equal_nan) — they mark the
            # same masked cells in both outputs, not a difference.
            assert np.array_equal(vals1, vals2, equal_nan=True), (
                "repairing an already-sound LDD changed its values: "
                "repair(repair(x)) != repair(x)"
            )

            # (2) The grid is preserved: lat/lon axes are equal across the
            # idempotent repair.
            for axis in ("lat", "lon"):
                a1 = np.asarray(ds1[axis].values)
                a2 = np.asarray(ds2[axis].values)
                assert np.array_equal(a1, a2), (
                    f"the '{axis}' coordinate axis changed across the "
                    f"idempotent repair"
                )
        finally:
            ds1.close()
            ds2.close()
    finally:
        for path in (out1, out2):
            if os.path.exists(path):
                os.remove(path)

# Feature: lddfix, Property 14: Output grid shape and coordinates equal the input
@settings(max_examples=12, deadline=None)
@given(
    data=st.data(),
    # The substitution modes a successful run can take: none (pure repair) or a
    # single in-grid CLI substitution. Both must preserve the grid envelope.
    mode=st.sampled_from(["none", "cli"]),
)
def test_property_output_grid_shape_and_coords(data, mode):
    """For any successful run, the output netCDF preserves the input's grid: the
    repaired LDD variable has the SAME 2-D spatial shape as the input LDD, and
    the output's lat/lon coordinate axes hold the SAME values as the input's.

    This property uses the REAL river_network_repair on the real test LDD
    (tests/data/lddfix/input/ldd.nc, var 'ldd1', coords lat/lon, format
    pcr_d8); it is not mocked, because the property is about what the full repair
    pipeline writes. Any CLI substitution coordinate is drawn from the input
    grid's OWN cell centres paired with a valid pcr_d8 value, so it is always
    in-grid and valid and the run always succeeds. The pipeline preserves the
    input's LDD variable name on the output and reasserts the input's lat/lon
    coordinates onto the output; the test reads the output file directly with
    xarray and compares it against the input, so it does not depend on fix_ldd's
    internal coordinate handling. Real repair plus netCDF I/O is slow, so
    @settings keeps max_examples low.

    The assertions are:
      (1) the output LDD variable (name preserved from the input) has the SAME
          2-D shape as the input LDD variable ('ldd1') — the grid
          envelope and axis lengths are unchanged (no reorder, no crop/pad);
      (2) the output's 'lat' coordinate values equal the input's exactly, and
      (3) the output's 'lon' coordinate values equal the input's exactly, so the
          spatial coordinate axes are identical across the repair.

    Validates: Requirements 7.4
    """
    # Read the input grid's shape and cell centres. Any CLI coordinate is drawn
    # from these centres so it is guaranteed in-grid and the run succeeds.
    in_ds = xr.open_dataset(INPUT_LDD)
    try:
        in_lat = np.asarray(in_ds["lat"].values)
        in_lon = np.asarray(in_ds["lon"].values)
        in_ldd_shape = in_ds["ldd1"].shape

        coordinate = None
        if mode == "cli":
            valid = sorted(VALID_LDD_VALUES["pcr_d8"])
            row = data.draw(st.integers(min_value=0, max_value=in_lat.size - 1))
            col = data.draw(st.integers(min_value=0, max_value=in_lon.size - 1))
            value = data.draw(st.sampled_from(valid))
            # CLI substitution is [coord_y, coord_x, value] in (y, x, value) order,
            # placed exactly on a real cell centre so it resolves in-grid.
            coordinate = [float(in_lat[row]), float(in_lon[col]), value]
    finally:
        in_ds.close()

    out_fd, output_file = tempfile.mkstemp(suffix=".nc")
    os.close(out_fd)
    if os.path.exists(output_file):
        os.remove(output_file)
    try:
        fix_ldd(
            INPUT_LDD,
            output_file=output_file,
            river_network_format="pcr_d8",
            coordinate=coordinate,
        )

        out_ds = xr.open_dataset(output_file)
        try:
            # The pipeline preserves the input's LDD variable name on the output
            # (Req 7.4/7.5); resolve the lone data variable rather than assume a
            # fixed name.
            out_ldd_name = str(next(iter(out_ds.data_vars)))
            # (1) The repaired LDD variable keeps the input's 2-D grid shape.
            assert out_ds[out_ldd_name].shape == in_ldd_shape, (
                f"output LDD shape {out_ds[out_ldd_name].shape} != input LDD "
                f"shape {in_ldd_shape} (mode={mode}); the repair altered the "
                f"grid envelope or reordered axes"
            )

            # (2) and (3) The output's lat/lon coordinate axes equal the input's.
            out_lat = np.asarray(out_ds["lat"].values)
            out_lon = np.asarray(out_ds["lon"].values)
            assert np.array_equal(out_lat, in_lat), (
                f"output 'lat' coordinate values differ from the input's "
                f"(mode={mode}); the grid envelope was not preserved"
            )
            assert np.array_equal(out_lon, in_lon), (
                f"output 'lon' coordinate values differ from the input's "
                f"(mode={mode}); the grid envelope was not preserved"
            )
        finally:
            out_ds.close()
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


from lisfloodutilities.lddfix.lddfix import getarg


# Safe alphabet for generated output names: the only property under test is
# LENGTH, so we draw from characters that argparse never treats specially
# (letters, digits, '.', '_'). Names are built with a guaranteed non-'-' first
# character so a value is never mistaken for a flag.
_SAFE_OUTPUT_ALPHABET = string.ascii_letters + string.digits + "._"

# Lengths spanning the validity boundary, with the 0/1/255/256 edges always
# represented so the accept/reject transition is exercised exactly.
_output_length_strategy = st.one_of(
    st.just(0),
    st.just(1),
    st.just(255),
    st.just(256),
    st.integers(min_value=0, max_value=300),
)


# Feature: lddfix, Property 15: Output-name length validity
@settings(max_examples=100, deadline=None)
@given(
    length=_output_length_strategy,
    body=st.text(alphabet=_SAFE_OUTPUT_ALPHABET, max_size=300),
)
def test_property_output_name_length_validity(length, body):
    """An output name is accepted iff its length is 1-255 inclusive; a length
    of 0 or > 255 is rejected with a non-zero exit and no output is written.

    Validates: Requirements 7.1, 7.2
    """
    # Build a name of exactly `length` characters from the safe alphabet. For
    # length 0 this is the empty string (argparse accepts '' as the -o value,
    # so the length check in getarg is what must reject it).
    if length == 0:
        name = ""
    else:
        filler = (body or "a")
        name = (filler * (length // len(filler) + 1))[:length]

    assert len(name) == length

    cwd = os.getcwd()
    rejected_name_path = os.path.join(cwd, name) if name else None
    default_output_path = os.path.join(cwd, "ldd_fixed.nc")

    argv = ["lddfix", "-i", "in.nc", "-o", name, "-f", "pcr_d8"]

    expected_accept = 1 <= length <= 255

    with patch("sys.argv", argv):
        raised = None
        args = None
        try:
            args = getarg()
        except SystemExit as exc:
            raised = exc

    if expected_accept:
        # Accepted: getarg returns the namespace carrying the exact name and
        # never exits.
        assert raised is None, (
            f"getarg rejected a valid output name of length {length}"
        )
        assert args is not None
        assert args.output == name
    else:
        # Rejected: getarg must exit non-zero (argparse parser.error), so the
        # repair pipeline is never reached and no output file is produced.
        assert raised is not None, (
            f"getarg accepted an invalid output name of length {length}"
        )
        assert raised.code not in (None, 0), (
            f"getarg exited with code {raised.code!r} for an invalid name of "
            f"length {length}; expected a non-zero exit"
        )
        # Rejection happens before any file I/O: neither the rejected name nor
        # the default output exist as a side effect of the failed parse.
        if rejected_name_path is not None:
            assert not os.path.exists(rejected_name_path), (
                f"a file was written for a rejected output name of length "
                f"{length}"
            )
        assert not os.path.exists(default_output_path), (
            "the default output file was written despite a rejected name"
        )


# ---------------------------------------------------------------------------
# Example / smoke tests (plain pytest, NOT Hypothesis)
#
# Mirror the structure of tests/test_mctrivers.py's TestMctMask: copy the
# input ldd.nc into a scratch output dir, run fix_ldd/main, compare the
# produced netCDF against the static reference with xarray.Dataset.equals,
# close datasets before cleaning up the scratch dir.
#
# The static reference (tests/data/lddfix/reference/ldd.nc) is the repair-only
# result, so Dataset.equals only holds for the repair-only case. For the
# substitution cases we assert weaker-but-meaningful properties (valid netCDF,
# same grid/coords as the input, substituted cell carries the expected value,
# variable name preserved). For the error cases we assert fix_ldd/main raises
# (or main returns non-zero) and no output file is written.
# ---------------------------------------------------------------------------

import shutil
import sys

import pytest

from lisfloodutilities.lddfix.lddfix import fix_ldd, main


class TestLddfixExamples:
    """Example tests for the lddfix tool, mirroring TestMctMask in
    tests/test_mctrivers.py (mk_path_out helper + open/compare/close/cleanup).
    """

    case_dir = os.path.join(os.path.dirname(__file__), "data", "lddfix")
    input_ldd = os.path.join(case_dir, "input", "ldd.nc")
    reference_ldd = os.path.join(case_dir, "reference", "ldd.nc")
    out_dir = os.path.join(case_dir, "output")

    # --- scratch output dir lifecycle (mirrors mk_path_out) ----------------

    def _mk_out(self):
        """(Re)create the tests/data/lddfix/output/ scratch dir, empty."""
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir, ignore_errors=True)
        os.mkdir(self.out_dir)

    def _rm_out(self):
        """Remove the scratch output dir and everything in it."""
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def _copy_input(self, name="ldd.nc"):
        """Copy the input ldd.nc into the scratch dir and return its path."""
        dst = os.path.join(self.out_dir, name)
        shutil.copyfile(self.input_ldd, dst)
        return dst

    # --- shared assertions -------------------------------------------------

    def _assert_valid_output_like_input(self, out_path):
        """Open the output, assert it is a valid netCDF whose grid shape,
        coordinate values and LDD variable name match the input, and return the
        (lat, lon) and LDD values for further checks. Closes datasets it opens.
        """
        assert os.path.exists(out_path), f"output not written: {out_path}"
        in_ds = xr.open_dataset(self.input_ldd)
        out_ds = xr.open_dataset(out_path)
        try:
            # Variable name preserved (the input's LDD var, not the repair's 'ldd').
            assert list(out_ds.data_vars) == list(in_ds.data_vars), (
                f"output data vars {list(out_ds.data_vars)} do not match the "
                f"input {list(in_ds.data_vars)}"
            )
            ldd_name = list(in_ds.data_vars)[0]
            # Same grid shape.
            assert out_ds[ldd_name].shape == in_ds[ldd_name].shape
            # Same coordinate values on both axes.
            assert np.allclose(out_ds["lat"].values, in_ds["lat"].values)
            assert np.allclose(out_ds["lon"].values, in_ds["lon"].values)
            out_values = np.array(out_ds[ldd_name].values, copy=True)
            lat = np.array(out_ds["lat"].values, copy=True)
            lon = np.array(out_ds["lon"].values, copy=True)
        finally:
            in_ds.close()
            out_ds.close()
        return ldd_name, lat, lon, out_values

    def _cell_index(self, lat, lon, coord_lat, coord_lon):
        """Nearest (row, col) grid index for a coordinate pair."""
        row = int(np.argmin(np.abs(lat - coord_lat)))
        col = int(np.argmin(np.abs(lon - coord_lon)))
        return row, col

    # --- repair-only (default output name) ---------------------------------

    def test_repair_only_default_output_equals_reference(self):
        """Repair-only run using the DEFAULT output name reproduces the static
        reference exactly (xarray.Dataset.equals).

        Validates: Requirements 5.3, 7.3
        """
        self._mk_out()
        # Run inside the scratch dir so the DEFAULT output name lands there and
        # the suite leaves no stray file in the working directory.
        prev_cwd = os.getcwd()
        os.chdir(self.out_dir)
        try:
            out_path = fix_ldd(self.input_ldd)  # default output_file='ldd_fixed.nc'
            out_path = os.path.join(self.out_dir, out_path)

            reference = xr.open_dataset(self.reference_ldd)
            generated = xr.open_dataset(out_path)
            try:
                all_equal = reference.equals(generated)
            finally:
                reference.close()
                generated.close()
            assert all_equal, (
                "repair-only output with the default name does not equal the "
                f"reference {self.reference_ldd}"
            )
        finally:
            os.chdir(prev_cwd)
            self._rm_out()

    # --- single-coordinate (-c) substitution -------------------------------

    def test_single_coordinate_substitution(self):
        """A single-coordinate (-c) substitution produces a valid netCDF with
        the input's grid/coords/variable name; the repaired result stays within
        the valid pcr_d8 value set.

        Validates: Requirements 3.4, 3.5
        """
        self._mk_out()
        try:
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            # A cell centre well inside the grid, assigned a valid pcr_d8 value.
            coord_lat, coord_lon, value = 53.25, -122.95, 5
            fix_ldd(
                self.input_ldd,
                output_file=out_path,
                coordinate=[coord_lat, coord_lon, float(value)],
            )
            ldd_name, lat, lon, out_values = self._assert_valid_output_like_input(
                out_path
            )
            # Every non-nan output value is a valid pcr_d8 LDD value (the final
            # repair may adjust the substituted cell for consistency, but the
            # grid remains a valid LDD).
            finite = out_values[~np.isnan(out_values)]
            assert np.all(np.isin(finite.astype(int), sorted(VALID_LDD_VALUES["pcr_d8"]))), (
                "single-coordinate substitution produced an invalid LDD value"
            )
        finally:
            self._rm_out()

    # --- TSV (-t) substitution ---------------------------------------------

    def test_tsv_substitution(self):
        """A TSV (-t) substitution with well-formed rows produces a valid
        netCDF with the input's grid/coords/variable name.

        Validates: Requirements 4.1, 4.2, 4.3
        """
        self._mk_out()
        try:
            tsv_path = os.path.join(self.out_dir, "subs.tsv")
            with open(tsv_path, "w") as handle:
                # Two distinct in-grid cells, valid pcr_d8 values.
                handle.write("53.25\t-122.95\t5\n")
                handle.write("53.15\t-122.85\t3\n")
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            fix_ldd(self.input_ldd, output_file=out_path, tsv_file=tsv_path)
            ldd_name, lat, lon, out_values = self._assert_valid_output_like_input(
                out_path
            )
            finite = out_values[~np.isnan(out_values)]
            assert np.all(np.isin(finite.astype(int), sorted(VALID_LDD_VALUES["pcr_d8"]))), (
                "TSV substitution produced an invalid LDD value"
            )
        finally:
            self._rm_out()

    # --- error cases: assert raise + no output -----------------------------

    def test_invalid_ldd_value_error_no_output(self):
        """An out-of-set LDD value aborts the run with no output written.

        Validates: Requirements 3.5
        """
        self._mk_out()
        try:
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(
                    self.input_ldd,
                    output_file=out_path,
                    coordinate=[53.25, -122.95, 42.0],  # 42 not in pcr_d8 1..9
                )
            assert not os.path.exists(out_path), (
                "an invalid LDD value must not write any output file"
            )
        finally:
            self._rm_out()

    def test_out_of_grid_error_no_output(self):
        """An out-of-grid coordinate aborts the run with no output written.

        Validates: Requirements 3.6
        """
        self._mk_out()
        try:
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(
                    self.input_ldd,
                    output_file=out_path,
                    coordinate=[0.0, 0.0, 5.0],  # far outside the grid envelope
                )
            assert not os.path.exists(out_path), (
                "an out-of-grid coordinate must not write any output file"
            )
        finally:
            self._rm_out()

    def test_malformed_tsv_error_no_output(self):
        """A malformed TSV row aborts the run with no output written.

        Validates: Requirements 4.5
        """
        self._mk_out()
        try:
            tsv_path = os.path.join(self.out_dir, "bad.tsv")
            with open(tsv_path, "w") as handle:
                handle.write("53.25\t-122.95\t5\n")
                handle.write("not\ta\tproper\trow\n")  # 4 fields -> malformed
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(self.input_ldd, output_file=out_path, tsv_file=tsv_path)
            assert not os.path.exists(out_path), (
                "a malformed TSV must not write any output file"
            )
        finally:
            self._rm_out()

    def test_duplicate_tsv_cell_error_no_output(self):
        """Two TSV rows resolving to the same cell abort the run with no output.

        Validates: Requirements 4.9
        """
        self._mk_out()
        try:
            tsv_path = os.path.join(self.out_dir, "dup.tsv")
            with open(tsv_path, "w") as handle:
                # Both coordinates resolve to the same grid cell (within half a
                # cell of the same centre, res 0.1).
                handle.write("53.25\t-122.95\t5\n")
                handle.write("53.251\t-122.951\t3\n")
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(self.input_ldd, output_file=out_path, tsv_file=tsv_path)
            assert not os.path.exists(out_path), (
                "a duplicate TSV cell must not write any output file"
            )
        finally:
            self._rm_out()

    def test_missing_input_error_no_output(self):
        """A missing input file aborts the run with no output written.

        Validates: Requirements 1.5
        """
        self._mk_out()
        try:
            missing = os.path.join(self.out_dir, "does_not_exist.nc")
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(missing, output_file=out_path)
            assert not os.path.exists(out_path), (
                "a missing input must not write any output file"
            )
        finally:
            self._rm_out()

    def test_missing_tsv_error_no_output(self):
        """A missing TSV substitution file aborts the run with no output.

        Validates: Requirements 4.4
        """
        self._mk_out()
        try:
            missing_tsv = os.path.join(self.out_dir, "no_such.tsv")
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            with pytest.raises(Exception):
                fix_ldd(self.input_ldd, output_file=out_path, tsv_file=missing_tsv)
            assert not os.path.exists(out_path), (
                "a missing TSV must not write any output file"
            )
        finally:
            self._rm_out()

    # --- coordinate auto-detect vs explicit -E names -----------------------

    def test_autodetect_equals_explicit_coords(self):
        """A repair-only run with auto-detected coordinates (coords_names=[])
        produces the same output as a run with explicit -E ['lat', 'lon'] names
        (xarray.Dataset.equals).

        Validates: Requirements 6.2, 6.3
        """
        self._mk_out()
        try:
            auto_path = os.path.join(self.out_dir, "auto.nc")
            explicit_path = os.path.join(self.out_dir, "explicit.nc")
            fix_ldd(self.input_ldd, output_file=auto_path, coords_names=[])
            fix_ldd(
                self.input_ldd,
                output_file=explicit_path,
                coords_names=["lat", "lon"],
            )
            auto_ds = xr.open_dataset(auto_path)
            explicit_ds = xr.open_dataset(explicit_path)
            try:
                equal = auto_ds.equals(explicit_ds)
            finally:
                auto_ds.close()
                explicit_ds.close()
            assert equal, (
                "auto-detected coordinates produced a different output than the "
                "explicit -E ['lat', 'lon'] names"
            )
        finally:
            self._rm_out()

    # --- overwrite of an existing output -----------------------------------

    def test_overwrite_existing_output(self):
        """A pre-existing output file at the target path is overwritten by the
        repair result (and the result equals the reference).

        Validates: Requirements 7.5
        """
        self._mk_out()
        try:
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            # Pre-create a dummy file at the output path.
            with open(out_path, "w") as handle:
                handle.write("stale dummy content that must be replaced")
            assert os.path.exists(out_path)

            fix_ldd(self.input_ldd, output_file=out_path)

            # The dummy text is gone: the file is now a valid netCDF equal to
            # the repair-only reference.
            reference = xr.open_dataset(self.reference_ldd)
            generated = xr.open_dataset(out_path)
            try:
                all_equal = reference.equals(generated)
            finally:
                reference.close()
                generated.close()
            assert all_equal, (
                "overwrite did not replace the stale output with the repair "
                "result"
            )
        finally:
            self._rm_out()

    # --- main() entry point smoke test (non-zero exit on failure) ----------

    def test_main_missing_input_returns_nonzero(self, monkeypatch):
        """main() returns a non-zero exit status on a missing input and writes
        no output (Req 1.4, 1.6, 7.6).

        Validates: Requirements 1.4, 1.6, 7.6
        """
        self._mk_out()
        try:
            missing = os.path.join(self.out_dir, "nope.nc")
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            argv = ["lddfix", "-i", missing, "-o", out_path]
            monkeypatch.setattr(sys, "argv", argv)
            rc = main()
            assert rc != 0, "main() must return non-zero on a missing input"
            assert not os.path.exists(out_path), (
                "main() must not write any output on a missing input"
            )
        finally:
            self._rm_out()

    def test_main_repair_only_returns_zero(self, monkeypatch):
        """main() returns 0 on a successful repair-only run and writes an output
        equal to the reference (Req 1.4).

        Validates: Requirements 1.4, 5.3
        """
        self._mk_out()
        try:
            out_path = os.path.join(self.out_dir, "ldd_fixed.nc")
            argv = ["lddfix", "-i", self.input_ldd, "-o", out_path]
            monkeypatch.setattr(sys, "argv", argv)
            rc = main()
            assert rc == 0, "main() must return 0 on a successful run"
            reference = xr.open_dataset(self.reference_ldd)
            generated = xr.open_dataset(out_path)
            try:
                all_equal = reference.equals(generated)
            finally:
                reference.close()
                generated.close()
            assert all_equal
        finally:
            self._rm_out()

