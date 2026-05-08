import xarray as xr
from typing import Tuple, Optional
import numpy as np

# Define common coordinate name patterns for automatic detection
X_COORD_NAMES = ('lon', 'x', 'rlon')
Y_COORD_NAMES = ('lat', 'y', 'rlat')

# LDD direction offsets for each direction value (1-9)
# Format: (row_offset, col_offset)
# Direction values follow PCRaster convention:
# N=8, NE=9, E=6, SE=3, S=2, SW=1, W=4, NW=7, 5=pit
LDD_OFFSETS = {
    9: (-1, 1),  # Northeast
    8: (-1, 0),  # North
    7: (-1, -1), # Northwest
    6: (0, 1),   # East
    4: (0, -1),  # West
    3: (1, 1),   # Southeast
    2: (1, 0),   # South
    1: (1, -1),  # Southwest
    5: (0, 0),   # Pit (no flow)
}

LDD_OFFSET_MIN = 1
LDD_OFFSET_MAX = 9
LDD_PIT_VALUE = 5  # LDD value for pits (no flow)
LDD_MISSING_VALUE = 0  # Value to use for missing/invalid LDD cells in the repaired LDD
VALID_LDD_VALUES = set(LDD_OFFSETS.keys()).difference({LDD_PIT_VALUE})  # Valid LDD values are 1-9 except 5 which is a pit


def copy_coordinates_and_attributes(source: xr.Dataset, target: xr.DataArray, y_proj: str, x_proj: str) -> xr.DataArray:
    """
    Copy coordinates and attributes from source Dataset to target DataArray.
    
    Parameters:
    -----------
    source : xr.Dataset
        The Dataset from which to copy attributes.
    target : xr.DataArray
        The DataArray to which attributes will be copied.
    
    Returns:
    --------
    xr.DataArray
        The target DataArray with copied attributes.
    """
    # target.attrs = source.attrs.copy()
    target = target.assign_coords({y_proj: source.coords[y_proj], x_proj: source.coords[x_proj]})
    
    # Get the main variable from the source Dataset (first data variable)
    if len(source.data_vars) > 0:
        main_var_name = list(source.data_vars)[0]
        main_var = source[main_var_name]
        
        # Copy attributes from the main variable (including grid_mapping, esri_pe_string, etc.)
        for attr_name in ['esri_pe_string', 'spatial_ref', 'crs_wkt', 'grid_mapping']:
            if attr_name in main_var.attrs:
                target.attrs[attr_name] = main_var.attrs[attr_name]
    
    # Also check for grid_mapping coordinate and copy if present
    if 'grid_mapping' in source.coords:
        target = target.assign_coords({'grid_mapping': source.coords['grid_mapping']})
        if 'grid_mapping' not in target.attrs:
            target.attrs['grid_mapping'] = source.coords['grid_mapping'].values
    return target


def get_downstream(i: int, j: int, ldd_arr: np.ndarray, rows: int, cols: int) -> Optional[Tuple[int, int]]:
    """
    Get the downstream cell coordinates for a given cell.
    Parameters:
    -----------
    i : int
        Row index of the cell.
    j : int
        Column index of the cell.
    ldd_arr : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    rows : int
        Number of rows in the LDD array.
    cols : int
        Number of columns in the LDD array.

    Returns:
    --------
    Optional[Tuple[int, int]]
        Coordinates of the downstream cell, or None if no valid downstream cell exists.
    """
    ldd_val = int(ldd_arr[i, j])
    if ldd_val < LDD_OFFSET_MIN or ldd_val > LDD_OFFSET_MAX or ldd_val == LDD_PIT_VALUE:
        return None  # Invalid or pit
    di, dj = LDD_OFFSETS[ldd_val]
    next_i, next_j = i + di, j + dj
    if 0 <= next_i < rows and 0 <= next_j < cols:
        return (next_i, next_j)
    return None  # Flows out of bounds


def remove_cycles(ldd_arr: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """
    Remove cycles from the LDD array by assigning missing values to all cells in a cycle.
    Parameters:
    -----------
    ldd_arr : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    rows : int
        Number of rows in the LDD array.
    cols : int
        Number of columns in the LDD array.
    Returns:
    --------
    np.ndarray
        LDD array with cycles removed (cells in cycles set to missing value).
    """
    # Track cells that are part of cycles
    cycle_cells = set()
    
    # For each cell, follow the downstream path to detect cycles
    for i in range(rows):
        for j in range(cols):
            ldd_val = int(ldd_arr[i, j])
            
            # Skip invalid cells and pits
            if ldd_val < LDD_OFFSET_MIN or ldd_val > LDD_OFFSET_MAX or ldd_val == LDD_PIT_VALUE:
                continue
            
            # Follow the downstream path, tracking visited cells
            path = []  # List of (i, j) tuples in order visited
            visited_in_path = set()  # Set of cells visited in current path
            
            current_i, current_j = i, j
            
            while current_i is not None and current_j is not None:
                # Check if we've hit a pit or invalid cell - no cycle here
                current_ldd = int(ldd_arr[current_i, current_j])
                if current_ldd < LDD_OFFSET_MIN or current_ldd > LDD_OFFSET_MAX or current_ldd == LDD_PIT_VALUE:
                    break
                
                # Check if we've already visited this cell in the current path - cycle detected!
                if (current_i, current_j) in visited_in_path:
                    # Found a cycle - mark all cells in the cycle as missing
                    # Find the start of the cycle in the path
                    cycle_start_idx = path.index((current_i, current_j))
                    for cycle_i, cycle_j in path[cycle_start_idx:]:
                        cycle_cells.add((cycle_i, cycle_j))
                    break
                
                # Add current cell to path
                path.append((current_i, current_j))
                visited_in_path.add((current_i, current_j))
                
                # Move to downstream cell
                downstream = get_downstream(current_i, current_j, ldd_arr, rows, cols)
                if downstream is None:
                    break
                current_i, current_j = downstream
    
    # Assign missing values to all cells in cycles
    for i, j in cycle_cells:
        ldd_arr[i, j] = LDD_MISSING_VALUE
    return ldd_arr


def correct_edge_pits(ldd_arr: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """
    Correct edge cells that flow out of bounds by assigning them the LDD code of a pit cell (5).
    Parameters:
    -----------
    ldd_arr : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    rows : int
        Number of rows in the LDD array.
    cols : int
        Number of columns in the LDD array.
    Returns:
    --------
    np.ndarray
        LDD array with edge cells that flow out of bounds corrected to pit code.
    """
    for i in range(rows):
        for j in range(cols):
            ldd_val = int(ldd_arr[i, j])
            
            # Skip if already a pit (5) or invalid
            if ldd_val == LDD_PIT_VALUE or ldd_val < LDD_OFFSET_MIN or ldd_val > LDD_OFFSET_MAX:
                continue
            
            # Check if this cell's immediate downstream is out of bounds
            di, dj = LDD_OFFSETS[ldd_val]
            next_i, next_j = i + di, j + dj
            
            if not (0 <= next_i < rows and 0 <= next_j < cols):
                # Flows out of bounds - convert to pit
                ldd_arr[i, j] = LDD_PIT_VALUE
    return ldd_arr


def correct_invalid_downstream(ldd_arr: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """
    Correct cells that flow to a cell with missing value (including cells in cycles) by assigning them the LDD code of a pit cell (5).
    Parameters:
    -----------
    ldd_arr : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    rows : int
        Number of rows in the LDD array.
    cols : int
        Number of columns in the LDD array.
    Returns:
    --------
    np.ndarray
        LDD array with cells that flow to invalid downstream cells corrected to pit code.
    """
    for i in range(rows):
        for j in range(cols):
            ldd_val = int(ldd_arr[i, j])
            
            # Skip if already a pit (5) or invalid
            if ldd_val == LDD_PIT_VALUE or ldd_val < LDD_OFFSET_MIN or ldd_val > LDD_OFFSET_MAX:
                continue
            
            # Check if this cell's downstream is a missing value
            di, dj = LDD_OFFSETS[ldd_val]
            next_i, next_j = i + di, j + dj
            
            if 0 <= next_i < rows and 0 <= next_j < cols:
                downstream_val = ldd_arr[next_i, next_j]
                if downstream_val == LDD_MISSING_VALUE:
                    # Drains to a missing value - convert to pit
                    ldd_arr[i, j] = LDD_PIT_VALUE
    return ldd_arr


def correct_invalid_values(ldd_arr: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """
    Correct cells with invalid LDD values (not in 1-9) by assigning them the LDD code of a pit cell (5).
    Parameters:
    -----------
    ldd_arr : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    rows : int
        Number of rows in the LDD array.
    cols : int
        Number of columns in the LDD array.
    Returns:
    --------
    np.ndarray
        LDD array with invalid values corrected to pit code.
    """
    # Find invalid values (not 1-9)
    invalid_mask = (ldd_arr < LDD_OFFSET_MIN) | (ldd_arr > LDD_OFFSET_MAX)
    for i in range(rows):
        for j in range(cols):
            if invalid_mask[i, j]:
                # Check if any neighbor flows into this cell
                is_outlet = False
                for direction_value in VALID_LDD_VALUES:
                    (di, dj) = LDD_OFFSETS[direction_value]
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols:
                        # Check if neighbor flows into this cell
                        # The neighbor's LDD should point to this cell
                        neighbor_ldd = ldd_arr[ni, nj]
                        if neighbor_ldd in VALID_LDD_VALUES:
                            # Check if the direction points to current cell
                            ndi, ndj = LDD_OFFSETS[neighbor_ldd]
                            if ni + ndi == i and nj + ndj == j:
                                is_outlet = True
                                break
                if is_outlet:
                    ldd_arr[i, j] = LDD_PIT_VALUE  # Set as pit
                else:
                    ldd_arr[i, j] = LDD_MISSING_VALUE  # Set as missing
    return ldd_arr


def lddrepair(ldd_array: np.ndarray) -> np.ndarray:
    """
    Repair LDD array - ensures all drainage paths end in a pit.
    Similar to pcraster's lddrepair function.
    
    The repair operation is done as follows:
    1. First, the cycles are removed by assigning missing values to all cells in a cycle.
    2. Second, cells with a local drain direction to the outside of the map or to a
       cell with a missing value (including cells that were in a cycle) are assigned
       the ldd code of a pit cell (code: 5).
    3. Third, cells with a local drain direction to a cell with a missing value
       (including cells that were in a cycle) are assigned the ldd code of a pit cell (code: 5).
    
    Parameters:
    -----------
    ldd_array : np.ndarray
        LDD array with values 1-9 (or -1/0 for missing/invalid)
    
    Returns:
    --------
    np.ndarray
        Repaired LDD array
    """
    # Create output array - start with copy of input
    ldd_repaired = ldd_array.copy()
    rows, cols = ldd_repaired.shape
    
    # Step 1: Remove cycles by assigning missing values to all cells in a cycle
    # A cycle is a set of cells that don't drain to a pit because they drain to each other
    # We detect cycles by following downstream paths and checking for visited cells
    ldd_repaired = remove_cycles(ldd_repaired, rows, cols)

    # Step 2: Find cells that are at the edge and flow out of bounds
    # These should be converted to pits
    ldd_repaired = correct_edge_pits(ldd_repaired, rows, cols)

    # Step 2b: Find cells that drain to a cell with missing value (including cells in cycles)
    # These should be converted to pits
    ldd_repaired = correct_invalid_downstream(ldd_repaired, rows, cols)
    
    # Step 3: For each invalid cell, check if any neighbor flows into it
    # If so, make this cell an outlet (pit = 5)
    ldd_repaired = correct_invalid_values(ldd_repaired, rows, cols)

    return ldd_repaired