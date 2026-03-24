__author__="Goncalo Gomes"
__date__="$Mar 24, 2026 12:01:00$"
__version__="0.1"
__updated__="$Mar 24, 2026 16:01:00$"

"""
This script downloads KIWI files from the WISKI API for a variable and time frame.
It downloads metadata, filters stations in the EFAS domain, downloads timeseries
data for each station, and merges the data with metadata to create KIWI files.

Usage:
    python download_timeseries.py <variable> [options]

Examples:
    python download_timeseries.py pr6
    python download_timeseries.py pr6 --start 2024-12-31 --end 2026-01-02

Copyright 2019-2020 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.
"""

import argparse
import os
import sys
import time
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import urllib.request
import urllib.error
import urllib.parse
import base64
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from lisfloodutilities.gridding.lib.utils import Config, FileUtils

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Columns for metadata file
METADATA_MAX_FIELDS_KEY = 'max_fields'
METADATA_COL_VALUE = 'ts_value'
METADATA_COL_QCODE = 'q_code'
METADATA_COL_STATION_ID = 'station_id'
METADATA_COL_NOGRIDDING = 'EFAS-ADDATTR-NOGRIDDING'
METADATA_COL_ISINARCMINDOMAIN = 'EFAS-ADDATTR-ISINARCMINDOMAIN'

metadata_indices = {}

# Skip header (first 9 lines based on the bash script)
# The data starts after line 9 (tail -n +10 in bash)
TIMESERIES_FILE_HEADER_LINES = 9

# TIMESERIES_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S"
TIMESERIES_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S.000Z"

# Column indices for timeseries files (0-based)
MAX_TIMESERIES_FIELDS = 3
TIMESERIES_IDX_TIMESTAMP = 0
TIMESERIES_IDX_VALUE = 1
TIMESERIES_IDX_QCODE = 2

# Constants
NEWLINE = '\n'
COLUMN_SEPARATOR = '\t'

# Date format for command-line arguments
ARG_DATE_FORMAT = "%Y-%m-%d"

# Environment variable for API key
ENV_KIWI_API_KEY = "KIWI_API_KEY"

# API access KEY (can be set via environment variable or command-line argument)
DEFAULT_KEY = os.getenv(ENV_KIWI_API_KEY, "")


class DownloadError(Exception):
    """Custom exception for download errors."""
    pass


class APIError(Exception):
    """Custom exception for API errors."""
    pass


def calculate_metadata_indices(header_line: str) -> dict:
    """
    Calculate the column indices from the metadata file header.
    
    This function parses the header line of a metadata file and finds the
    indices of the required columns based on their names.
    If the indices have already been calculated and stored in the
    global variable `metadata_indices`, it returns that instead of recalculating.
    
    Args:
        header_line: The header line from the metadata file (tab-separated)
    
    Returns:
        A dictionary with the following keys:
        - METADATA_COL_VALUE: index of the ts_value column
        - METADATA_COL_QCODE: index of the q_code column
        - METADATA_COL_STATION_ID: index of the station_id column
        - METADATA_COL_NOGRIDDING: index of the EFAS-ADDATTR-NOGRIDDING column
        - METADATA_COL_ISINARCMINDOMAIN: index of the EFAS-ADDATTR-ISINARCMINDOMAIN column
        - METADATA_MAX_FIELDS_KEY: total number of columns in the header
    """
    global metadata_indices
    if len(metadata_indices) > 0:
        return metadata_indices

    columns = header_line.strip().split(COLUMN_SEPARATOR)

    indices = {
        METADATA_COL_VALUE: columns.index(METADATA_COL_VALUE),
        METADATA_COL_QCODE: columns.index(METADATA_COL_QCODE),
        METADATA_COL_STATION_ID: columns.index(METADATA_COL_STATION_ID),
        METADATA_COL_NOGRIDDING: columns.index(METADATA_COL_NOGRIDDING),
        METADATA_COL_ISINARCMINDOMAIN: columns.index(METADATA_COL_ISINARCMINDOMAIN),
        METADATA_MAX_FIELDS_KEY: len(columns)
    }
    metadata_indices.update(indices)

    return metadata_indices


def create_auth_header(key: str) -> str:
    """Create the Authorization header for API requests."""
    # encoded = base64.b64encode(key.encode()).decode()
    encoded = key
    return f"Basic {encoded}"


def build_url(base_url: str, path: str, params: Dict[str, str]) -> str:
    """Build a URL with query parameters."""
    # Use safe=',' to keep commas unencoded in query string values
    query_string = urllib.parse.urlencode(params, safe=',: */.-')

    return f"https://{base_url}/{path}?{query_string}"


def download_file(url: str, output_path: str, auth_header: str,
                  max_retries: int = 5, retry_delay: int = 20) -> bool:
    """
    Download a file from URL to output_path with retry logic.
    
    Args:
        url: The URL to download from
        output_path: Local path to save the file
        auth_header: Authorization header
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retries in seconds
        
    Returns:
        True if download was successful, False otherwise
    """
    headers = {
        "Authorization": auth_header
    }
    
    # Get system proxy settings
    proxies = urllib.request.getproxies()

    # Build the opener with proxy support if proxies are configured
    if proxies:
        proxy_handler = urllib.request.ProxyHandler(proxies)
        opener = urllib.request.build_opener(proxy_handler)
    else:
        opener = urllib.request.build_opener()
    
    for attempt in range(max_retries):
        try:
            request = urllib.request.Request(url, headers=headers)
            
            # Use the opener which respects system proxy settings
            with opener.open(request, timeout=300) as response:
                content = response.read().decode('utf-8')
                
                # Check for HTML error responses
                if content.strip().startswith('<html') or '<html>' in content:
                    raise APIError(f"API returned HTML error response: {content[:500]}")
                
                # Write content to file
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(content)
                
                logger.info(f"Successfully downloaded: {output_path}")
                return True
                
        except urllib.error.HTTPError as e:
            # Try to read the error response body for more details
            error_body = ""
            try:
                error_body = e.read().decode('utf-8')
            except:
                pass
            logger.warning(f"HTTP Error {e.code} on attempt {attempt + 1}/{max_retries}: {e.reason}")
            if error_body:
                logger.warning(f"Error response body: {error_body[:500]}")
            if attempt < max_retries - 1:
                logger.info(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                logger.error(f"Failed to download after {max_retries} attempts: {url}")
                raise DownloadError(f"HTTP Error {e.code}: {e.reason}. Response: {error_body[:500]}") from e
                
        except urllib.error.URLError as e:
            logger.warning(f"URL Error on attempt {attempt + 1}/{max_retries}: {e.reason}")
            if attempt < max_retries - 1:
                logger.info(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                logger.error(f"Failed to download after {max_retries} attempts: {url}")
                raise DownloadError(f"URL Error: {e.reason}") from e
                
        except APIError as e:
            logger.error(f"API Error: {e}")
            raise
            
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise DownloadError(f"Unexpected error: {e}") from e
    
    return False


def check_file_for_errors(file_path: str) -> bool:
    """
    Check if a file contains HTML error content.
    
    Args:
        file_path: Path to the file to check
        
    Returns:
        True if file contains errors, False otherwise
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read(500)  # Read first 500 chars to check for HTML
            return '<html>' in content.lower() or content.strip().startswith('<!DOCTYPE')
    except Exception as e:
        logger.warning(f"Error checking file {file_path}: {e}")
        return True


def download_metadata(conf: Config, variable: str, base_path: str, yyyy: str, mm: str, dd: str,
                      auth_header: str, do_download: bool = True) -> str:
    """
    Download the metadata file from the API.
    
    Args:
        conf: Configuration object
        variable: The variable to download (e.g., 'pr6', 'ta6')
        base_path: Base path for output files
        yyyy: Year
        mm: Month
        dd: Day
        auth_header: Authorization header
        do_download: Whether to actually download the file
        
    Returns:
        Path to the metadata file
    """
    timeseries_id = conf.api_timeseries_id
    if timeseries_id is None or len(timeseries_id) == 0:
        logger.error(f"No timeseries ID configured for variable '{variable}'")
        raise ValueError(f"No timeseries ID configured for variable '{variable}'")
    hh = conf.api_timeseries_time_to_download
    if hh is None or len(hh) == 0:
        logger.error(f"No timeseries time to download configured for variable '{variable}'")
        raise ValueError(f"No timeseries time to download configured for variable '{variable}'")
    
    stations_metadata_file = os.path.join(base_path, f"{variable}_timeseries_metadata.tsv")
    
    if not do_download:
        return stations_metadata_file
    
    # Build the metadata URL
    params = {
        "datasource": "1",
        "service": "kisters",
        "type": "queryServices",
        "request": "getTimeseriesValuelayer",
        "format": "ascii",
        "timeseriesgroup_id": timeseries_id,
        "date": f"{yyyy}-{mm}-{dd}T{hh}:00:00",
        "metadata": "true",
        "invalidValue": "-999",
        "invalidPeriod": "PT1H",
        "md_returnfields": f"site_no,site_name,station_name,station_no,{METADATA_COL_STATION_ID},ts_unitsymbol,ts_shortname,parametertype_name,ca_sta,ca_par",
        "crs": "global",
        "ca_sta_returnfields": f"station_status,station_diary_status,{METADATA_COL_NOGRIDDING},EFAS-ADDATTR-COUNTRY,{METADATA_COL_ISINARCMINDOMAIN}",
        "ca_par_returnfields": "EXCLUDE,INACTIVE_histattr",
        "returnfields": f"sta_location,{METADATA_COL_VALUE},{METADATA_COL_QCODE}",
        "valueType": "matchingValue_nocalc"
    }
    
    url = build_url(conf.api_url, "EFASDATA/KiWIS", params)
    
    logger.info(f"DOWNLOADING METADATA FILE: {stations_metadata_file}")
    download_file(url, stations_metadata_file, auth_header)
    logger.info(f"DOWNLOADED METADATA FILE: {stations_metadata_file}")
    
    # Process the metadata file - replace value and qcode columns with placeholders
    process_metadata_file(stations_metadata_file)
    
    return stations_metadata_file


def process_metadata_file(metadata_file: str) -> None:
    """
    Process the metadata file to replace value and qcode with placeholders.
    
    Args:
        metadata_file: Path to the metadata file
    """
    try:
        with open(metadata_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        if not lines:
            return
            
        # Process each line (skip header)
        processed_lines = []
        indices = None
        for i, line in enumerate(lines):
            fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            if i == 0:
                # Header line - keep as is and calculate indices
                processed_lines.append(line.rstrip(NEWLINE))
                indices = calculate_metadata_indices(line)
            elif indices is not None:
                # Data lines - replace value and qcode columns with placeholders
                if len(fields) >= indices[METADATA_MAX_FIELDS_KEY]:
                    fields[indices[METADATA_COL_VALUE]] = "{value}"
                    fields[indices[METADATA_COL_QCODE]] = "{qcode}"
                processed_lines.append(COLUMN_SEPARATOR.join(fields))
            else:
                # No header found, just append the line as-is
                processed_lines.append(line.rstrip(NEWLINE))

        # Write back
        with open(metadata_file, 'w', encoding='utf-8') as f:
            f.write(NEWLINE.join(processed_lines) + NEWLINE)
            
        logger.info(f"Processed metadata file: {metadata_file}")
        
    except Exception as e:
        logger.error(f"Error processing metadata file: {e}")
        raise


def download_station_list(conf: Config, variable: str, base_path: str, auth_header: str,
                          do_download: bool = True) -> str:
    """
    Download the list of stations from the API.
    
    Args:
        conf: Configuration object
        variable: The variable to download
        base_path: Base path for output files
        auth_header: Authorization header
        do_download: Whether to actually download the file
        
    Returns:
        Path to the stations file
    """
    stations_filename = os.path.join(base_path, f"{variable}_stations_to_download.tsv")
    
    if not do_download:
        return stations_filename
    
    ts_path_wildcard = conf.api_timeseries_path_wildcard
    if ts_path_wildcard is None or len(ts_path_wildcard) == 0:
        logger.error(f"No timeseries path wildcard configured for variable '{variable}'")
        raise ValueError(f"No timeseries path wildcard configured for variable '{variable}'")
    
    # Build the station list URL
    params = {
        "datasource": "1",
        "service": "kisters",
        "type": "queryServices",
        "request": "getTimeseriesList",
        "format": "ascii",
        "ts_path": ts_path_wildcard,
        "returnfields": "station_id,ts_path"
    }
    
    url = build_url(conf.api_url, "EFASDATA/KiWIS", params)
    
    stations_to_download_file = os.path.join(base_path, f"{variable}_stations_to_download_raw.tsv")
    
    logger.info(f"DOWNLOADING STATIONS FILE: {stations_to_download_file}")
    download_file(url, stations_to_download_file, auth_header)
    logger.info(f"DOWNLOADED STATIONS FILE: {stations_to_download_file}")
    
    # Filter stations in EFAS domain
    stations_metadata_file = os.path.join(base_path, f"{variable}_timeseries_metadata.tsv")
    station_ids_efas_domain_file = os.path.join(base_path, f"{variable}_station_ids_efas_domain.tsv")
    
    filter_stations_efas_domain(stations_metadata_file, stations_to_download_file, 
                                 station_ids_efas_domain_file, stations_filename)
    
    return stations_filename


def filter_stations_efas_domain(metadata_file: str, stations_file: str,
                                 efas_domain_file: str, output_file: str) -> None:
    """
    Filter stations to only include those in the EFAS domain.
    
    Args:
        metadata_file: Path to the metadata file
        stations_file: Path to the raw stations file
        efas_domain_file: Path to output EFAS domain station IDs
        output_file: Path to output filtered stations
    """
    try:
        # Get station_ids in EFAS domain (column 15 != 'yes' and column 17 != 'no')
        efas_station_ids = set()
        with open(metadata_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        # Calculate column indices from header
        if not lines:
            logger.warning(f"No lines found in metadata file: {metadata_file}")
            return
        indices = calculate_metadata_indices(lines[0])
        for line in lines[1:]:  # Skip header
            fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            if len(fields) >= indices[METADATA_MAX_FIELDS_KEY]:
                nogridding = fields[indices[METADATA_COL_NOGRIDDING]] if len(fields) > indices[METADATA_COL_NOGRIDDING] else ""
                isinarcmindomain = fields[indices[METADATA_COL_ISINARCMINDOMAIN]] if len(fields) > indices[METADATA_COL_ISINARCMINDOMAIN] else ""
                
                if nogridding != "yes" and isinarcmindomain != "no":
                    if len(fields) > indices[METADATA_COL_STATION_ID]:
                        efas_station_ids.add(fields[indices[METADATA_COL_STATION_ID]])
        
        # Write EFAS domain station IDs
        with open(efas_domain_file, 'w', encoding='utf-8') as f:
            for station_id in efas_station_ids:
                f.write(f"{station_id}{NEWLINE}")
        
        logger.info(f"Found {len(efas_station_ids)} stations in EFAS domain")
        
        # Filter stations file
        with open(stations_file, 'r', encoding='utf-8') as f:
            station_lines = f.readlines()
        
        filtered_stations = []
        for line in station_lines[1:]:  # Skip header
            fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            if len(fields) >= 2 and fields[0] in efas_station_ids:
                filtered_stations.append(line.rstrip(NEWLINE))
        
        # Write filtered stations
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"station_id{COLUMN_SEPARATOR}ts_path{NEWLINE}")  # Header
            for station in filtered_stations:
                f.write(f"{station}{NEWLINE}")
        
        logger.info(f"Filtered to {len(filtered_stations)} stations to download")
        
    except Exception as e:
        logger.error(f"Error filtering stations: {e}")
        raise


def download_station_data(conf: Config, variable: str, base_path: str, start_period: str, 
                          end_period: str, auth_header: str, 
                          do_download: bool = True) -> List[Tuple[str, str]]:
    """
    Download timeseries data for each station.
    
    Args:
        conf: Configuration object
        variable: The variable to download
        base_path: Base path for output files
        start_period: Start date (YYYY-MM-DD)
        end_period: End date (YYYY-MM-DD)
        auth_header: Authorization header
        do_download: Whether to actually download the files
        
    Returns:
        List of (station_id, file_path) tuples
    """
    base_folder_timeseries = os.path.join(base_path, "timeseries", variable)
    os.makedirs(base_folder_timeseries, exist_ok=True)
    
    stations_filename = os.path.join(base_path, f"{variable}_stations_to_download.tsv")
    
    if not do_download:
        # Return list of existing files
        files = []
        if os.path.exists(base_folder_timeseries):
            for f in os.listdir(base_folder_timeseries):
                if f.endswith('_timeseries.tsv'):
                    # Extract station_id from filename
                    match = re.search(r'station_(\d+)_timeseries', f)
                    if match:
                        files.append((match.group(1), os.path.join(base_folder_timeseries, f)))
        return files
    
    # Read stations to download
    with open(stations_filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    stations = []
    for line in lines[1:]:  # Skip header
        fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
        if len(fields) >= 2:
            stations.append((fields[0], fields[1]))
    
    logger.info(f"Downloading data for {len(stations)} stations...")
    
    downloaded_files = []
    errors = []
    
    for station_id, ts_path in stations:
        output = os.path.join(base_folder_timeseries, 
                              f"{variable}_station_{station_id}_timeseries.tsv")
        
        # Build the data URL
        params = {
            "datasource": "1",
            "service": "kisters",
            "type": "queryServices",
            "request": "getTimeseriesValues",
            "format": "ascii",
            "ts_path": ts_path,
            "from": start_period,
            "to": end_period,
            "metadata": "true",
            "returnfields": "Timestamp,Value,Quality Code",
            "md_returnfields": "station_id,site_no,site_name,station_no,station_name,station_latitude,station_longitude"
        }
        
        url = build_url(conf.api_url, "EFASDATA/KiWIS", params)
        
        logger.info(f"DOWNLOADING FILE: {output}")
        
        try:
            download_file(url, output, auth_header)
            
            # Check for errors in the downloaded file
            if check_file_for_errors(output):
                logger.error(f"File contains error: {output}")
                errors.append(output)
            else:
                downloaded_files.append((station_id, output))
                
        except DownloadError as e:
            logger.error(f"Failed to download {output}: {e}")
            errors.append(output)
    
    # Report files with errors
    if errors:
        logger.warning(f"Detected {len(errors)} files with errors:")
        for f in errors:
            logger.warning(f"  File with error: {f}")
    
    return downloaded_files


def merge_timeseries_with_metadata(conf: Config, variable: str, base_path: str, 
                                    start_period: str, end_period: str,
                                    do_merge: bool = True) -> None:
    """
    Merge all timeseries with the metadata to create KIWI files.
    
    Args:
        conf: Configuration object
        variable: The variable to process
        base_path: Base path for files
        start_period: Start date (YYYY-MM-DD)
        end_period: End date (YYYY-MM-DD)
        do_merge: Whether to actually perform the merge
    """
    if not do_merge:
        return
    
    logger.info("MERGE TIMESERIES INTO KIWI FILES")
    
    base_folder_timeseries = os.path.join(base_path, "timeseries", variable)
    base_folder_meteo = os.path.join(base_path, "meteo", variable)
    
    start_year = int(start_period[:4])
    end_year = int(end_period[:4])
    
    stations_metadata_file = os.path.join(base_path, f"{variable}_timeseries_metadata.tsv")
    
    # Read metadata into a dictionary keyed by station_id
    metadata_dict = {}
    with open(stations_metadata_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header = lines[0].rstrip(NEWLINE)
    indices = calculate_metadata_indices(header)
    
    for line in lines[1:]:
        fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
        if len(fields) > indices[METADATA_COL_STATION_ID]:
            station_id = fields[indices[METADATA_COL_STATION_ID]]
            metadata_dict[station_id] = line.rstrip(NEWLINE)
    
    # Process each timeseries file
    timeseries_files = []
    if os.path.exists(base_folder_timeseries):
        for f in os.listdir(base_folder_timeseries):
            if f.endswith('_timeseries.tsv'):
                timeseries_files.append(os.path.join(base_folder_timeseries, f))
    
    for cur_filename in sorted(timeseries_files):
        base_name = os.path.basename(cur_filename)
        # Extract station_id from filename: pr6_station_12345_timeseries.tsv
        match = re.search(r'station_(\d+)_timeseries', base_name)
        if not match:
            continue
            
        station_id = match.group(1)
        metadata_row = metadata_dict.get(station_id, "")
        
        if not metadata_row:
            logger.warning(f"No metadata found for station {station_id}")
            continue
        
        logger.info(f"Processing file: {cur_filename}")
        
        # Read the timeseries data
        try:
            with open(cur_filename, 'r', encoding='utf-8') as f:
                data_lines = f.readlines()
            
            # Skip header (first 9 lines) data starts from line 10 (index 9)
            data_start = TIMESERIES_FILE_HEADER_LINES
            
            for line in data_lines[data_start:]:
                fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)

                if len(fields) < MAX_TIMESERIES_FIELDS:
                    continue

                cur_timestep = fields[TIMESERIES_IDX_TIMESTAMP]
                cur_value = fields[TIMESERIES_IDX_VALUE]
                cur_qcode = fields[TIMESERIES_IDX_QCODE]
                
                # Skip empty values
                if not cur_value or not cur_timestep:
                    continue

                cur_datetime = datetime.strptime(cur_timestep, TIMESERIES_TIMESTAMP_FORMAT)
                
                # Parse the timestamp
                cur_year = cur_datetime.strftime("%Y")
                cur_month = cur_datetime.strftime("%m")
                cur_day = cur_datetime.strftime("%d")

                if int(cur_year) < start_year or int(cur_year) > end_year:
                    continue

                # Create the output folder and file
                working_folder = os.path.join(base_folder_meteo, cur_year, cur_month, cur_day)
                os.makedirs(working_folder, exist_ok=True)
                
                kiwi_filename = cur_datetime.strftime(conf.input_timestamp_pattern)

                kiwi_file = os.path.join(working_folder, kiwi_filename)
                
                # Write header if file doesn't exist
                if not os.path.exists(kiwi_file):
                    with open(kiwi_file, 'w', encoding='utf-8') as f:
                        f.write(header + NEWLINE)
                
                # Replace placeholders and append data
                modified_row = metadata_row.replace("{value}", cur_value).replace("{qcode}", cur_qcode)
                with open(kiwi_file, 'a', encoding='utf-8') as f:
                    f.write(modified_row + NEWLINE)
        
        except Exception as e:
            logger.error(f"Error processing {cur_filename}: {e}")
            continue
    
    logger.info("FINISHED MERGING TIMESERIES")


def create_parser() -> argparse.ArgumentParser:
    """Main function to run the download script."""
    parser = argparse.ArgumentParser(
        description="Download timeseries data from WISKI API"
    )
    parser.add_argument(
        "variable",
        type=str,
        help="Variable to download (e.g., pr6, ta6, pr, tx, tn, pd, ws, rg, wx)"
    )
    parser.add_argument(
        "-c", "--conf",
        dest="config_type",
        required=True,
        help="Set the grid configuration type to use.",
        metavar="{5x5km, 1arcmin,...}"
    )
    parser.add_argument(
        "-p", "--pathconf",
        dest="config_base_path",
        required=False,
        type=FileUtils.folder_type,
        help="Overrides the base path where the configurations are stored.",
        metavar="/path/to/config"
    )
    parser.add_argument(
        "--base-path",
        type=FileUtils.folder_type,
        default="/tmp/download_timeseries",
        help="Base path for output files"
    )
    parser.add_argument(
        "--start",
        type=str,
        default="1989-12-31",
        help="Start date (YYYY-MM-DD)"
    )
    parser.add_argument(
        "--end",
        type=str,
        default=None,
        help="End date (YYYY-MM-DD). Defaults to current date."
    )
    parser.add_argument(
        "--no-metadata",
        action="store_true",
        help="Skip downloading metadata"
    )
    parser.add_argument(
        "--no-stations",
        action="store_true",
        help="Skip downloading station list"
    )
    parser.add_argument(
        "--no-data",
        action="store_true",
        help="Skip downloading station data"
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        help="Skip merging timeseries with metadata"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help=f"API key for authentication (if not provided, uses the environment variable {ENV_KIWI_API_KEY})"
    )
    return parser


def main():
    """Main function to run the download script."""
    parser = create_parser()

    args = parser.parse_args()

    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Parse start date for metadata download
        start_date = datetime.strptime(args.start, ARG_DATE_FORMAT)
        start_year = start_date.strftime("%Y")
        start_month = start_date.strftime("%m")
        start_day = start_date.strftime("%d")
    except ValueError:
        logger.error("Invalid start date format. Use YYYY-MM-DD.")
        sys.exit(1)

    # Set end date
    if args.end is None:
        # If end date is not provided, use current date
        end_date = datetime.now()
        end_year = end_date.strftime("%Y")
        end_month = end_date.strftime("%m")
        end_day = end_date.strftime("%d")
        args.end = f"{end_year}-{end_month}-{end_day}"
    else:
        # Validate end date format
        try:
            end_date = datetime.strptime(args.end, ARG_DATE_FORMAT)
            end_year = end_date.strftime("%Y")
            end_month = end_date.strftime("%m")
            end_day = end_date.strftime("%d")
        except ValueError:
            logger.error("Invalid end date format. Use YYYY-MM-DD.")
            sys.exit(1)

    # Create base path
    os.makedirs(args.base_path, exist_ok=True)

    # Create authorization header (use provided key or default)
    api_key = args.api_key if args.api_key else DEFAULT_KEY
    if len(api_key) == 0:
        logger.error("No API key provided. Set it via --api-key argument or environment variable KIWI_API_KEY.")
        sys.exit(1)

    auth_header = create_auth_header(api_key)

    # Determine what to do
    do_download_metadata = not args.no_metadata
    do_download_stations = not args.no_stations
    do_download_data = not args.no_data
    do_merge = not args.no_merge

    # Setup the configuration for the specified variable and config type
    configuration_base_folder = os.path.join(program_path, '../src/lisfloodutilities/gridding/configuration')
    if args.config_base_path is not None and len(args.config_base_path) > 0:
        configuration_base_folder = args.config_base_path
    file_utils = FileUtils(args.variable, quiet_mode=(not args.verbose))
    config_type_path = file_utils.get_config_type_path(configuration_base_folder, args.config_type)
    config_filename = file_utils.get_config_file(config_type_path)

    conf = Config(config_filename, start_date, end_date, quiet_mode=(not args.verbose))

    logger.info(f"Starting download for variable: {args.variable}")
    logger.info(f"Date range: {args.start} to {args.end}")
    logger.info(f"Base path: {args.base_path}")

    try:
        # Step 1: Download metadata
        # (use end date to ensure the most recent metadata is downloaded,
        # as it may contain station status updates)
        if do_download_metadata:
            download_metadata(
                conf,
                args.variable,
                args.base_path,
                end_year,
                end_month,
                end_day,
                auth_header,
                do_download=True
            )

        # Step 2: Download station list
        if do_download_stations:
            download_station_list(
                conf,
                args.variable,
                args.base_path,
                auth_header,
                do_download=True
            )

        # Step 3: Download station data
        if do_download_data:
            download_station_data(
                conf,
                args.variable,
                args.base_path,
                args.start,
                args.end,
                auth_header,
                do_download=True
            )

        # Step 4: Merge timeseries with metadata
        if do_merge:
            merge_timeseries_with_metadata(
                conf,
                args.variable,
                args.base_path,
                args.start,
                args.end,
                do_merge=True
            )

        logger.info("Download completed successfully!")

    except Exception as e:
        logger.error(f"Error during download: {e}")
        sys.exit(1)


def main_script():
    sys.exit(main())


if __name__ == "__main__":
    main_script()
