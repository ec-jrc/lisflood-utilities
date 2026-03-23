#!/usr/bin/env python
"""
This script downloads KIWI files from the WISKI API for a variable and time frame.
It downloads metadata, filters stations in the EFAS domain, downloads timeseries
data for each station, and merges the data with metadata to create KIWI files.

Usage:
    python download_timeseries.py <variable> [options]

Examples:
    python download_timeseries.py pr6
    python download_timeseries.py pr6 --start 2024-12-31 --end 2026-01-02
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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Column indices for metadata file (0-based)
MAX_METADATA_FIELDS = 18
METADATA_IDX_VALUE = 2
METADATA_IDX_QCODE = 3
METADATA_IDX_STATION_ID = 8
METADATA_IDX_NOGRIDDING = 14
METADATA_IDX_ISINARCMINDOMAIN = 16

# Column indices for timeseries files (0-based)
MAX_TIMESERIES_FIELDS = 3
TIMESERIES_IDX_TIMESTAMP = 0
TIMESERIES_IDX_VALUE = 1
TIMESERIES_IDX_QCODE = 2

# Constants
NEWLINE = '\n'
COLUMN_SEPARATOR = '\t'

DATE_FORMAT = "%Y-%m-%d"
DATE_SEPARATOR = "-"

KEY_TIMESERIES_ID = "id"
KEY_TIMESERIES_TIME_TO_DOWNLOAD = "time"
KEY_TIMESERIES_PATH_WILDCARD = "path"

# Configuration dictionary for timeseries variables
TIMESERIES_CONFIG: Dict[str, Dict[str, str]] = {
    "pr6": {KEY_TIMESERIES_ID: "1417495", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "06", KEY_TIMESERIES_PATH_WILDCARD: "*/*/Precip/6h.Total"},
    "ta6": {KEY_TIMESERIES_ID: "1417496", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "06", KEY_TIMESERIES_PATH_WILDCARD: "*/*/AT/6h.Mean"},
    "pr": {KEY_TIMESERIES_ID: "236114", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "06", KEY_TIMESERIES_PATH_WILDCARD: "*/*/Precip/D6Day.Total"},
    "tx": {KEY_TIMESERIES_ID: "236116", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "18", KEY_TIMESERIES_PATH_WILDCARD: "*/*/AT/Day.Max"},
    "tn": {KEY_TIMESERIES_ID: "236115", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "06", KEY_TIMESERIES_PATH_WILDCARD: "*/*/AT/Day.Min"},
    "pd": {KEY_TIMESERIES_ID: "236119", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "00", KEY_TIMESERIES_PATH_WILDCARD: "*/*/VP/1440min.Cmd.P"},
    "ws": {KEY_TIMESERIES_ID: "236117", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "00", KEY_TIMESERIES_PATH_WILDCARD: "*/*/WSpeed/Day.Mean"},
    "rg": {KEY_TIMESERIES_ID: "236118", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "00", KEY_TIMESERIES_PATH_WILDCARD: "*/*/SunRad/Day.Energy"},
    "wx": {KEY_TIMESERIES_ID: "236117", KEY_TIMESERIES_TIME_TO_DOWNLOAD: "00", KEY_TIMESERIES_PATH_WILDCARD: "*/*/WSpeed/Day.Max"}
}

# Environment variable for API key
ENV_KIWI_API_KEY = "KIWI_API_KEY"

# API Configuration
API_URL = "cems-mdcc.kisterscloud.de"
DEFAULT_KEY = os.getenv(ENV_KIWI_API_KEY, "")  # Allow API key to be set via environment variable


class DownloadError(Exception):
    """Custom exception for download errors."""
    pass


class APIError(Exception):
    """Custom exception for API errors."""
    pass


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


def download_metadata(variable: str, base_path: str, yyyy: str, mm: str, dd: str,
                      auth_header: str, do_download: bool = True) -> str:
    """
    Download the metadata file from the API.
    
    Args:
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
    timeseries_id = TIMESERIES_CONFIG[variable][KEY_TIMESERIES_ID]
    hh = TIMESERIES_CONFIG[variable][KEY_TIMESERIES_TIME_TO_DOWNLOAD]
    
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
        "md_returnfields": "site_no,site_name,station_name,station_no,station_id,ts_unitsymbol,ts_shortname,parametertype_name,ca_sta,ca_par",
        "crs": "global",
        "ca_sta_returnfields": "station_status,station_diary_status,EFAS-ADDATTR-NOGRIDDING,EFAS-ADDATTR-COUNTRY,EFAS-ADDATTR-ISINARCMINDOMAIN",
        "ca_par_returnfields": "EXCLUDE,INACTIVE_histattr",
        "returnfields": "sta_location,ts_value,q_code",
        "valueType": "matchingValue_nocalc"
    }
    
    url = build_url(API_URL, "EFASDATA/KiWIS", params)
    
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
        for i, line in enumerate(lines):
            fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            if i == 0:
                # Header line - keep as is
                processed_lines.append(line.rstrip(NEWLINE))
            else:
                # Data lines - replace column 3 with {value} and column 4 with {qcode}
                if len(fields) >= 4:
                    fields[METADATA_IDX_VALUE] = "{value}"
                    fields[METADATA_IDX_QCODE] = "{qcode}"
                processed_lines.append(COLUMN_SEPARATOR.join(fields))
        
        # Write back
        with open(metadata_file, 'w', encoding='utf-8') as f:
            f.write(NEWLINE.join(processed_lines) + NEWLINE)
            
        logger.info(f"Processed metadata file: {metadata_file}")
        
    except Exception as e:
        logger.error(f"Error processing metadata file: {e}")
        raise


def download_station_list(variable: str, base_path: str, auth_header: str,
                          do_download: bool = True) -> str:
    """
    Download the list of stations from the API.
    
    Args:
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
    
    ts_path_wildcard = TIMESERIES_CONFIG[variable][KEY_TIMESERIES_PATH_WILDCARD]
    
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
    
    url = build_url(API_URL, "EFASDATA/KiWIS", params)
    
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
        
        for line in lines[1:]:  # Skip header
            fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            if len(fields) >= MAX_METADATA_FIELDS:
                nogridding = fields[METADATA_IDX_NOGRIDDING] if len(fields) > METADATA_IDX_NOGRIDDING else ""
                isinarcmindomain = fields[METADATA_IDX_ISINARCMINDOMAIN] if len(fields) > METADATA_IDX_ISINARCMINDOMAIN else ""
                
                if nogridding != "yes" and isinarcmindomain != "no":
                    if len(fields) > METADATA_IDX_STATION_ID:
                        efas_station_ids.add(fields[METADATA_IDX_STATION_ID])
        
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


def download_station_data(variable: str, base_path: str, start_period: str, 
                          end_period: str, auth_header: str, 
                          do_download: bool = True) -> List[Tuple[str, str]]:
    """
    Download timeseries data for each station.
    
    Args:
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
        
        url = build_url(API_URL, "EFASDATA/KiWIS", params)
        
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


def merge_timeseries_with_metadata(variable: str, base_path: str, 
                                    start_period: str, end_period: str,
                                    do_merge: bool = True) -> None:
    """
    Merge all timeseries with the metadata to create KIWI files.
    
    Args:
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
    
    for line in lines[1:]:
        fields = line.rstrip(NEWLINE).split(COLUMN_SEPARATOR)
        if len(fields) > METADATA_IDX_STATION_ID:
            station_id = fields[METADATA_IDX_STATION_ID]
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
            
            # Skip header (first 9 lines based on the bash script)
            # The data starts after line 9 (tail -n +10 in bash)
            data_start = 9
            
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
                
                # Filter by year
                cur_year = cur_timestep[:4]
                if not cur_year.isdigit():
                    continue
                    
                if int(cur_year) < start_year or int(cur_year) > end_year:
                    continue
                
                # Parse the timestamp
                try:
                    yyyy = cur_timestep[0:4]
                    mm = cur_timestep[5:7]
                    dd = cur_timestep[8:10]
                    hh = cur_timestep[11:13]
                except (IndexError, ValueError):
                    continue
                
                # Create the output folder and file
                working_folder = os.path.join(base_folder_meteo, yyyy, mm, dd)
                os.makedirs(working_folder, exist_ok=True)
                
                kiwi_file = os.path.join(working_folder, 
                                         f"{variable}{yyyy}{mm}{dd}{hh}00.kiwi")
                
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
        choices=list(TIMESERIES_CONFIG.keys()),
        help="Variable to download (e.g., pr6, ta6, pr, tx, tn, pd, ws, rg, wx)"
    )
    parser.add_argument(
        "--base-path",
        type=str,
        default="/mnt/nahaUsers/gomesgo/CALIBRATION_6.0/download_timeseries",
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

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Parse start date for metadata download
        start_date = datetime.strptime(args.start, DATE_FORMAT)
        start_year, start_month, start_day = start_date.strftime(DATE_FORMAT).split(DATE_SEPARATOR)
    except ValueError:
        logger.error("Invalid start date format. Use YYYY-MM-DD.")
        sys.exit(1)

    # Set end date
    if args.end is None:
        # If end date is not provided, use current date
        end_date = datetime.now()
        end_year, end_month, end_day = end_date.strftime(DATE_FORMAT).split(DATE_SEPARATOR)
        args.end = f"{end_year}-{end_month}-{end_day}"
    else:
        # Validate end date format
        try:
            end_date = datetime.strptime(args.end, DATE_FORMAT)
            end_year, end_month, end_day = end_date.strftime(DATE_FORMAT).split(DATE_SEPARATOR)
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

    logger.info(f"Starting download for variable: {args.variable}")
    logger.info(f"Date range: {args.start} to {args.end}")
    logger.info(f"Base path: {args.base_path}")

    try:
        # Step 1: Download metadata (use start date to ensure data availability)
        if do_download_metadata:
            download_metadata(
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
                args.variable,
                args.base_path,
                auth_header,
                do_download=True
            )

        # Step 3: Download station data
        if do_download_data:
            download_station_data(
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
