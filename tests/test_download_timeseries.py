"""
Unit tests for download_timeseries.py

These tests cover the main functions in the download_timeseries module.
Since the module interacts with external APIs, most tests use mocking.

Usage:
    pytest tests/test_download_timeseries.py -v

Copyright 2019-2020 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission 
subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
"""

import os
import sys
import tempfile
import shutil
import pytest
from pathlib import Path
from datetime import datetime


@pytest.fixture(autouse=True)
def reset_metadata_indices():
    """
    Reset the global metadata_indices variable before each test.
    
    This ensures that tests don't interfere with each other due to
    cached column indices from previous tests.
    """
    import lisfloodutilities.gridding.download_timeseries as download_timeseries
    
    # Reset the global metadata_indices dictionary
    download_timeseries.metadata_indices = {}
    
    yield
    
    # Also reset after the test to ensure clean state
    download_timeseries.metadata_indices = {}


# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lisfloodutilities.gridding.download_timeseries import (
    calculate_metadata_indices,
    create_auth_header,
    build_url,
    check_file_for_errors,
    process_metadata_file,
    filter_stations_efas_domain,
    merge_timeseries_with_metadata,
    METADATA_COL_VALUE,
    METADATA_COL_QCODE,
    METADATA_COL_STATION_ID,
    METADATA_COL_NOGRIDDING,
    METADATA_COL_ISINARCMINDOMAIN,
    METADATA_MAX_FIELDS_KEY,
    COLUMN_SEPARATOR,
    NEWLINE,
)
from lisfloodutilities.gridding.lib.utils import Config


# Test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'download_timeseries', 'reference')


class TestCheckFileForErrors:
    """Tests for check_file_for_errors function."""

    def test_check_file_with_html_error(self):
        """Test detection of HTML error content."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False) as f:
            f.write("<html><body>Error 404</body></html>")
            temp_file = f.name
        
        try:
            result = check_file_for_errors(temp_file)
            assert result is True
        finally:
            os.unlink(temp_file)

    def test_check_file_with_valid_content(self):
        """Test that valid content is not flagged as error."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("station_id\tvalue\n123\t10.5")
            temp_file = f.name
        
        try:
            result = check_file_for_errors(temp_file)
            assert result is False
        finally:
            os.unlink(temp_file)

    def test_check_nonexistent_file(self):
        """Test handling of nonexistent file."""
        result = check_file_for_errors("/nonexistent/file.txt")
        assert result is True  # Should return True (error) for missing file


class TestProcessMetadataFile:
    """Tests for process_metadata_file function."""

    def test_process_metadata_file(self):
        """Test that metadata file is processed correctly."""
        # Create a temporary metadata file
        header = "station_latitude\tstation_longitude\tts_value\tq_code\tsite_no\tstation_id\tEFAS-ADDATTR-NOGRIDDING\tEFAS-ADDATTR-ISINARCMINDOMAIN"
        data_line = "44.702144\t10.956703\t15.73\t40\t1221\t173904\tno\tyes"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_file = os.path.join(tmpdir, "metadata.tsv")
            
            with open(metadata_file, 'w') as f:
                f.write(header + NEWLINE)
                f.write(data_line + NEWLINE)
            
            # Process the file
            process_metadata_file(metadata_file)
            
            # Read and verify
            with open(metadata_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2
            assert lines[0].strip() == header
            # Value and qcode should be replaced with placeholders
            assert "{value}" in lines[1]
            assert "{qcode}" in lines[1]
            assert "15.73" not in lines[1]
            assert "40" not in lines[1]


class TestFilterStationsEfasDomain:
    """Tests for filter_stations_efas_domain function."""

    def test_filter_stations_efas_domain(self):
        """Test filtering stations in EFAS domain."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create metadata file
            metadata_file = os.path.join(tmpdir, "metadata.tsv")
            metadata_content = (
                "station_latitude\tstation_longitude\tts_value\tq_code\tstation_id\t"
                "EFAS-ADDATTR-NOGRIDDING\tEFAS-ADDATTR-ISINARCMINDOMAIN\n"
                "44.702144\t10.956703\t15.73\t40\t173904\tno\tyes\n"  # In EFAS domain
                "44.089031\t12.27459\t16.57\t40\t173905\tno\tyes\n"  # In EFAS domain
                "50.0\t10.0\t10.0\t40\t173906\tyes\tno\n"  # Not in EFAS domain (nogridding)
            )
            with open(metadata_file, 'w') as f:
                f.write(metadata_content)
            
            # Create stations file
            stations_file = os.path.join(tmpdir, "stations.tsv")
            stations_content = (
                "station_id\tts_path\n"
                "173904\tpath/to/173904\n"
                "173905\tpath/to/173905\n"
                "173906\tpath/to/173906\n"
            )
            with open(stations_file, 'w') as f:
                f.write(stations_content)
            
            # Output files
            efas_domain_file = os.path.join(tmpdir, "efas_domain.tsv")
            output_file = os.path.join(tmpdir, "filtered_stations.tsv")
            
            # Run filter
            filter_stations_efas_domain(metadata_file, stations_file, efas_domain_file, output_file)

            # Verify output
            with open(output_file, 'r') as f:
                lines = f.readlines()

            # Should have header + 2 stations (173904 and 173905)
            assert len(lines) == 3
            assert "173904" in lines[1].rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            assert "173905" in lines[2].rstrip(NEWLINE).split(COLUMN_SEPARATOR)
            # 173906 should be filtered out
            assert "173906" not in "".join(lines)


class TestMergeTimeseriesWithMetadata:
    """Tests for merge_timeseries_with_metadata function."""

    @pytest.fixture
    def temp_workspace(self):
        """Create a temporary workspace with test data."""
        tmpdir = tempfile.mkdtemp()
        
        # Copy reference data to temp directory
        reference_dir = TEST_DATA_DIR
        
        # Create metadata file
        metadata_src = os.path.join(reference_dir, "tx_timeseries_metadata.tsv")
        metadata_dst = os.path.join(tmpdir, "tx_timeseries_metadata.tsv")
        shutil.copy(metadata_src, metadata_dst)
        
        # Create timeseries directory and copy files
        timeseries_dir = os.path.join(tmpdir, "timeseries", "tx")
        os.makedirs(timeseries_dir, exist_ok=True)
        
        ts_src1 = os.path.join(reference_dir, "timeseries", "tx", "tx_station_173904_timeseries.tsv")
        ts_src2 = os.path.join(reference_dir, "timeseries", "tx", "tx_station_173905_timeseries.tsv")
        shutil.copy(ts_src1, os.path.join(timeseries_dir, "tx_station_173904_timeseries.tsv"))
        shutil.copy(ts_src2, os.path.join(timeseries_dir, "tx_station_173905_timeseries.tsv"))
        
        yield tmpdir
        
        # Cleanup
        shutil.rmtree(tmpdir)

    def test_merge_timeseries_with_metadata(self, temp_workspace):
        """Test merging timeseries with metadata creates correct KIWI files."""
        # Create a mock config object
        class MockConfig:
            input_timestamp_pattern = "tx%Y%m%d%H%M.kiwi"
        
        conf = MockConfig()
        variable = "tx"
        base_path = temp_workspace
        start_period = "2026-03-20"
        end_period = "2026-03-24"
        
        # Run merge
        merge_timeseries_with_metadata(conf, variable, base_path, start_period, end_period, do_merge=True)
        
        # Verify output files exist
        expected_files = [
            os.path.join(base_path, "meteo", "tx", "2026", "03", "20", "tx202603201800.kiwi"),
            os.path.join(base_path, "meteo", "tx", "2026", "03", "21", "tx202603211800.kiwi"),
            os.path.join(base_path, "meteo", "tx", "2026", "03", "22", "tx202603221800.kiwi"),
            os.path.join(base_path, "meteo", "tx", "2026", "03", "23", "tx202603231800.kiwi"),
            os.path.join(base_path, "meteo", "tx", "2026", "03", "24", "tx202603241800.kiwi"),
        ]
        
        for expected_file in expected_files:
            assert os.path.exists(expected_file), f"Expected file not found: {expected_file}"
        
        # Verify content matches reference
        reference_dir = TEST_DATA_DIR
        for expected_file in expected_files:
            # Get relative path from meteo/tx
            rel_path = os.path.relpath(expected_file, os.path.join(base_path, "meteo", "tx"))
            reference_file = os.path.join(reference_dir, "meteo", "tx", rel_path)
            
            if os.path.exists(reference_file):
                with open(expected_file, 'r') as f:
                    actual_content = f.read()
                with open(reference_file, 'r') as f:
                    expected_content = f.read()
                
                assert actual_content == expected_content, f"Content mismatch for {expected_file}"

    def test_merge_no_data_when_disabled(self, temp_workspace):
        """Test that merge is skipped when do_merge is False."""
        class MockConfig:
            input_timestamp_pattern = "tx%Y%m%d%H%M.kiwi"
        
        conf = MockConfig()
        variable = "tx"
        base_path = temp_workspace
        start_period = "2026-03-20"
        end_period = "2026-03-24"
        
        # Run merge with do_merge=False
        merge_timeseries_with_metadata(conf, variable, base_path, start_period, end_period, do_merge=False)
        
        # Verify no meteo directory was created
        meteo_dir = os.path.join(base_path, "meteo")
        assert not os.path.exists(meteo_dir)


class TestIntegrationWithReferenceData:
    """Integration tests using reference data."""

    def test_full_workflow_simulation(self):
        """Simulate the full workflow using reference data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            reference_dir = TEST_DATA_DIR
            
            # Step 1: Process metadata file
            metadata_src = os.path.join(reference_dir, "tx_timeseries_metadata.tsv")
            metadata_dst = os.path.join(tmpdir, "tx_timeseries_metadata.tsv")
            shutil.copy(metadata_src, metadata_dst)
            
            # Verify metadata processing
            process_metadata_file(metadata_dst)
            
            with open(metadata_dst, 'r') as f:
                content = f.read()
            assert "{value}" in content
            assert "{qcode}" in content
            
            # Step 2: Filter stations
            stations_raw = os.path.join(reference_dir, "tx_stations_to_download_raw.tsv")
            stations_filtered = os.path.join(reference_dir, "tx_stations_to_download.tsv")
            efas_domain = os.path.join(reference_dir, "tx_station_ids_efas_domain.tsv")
            
            stations_raw_dst = os.path.join(tmpdir, "tx_stations_to_download_raw.tsv")
            stations_filtered_dst = os.path.join(tmpdir, "tx_stations_to_download.tsv")
            efas_domain_dst = os.path.join(tmpdir, "tx_station_ids_efas_domain.tsv")
            
            shutil.copy(stations_raw, stations_raw_dst)
            shutil.copy(efas_domain, efas_domain_dst)
            
            filter_stations_efas_domain(
                metadata_dst,
                stations_raw_dst,
                efas_domain_dst,
                stations_filtered_dst
            )
            
            # Verify filtered stations match reference
            with open(stations_filtered_dst, 'r') as f:
                actual = f.read()
            with open(stations_filtered, 'r') as f:
                expected = f.read()
            assert actual == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
