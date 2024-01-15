from pathlib import Path
import numpy as np
import pandas as pd
import re
from datetime import datetime as dt
from typing import List, Tuple
from scipy.spatial import cKDTree
import geopy.distance


class KiwisFilter:
    """
    Class to filter Kiwis files metadata and obtain a dataframe containing only the coordinates and values
    to be used for interpolation. 
    """

    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}):
        self.args = filter_args
        self.filter_columns = filter_columns
        self.stati = {"Active": 1, "Inactive": 0, "yes": 0, "no": 1, "Closed": 0, "Under construction": 0}
        self.defaultReturn = 1
        self.cur_timestamp = ''
        self.setup_column_names()
        self.OUTPUT_COLUMNS = [self.COL_LON, self.COL_LAT, self.COL_VALUE]

    def filter(self, kiwis_files: List[Path], kiwis_timestamps: List[str], kiwis_data_frames: List[pd.DataFrame]) -> List[pd.DataFrame]:
        """
        Filter all kiwis files in the list and returns a list of the corresponding filtered pandas data frames.
        If the kiwis_data_frames is not empty then filter the kiwis dataframes instead of the kiwis_files 
        """
        filtered_data_frames = []
        i = 0
        for file_path in kiwis_files:
            if len(kiwis_data_frames) > 0:
                df_kiwis = kiwis_data_frames[i]
            else:
                df_kiwis = pd.read_csv(file_path, sep="\t")
            self.cur_timestamp = dt.strptime(f'{kiwis_timestamps[i]}', "%Y%m%d%H%M")
            df_kiwis = self.apply_filter(df_kiwis)
            filtered_data_frames.append(df_kiwis)
            i += 1
        return filtered_data_frames

    def apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter kiwis leaving only the rows to be used for point file creation
        """
        # convert to string to make it easy to compare columns that have a mixture of string and number data
        df = df.astype(str)
        df = df.replace('nan', '')
        # Translate status columns
        df[f'{self.COL_STATION_DIARY_STATUS}_INTERNAL'] = df[self.COL_STATION_DIARY_STATUS].apply(self.__rewrite_column)
        df[f'{self.COL_INACTIVE_HISTORIC}_INTERNAL'] = df[self.COL_INACTIVE_HISTORIC].apply(self.__rewrite_column)
        # Apply filtering rules
        df = df.loc[((df[f'{self.COL_QUALITY_CODE}'] == '40') | (df[f'{self.COL_QUALITY_CODE}'] == '120')) & 
                    (df[f'{self.COL_NO_GRIDDING}'] == 'no') & (df[f'{self.COL_IS_IN_DOMAIN}'] == 'yes') &
                    (df[f'{self.COL_EXCLUDE}'] != 'yes') & (df[f'{self.COL_STATION_DIARY_STATUS}_INTERNAL'] == 1) &
                    (df[f'{self.COL_INACTIVE_HISTORIC}_INTERNAL'] == 1)]
        return df

    def get_dataframe_output_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[self.OUTPUT_COLUMNS]

    def setup_column_names(self):
        self.COL_LAT = self.__get_column_name('COL_LAT', 'station_latitude')
        self.COL_LON = self.__get_column_name('COL_LON', 'station_longitude')
        self.COL_VALUE = self.__get_column_name('COL_VALUE', 'ts_value')
        self.COL_QUALITY_CODE = self.__get_column_name('COL_QUALITY_CODE', 'q_code')
        self.COL_PROVIDER_ID = self.__get_column_name('COL_PROVIDER_ID', 'site_no')
        self.COL_STATION_NUM = self.__get_column_name('COL_STATION_NUM', 'station_no')
        self.COL_STATION_ID = self.__get_column_name('COL_STATION_ID', 'station_id')
        self.COL_STATION_DIARY_STATUS = self.__get_column_name('COL_STATION_DIARY_STATUS', 'station_diary_status')
        self.COL_NO_GRIDDING = self.__get_column_name('COL_NO_GRIDDING', 'EFAS-ADDATTR-NOGRIDDING')
        self.COL_IS_IN_DOMAIN = self.__get_column_name('COL_IS_IN_DOMAIN', 'EFAS-ADDATTR-ISINARCMINDOMAIN')
        self.COL_EXCLUDE = self.__get_column_name('COL_EXCLUDE', 'EXCLUDE')
        self.COL_INACTIVE_HISTORIC = self.__get_column_name('COL_INACTIVE_HISTORIC', 'INACTIVE_histattr')
    
    def get_class_name(self) -> str:
        return self.__class__.__name__
    
    def get_class_description(self) -> str:
        return f'Filter: {self.get_class_name()} Args: {self.args}'

    def __get_column_name(self, column_arg_key: str, column_default_name: str) -> str:
        return column_default_name if column_arg_key not in self.filter_columns else self.filter_columns[column_arg_key]

    def __getStatus(self, status) -> int:
        try:
            curstate = int(self.stati[status])
            return curstate
        except:
            return 0

    def __getNotStatus(self, status) -> int:
        try:
            curstate = int(not self.stati[status])
            return curstate
        except:
            return 0

    def __rewrite_column(self, cell_value: str) -> int:
        return_code = self.defaultReturn
        to_eval_timestamp = self.cur_timestamp
        status_strings = None
        
        if cell_value:
            status_strings = cell_value.split("<br>")
        datetime_list = []
        status_list = []
        if status_strings:
            for crnt_string in status_strings:
                date_str, status = re.match("(\d\d\d\d-\d\d-\d\d \d\d:\d\d:\d\d): (.*)", crnt_string).groups()
                datetime_list.append(dt.strptime(date_str, "%Y-%m-%d %H:%M:%S"))
                status_list.append(status)
        if datetime_list:
            if to_eval_timestamp < datetime_list[0]:
                """Timestamp is before the first timestamp, return code should be the NOT of status"""
                return_code = self.__getNotStatus(status_list[0])
            elif to_eval_timestamp >= datetime_list[-1]:
                """Timestamp is after the last timestamp, using last timestamp as value"""
                return_code = self.__getStatus(status_list[-1])
            else:
                for i in range(0, len(datetime_list)):
                    if i != len(datetime_list) - 1:
                        if datetime_list[i] <= to_eval_timestamp < datetime_list[i+1]:
                            return_code = self.__getStatus(status_list[i])
        return return_code


class ObservationsKiwisFilter(KiwisFilter):
    """
    Class to filter Kiwis files metadata for stations that contain another station in the vicinity.
    Expects to have in filter_args a dictionary containing the provider ID whose stations we want to
    filter (as key) and the radius (in decimal degrees) to find the vicinity station from other providers (as value).
    """
    
    CLUSTER_COLLAPSE_RADIUS = 0.011582073434000193 # decimal degrees (1287 m)
    
    @staticmethod
    def kilometers2degrees(km: float) -> float:
        # Convert km to degrees of latitude
        delta_lat = km * 0.00899928005
        return delta_lat

    def apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().apply_filter(df)
        # Apply filtering rules for each provider
        for provider_id in self.args:
            radius_km = self.args[provider_id]
            radius = self.kilometers2degrees(radius_km)
            other_providers = df[(df[self.COL_PROVIDER_ID] != provider_id)]
            tree = cKDTree(other_providers[[self.COL_LON, self.COL_LAT]])
            df['has_neighbor_within_radius'] = df.apply(self.has_neighbor_within_radius_from_other_providers,
                                                        axis=1, tree=tree, provider_id=provider_id, radius=radius)
            df = df.loc[~(df['has_neighbor_within_radius'])]
        return df
    
    def has_neighbor_within_radius_from_other_providers(self, row: pd.Series, tree: cKDTree = None, provider_id: int = 0,
                                                        radius: float = CLUSTER_COLLAPSE_RADIUS) -> bool:
        cur_provider_id = row[self.COL_PROVIDER_ID]
        if cur_provider_id == provider_id:
            location = (row[self.COL_LON], row[self.COL_LAT])
            nearby_points = tree.query_ball_point(location, radius)
            return len(nearby_points) > 0
        return False


class DowgradedObservationsKiwisFilter(ObservationsKiwisFilter):
    """
    Class to filter Kiwis files metadata for stations whose daily data was down graded to 6hourly data
    by dividing the daily value by 4.
    Expects to have in filter_args a dictionary containing the provider ID whose stations we want to
    filter (as key) and the radius to find the real 6h observations from other providers (as value). 
    It will filter the station form all the 4 files if and only if at least one station in one of the
    files contains a real observation from other providers in the indicated radius.
    In the presence of a nearby real 6h station, this filter will remove all instances of the
    downgraded station from all the dataframes. 
    """

    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}):
        super().__init__(filter_columns, filter_args)
        self.filtered_station_ids = {}

    def filter(self, kiwis_files: List[Path], kiwis_timestamps: List[str], kiwis_data_frames: List[pd.DataFrame]) -> List[pd.DataFrame]:
        """
        Filter all kiwis files in the list and returns a list of the corresponding filtered pandas data frames.
        If the kiwis_data_frames is not empty then filter the kiwis dataframes instead of the kiwis_files 
        """
        # Filter all the dataframes
        filtered_data_frames = super().filter(kiwis_files, kiwis_timestamps, kiwis_data_frames)
        # Remove all stations to be filtered from all the dataframes
        return_data_frames = []
        for df in filtered_data_frames:
            df_filtered = df[~df[self.COL_STATION_ID].isin(self.filtered_station_ids)]
            return_data_frames.append(df_filtered)
        return return_data_frames

    def apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().apply_filter(df)
        # Apply filtering rules for each provider
        for provider_id in self.args:
            radius_km = self.args[provider_id]
            radius = self.kilometers2degrees(radius_km)
            other_providers = df[(df[self.COL_PROVIDER_ID] != provider_id)]
            tree = cKDTree(other_providers[[self.COL_LON, self.COL_LAT]])
            df['has_neighbor_within_radius'] = df.apply(self.has_neighbor_within_radius_from_other_providers,
                                                        axis=1, tree=tree, provider_id=provider_id, radius=radius)
            self.set_filtered_stations(df)
            df = df.loc[~(df['has_neighbor_within_radius'])]
        return df
    
    def set_filtered_stations(self, df: pd.DataFrame):
        # Get the stations to filter
        df_filtered = df.loc[(df['has_neighbor_within_radius'])]
        for station_id in df_filtered[self.COL_STATION_ID].values:
            self.filtered_station_ids[station_id] = ''
