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

    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}, var_code: str = '', quiet_mode: bool = False):
        self.args = filter_args
        self.filter_columns = filter_columns
        self.var_code = var_code
        self.quiet_mode = quiet_mode
        self.QUALITY_CODE_VALID = '40'
        self.QUALITY_CODE_SUSPICIOUS = '120'
        self.QUALITY_CODE_WRONG = '160'
        self.stati = {"Active": 1, "Inactive": 0, "yes": 0, "no": 1, "Closed": 0, "Under construction": 0}
        self.defaultReturn = 1
        self.cur_timestamp = ''
        self.setup_column_names()
        self.OUTPUT_COLUMNS = [self.COL_LON, self.COL_LAT, self.COL_VALUE]
        self.INTERNAL_COLUMNS = [f'{self.COL_STATION_DIARY_STATUS}_INTERNAL', f'{self.COL_INACTIVE_HISTORIC}_INTERNAL']

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
        # Apply filtering rules but get all quality codes for the statistic
        df = df.loc[((df[f'{self.COL_QUALITY_CODE}'] == self.QUALITY_CODE_VALID) |
                     (df[f'{self.COL_QUALITY_CODE}'] == self.QUALITY_CODE_SUSPICIOUS) |
                     (df[f'{self.COL_QUALITY_CODE}'] == self.QUALITY_CODE_WRONG)) & 
                    (df[f'{self.COL_NO_GRIDDING}'] == 'no') & (df[f'{self.COL_IS_IN_DOMAIN}'] == 'yes') &
                    (df[f'{self.COL_EXCLUDE}'] != 'yes') & (df[f'{self.COL_STATION_DIARY_STATUS}_INTERNAL'] == 1) &
                    (df[f'{self.COL_INACTIVE_HISTORIC}_INTERNAL'] == 1)]
        self.print_statistics(df)
        # Show only codes valid and suspicious
        df = df.loc[((df[f'{self.COL_QUALITY_CODE}'] == self.QUALITY_CODE_VALID) |
                     (df[f'{self.COL_QUALITY_CODE}'] == self.QUALITY_CODE_SUSPICIOUS))]
        return df

    @staticmethod
    def get_totals_by_quality_code(row: pd.Series, column_quality_code: str, quality_code: str) -> int:
        cur_quality_code = row[column_quality_code]
        if cur_quality_code == quality_code:
            return row['count']
        return 0

    def print_statistics(self, df: pd.DataFrame):
        timestamp = self.cur_timestamp.strftime('%Y-%m-%d %H:%M:%S')
        new_df = df.groupby([self.COL_PROVIDER_ID, self.COL_QUALITY_CODE]).size().reset_index(name='count')
        # Transpose the quality codes
        new_df[self.QUALITY_CODE_VALID] = new_df.apply(KiwisFilter.get_totals_by_quality_code, axis=1,
                                                       column_quality_code=self.COL_QUALITY_CODE,
                                                       quality_code=self.QUALITY_CODE_VALID)
        new_df[self.QUALITY_CODE_SUSPICIOUS] = new_df.apply(KiwisFilter.get_totals_by_quality_code, axis=1,
                                                            column_quality_code=self.COL_QUALITY_CODE,
                                                            quality_code=self.QUALITY_CODE_SUSPICIOUS)
        new_df[self.QUALITY_CODE_WRONG] = new_df.apply(KiwisFilter.get_totals_by_quality_code, axis=1,
                                                       column_quality_code=self.COL_QUALITY_CODE,
                                                       quality_code=self.QUALITY_CODE_WRONG)
        new_df.drop(columns=[self.COL_QUALITY_CODE, 'count'], inplace=True)
        new_df = new_df.groupby(self.COL_PROVIDER_ID)[[self.QUALITY_CODE_VALID,
                                                       self.QUALITY_CODE_SUSPICIOUS,
                                                       self.QUALITY_CODE_WRONG]].sum()
        new_df.reset_index(inplace=True)
        for index, row in new_df.iterrows():
            provider_id = row[self.COL_PROVIDER_ID]
            quality_code_valid = row[self.QUALITY_CODE_VALID]
            quality_code_suspicious = row[self.QUALITY_CODE_SUSPICIOUS]
            quality_code_wrong = row[self.QUALITY_CODE_WRONG]
            total = quality_code_valid + quality_code_suspicious + quality_code_wrong
            stats_string = (
                f'#KIWIS_STATS: {{"TIMESTAMP": "{timestamp}", "VAR_CODE": "{self.var_code}", '
                f'"PROVIDER_ID": {provider_id}, "QUALITY_CODE_VALID": {quality_code_valid}, '
                f'"QUALITY_CODE_SUSPICIOUS": {quality_code_suspicious}, "QUALITY_CODE_WRONG": {quality_code_wrong}, '
                f'"TOTAL_OBSERVATIONS": {total}}}'
            )
            self.print_msg(stats_string)

    def print_msg(self, msg: str = ''):
        if not self.quiet_mode:
            print(msg)

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
    
    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}, var_code: str = '', quiet_mode: bool = False):
        super().__init__(filter_columns, filter_args, var_code, quiet_mode)
        self.INTERNAL_COLUMNS.append('has_neighbor_within_radius')
        # Calculating the radius in decimal degrees
        self.provider_radius = {}
        for provider_id in self.args:
            radius_km = self.args[provider_id]
            radius = self.kilometers2degrees(radius_km)
            self.provider_radius[provider_id] = radius
    
    @staticmethod
    def kilometers2degrees(km: float) -> float:
        # Convert km to degrees of latitude
        delta_lat = km * 0.00899928005
        return delta_lat

    def apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().apply_filter(df)
        # Apply filtering rules for each provider
        for provider_id in self.provider_radius:
            radius = self.provider_radius[provider_id]
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

    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}, var_code: str = '', quiet_mode: bool = False):
        super().__init__(filter_columns, filter_args, var_code, quiet_mode)
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
        for provider_id in self.provider_radius:
            radius = self.provider_radius[provider_id]
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


class DowgradedDailyTo6HourlyObservationsKiwisFilter(ObservationsKiwisFilter):
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

    def __init__(self, filter_columns: dict = {}, filter_args: dict = {}, var_code: str = '', quiet_mode: bool = False):
        super().__init__(filter_columns, filter_args, var_code, quiet_mode)
        self.INTERNAL_COLUMNS.append('has_neighbor_within_radius_complete_6h')
        self.kiwis_24h_dataframe = None
        self.df_24h_without_neighbors = None
        self.df_24h_with_neighbors = None
        self.all_6h_df = None
        self.num_6h_slots = 0
        self.PROVIDER_DWD_SYNOP = '1121'
        self.PROVIDER_MARS = '1295'
    
    def get_all_6hourly_station_values_df(self, kiwis_data_frames: List[pd.DataFrame]) -> pd.DataFrame:
        merged_df = pd.concat(kiwis_data_frames)
        merged_df = merged_df[[self.COL_LON, self.COL_LAT, self.COL_PROVIDER_ID, self.COL_STATION_NUM, self.COL_STATION_ID, self.COL_VALUE]]
        merged_df.reset_index(drop=True, inplace=True)
        result_df = merged_df.astype({self.COL_VALUE: 'float'}).groupby([self.COL_LON, self.COL_LAT,
                                                                         self.COL_PROVIDER_ID,
                                                                         self.COL_STATION_NUM,
                                                                         self.COL_STATION_ID])[self.COL_VALUE].agg(['sum','count']).reset_index()
        result_df.columns = [self.COL_LON, self.COL_LAT, self.COL_PROVIDER_ID, self.COL_STATION_NUM, self.COL_STATION_ID, 'sum_6h_values', 'count_6h_slots']
        result_df.reset_index(drop=True, inplace=True)
        return result_df

    def format_dwd_synop_wmo_num(self, station_num: str) -> str:
        """
        Used to get the DWD Synop sations with WMO Number to be formatted like the ones from MARS to allow comparison
        """
        try:
            station_num_float = float(station_num)
            return f'{station_num_float:.1f}'
        except:
            return station_num
    
    def get_decumulated_24h_value_for_missing_6h_values(self, row: pd.Series, tree: cKDTree = None, provider_id: int = 0,
                                                        radius: float = ObservationsKiwisFilter.CLUSTER_COLLAPSE_RADIUS,
                                                        stations_6h_df: pd.DataFrame = None) -> float:
        """
        DECUMULATED_PR = (PR - Sum(PR6)) / (number of missing values)
        If there are more than one 6h station in the radius, select according to the following rules by order:
        1. get the one with the highest number of slots
        2. get the one with the lowest positive difference to the 24h value
        3. get the one on the top
        """
        cur_provider_id = row[self.COL_PROVIDER_ID]
        cur_24h_value = row[self.COL_VALUE]
        decumulated_value = cur_24h_value
        if cur_provider_id == provider_id:
            location = (row[self.COL_LON], row[self.COL_LAT])
            nearby_points = tree.query_ball_point(location, radius)
            stations_6h = stations_6h_df.loc[stations_6h_df.index.isin(nearby_points)][['sum_6h_values', 'count_6h_slots']]
            if stations_6h.empty:
                return decumulated_value
            sum_6h_values, count_6h_slots = stations_6h.iloc[0][['sum_6h_values', 'count_6h_slots']]

            if len(nearby_points) > 1:
                stations_6h = stations_6h.assign(diff_to_24h=(cur_24h_value - stations_6h['sum_6h_values']))
                stations_6h = stations_6h.sort_values(by=['count_6h_slots', 'diff_to_24h'], ascending=[False, True]).reset_index(drop=True)
                stations_6h = stations_6h[stations_6h['diff_to_24h'] >= 0]
                if not stations_6h.empty:
                    stations_6h.reset_index(drop=True, inplace=True)
                    sum_6h_values, count_6h_slots = stations_6h.iloc[0][['sum_6h_values', 'count_6h_slots']]

            missing_6h_slots = self.num_6h_slots - count_6h_slots
            if missing_6h_slots > 0:
                decumulated_value = (cur_24h_value - sum_6h_values) / missing_6h_slots
            else:
                decumulated_value = cur_24h_value / self.num_6h_slots
        return decumulated_value
    
    def update_column_if_provider_stations_are_in_radius(self, df: pd.DataFrame, tree: cKDTree, column_to_update: str = 'has_neighbor_within_radius') -> pd.DataFrame:
        for provider_id in self.provider_radius:
            radius = self.provider_radius[provider_id]
            df.loc[df[self.COL_PROVIDER_ID] == provider_id, column_to_update] = df.apply(self.has_neighbor_within_radius_from_other_providers,
                                                                                         axis=1, tree=tree, provider_id=provider_id, radius=radius)
        return df

    def filter(self, kiwis_files: List[Path], kiwis_timestamps: List[str], kiwis_data_frames: List[pd.DataFrame]) -> List[pd.DataFrame]:
        """
        Filter all kiwis files in the list and returns a list of the corresponding filtered pandas data frames.
        If the kiwis_data_frames is not empty then filter the kiwis dataframes instead of the kiwis_files 
        """
        # Guarantee datatype of value column
        for i in range(len(kiwis_data_frames)):
            kiwis_data_frames[i] = kiwis_data_frames[i].astype({
                                                                # self.COL_LON: 'float',
                                                                # self.COL_LAT: 'float',
                                                                # self.COL_PROVIDER_ID: 'int',
                                                                # self.COL_STATION_ID: 'int',
                                                                self.COL_VALUE: 'float'})

        self.kiwis_24h_dataframe = kiwis_data_frames[0]
        kiwis_6h_dataframes = kiwis_data_frames[1:]
        self.num_6h_slots = len(kiwis_6h_dataframes)
        providers_ids = list(self.provider_radius.keys())

        # Filter all the dataframes
        filtered_data_frames = super().filter(kiwis_files[1:], kiwis_timestamps[1:], kiwis_6h_dataframes)

        # Get all 6hourly slots aggregated
        self.all_6h_df = self.get_all_6hourly_station_values_df(filtered_data_frames)

        # Select the daily stations from the providers that do not contain a 6hourly station within radius and decumulate
        # its value dividing by the number of 6hourly slots. These will be inserted directly in all the 6hourly slots
        tree_all_6h = cKDTree(self.all_6h_df[[self.COL_LON, self.COL_LAT]])
        
        self.kiwis_24h_dataframe['has_neighbor_within_radius'] = False
        self.kiwis_24h_dataframe = self.update_column_if_provider_stations_are_in_radius(df=self.kiwis_24h_dataframe, tree=tree_all_6h)
        df_24h_from_providers = self.kiwis_24h_dataframe[self.kiwis_24h_dataframe[self.COL_PROVIDER_ID].isin(providers_ids)]
        self.df_24h_without_neighbors = df_24h_from_providers[df_24h_from_providers['has_neighbor_within_radius'] == False]

        # From the station without neighbors remove MARS stations that have a DWDSynop with the same WMO Number because
        # they are actually the same station even though they might not be exactly within the defined radius.
        wmo_nums_of_6h_dwd_stations = self.all_6h_df.loc[self.all_6h_df[self.COL_PROVIDER_ID] == self.PROVIDER_DWD_SYNOP, self.COL_STATION_NUM].values
        # correct wmo num from DWD Synop into MARS format
        wmo_nums_of_6h_dwd_stations = map(self.format_dwd_synop_wmo_num, wmo_nums_of_6h_dwd_stations)
        # Remove MARS stations that do not have neighbors but have a 6h DWDSynop station with the same wmo num
        mask = ((self.df_24h_without_neighbors[self.COL_PROVIDER_ID] == self.PROVIDER_MARS) &
                (self.df_24h_without_neighbors[self.COL_STATION_NUM].isin(wmo_nums_of_6h_dwd_stations)))
        self.df_24h_without_neighbors = self.df_24h_without_neighbors[~mask]

        # Decumulate the values of stations without neighbors
        self.df_24h_without_neighbors[self.COL_VALUE] /= self.num_6h_slots

        # Remove the daily rows that have a complete 6hourly station in the radius because there is no need to decumulate in these cases
        complete_6h_df = self.all_6h_df.loc[self.all_6h_df['count_6h_slots'] == self.num_6h_slots]
        self.df_24h_with_neighbors = df_24h_from_providers[df_24h_from_providers['has_neighbor_within_radius'] == True]
        tree_complete_6h = cKDTree(complete_6h_df[[self.COL_LON, self.COL_LAT]])
        self.df_24h_with_neighbors['has_neighbor_within_radius_complete_6h'] = False
        self.df_24h_with_neighbors = self.update_column_if_provider_stations_are_in_radius(df=self.df_24h_with_neighbors, tree=tree_complete_6h,
                                                                                           column_to_update='has_neighbor_within_radius_complete_6h')
        self.df_24h_with_neighbors = self.df_24h_with_neighbors[self.df_24h_with_neighbors['has_neighbor_within_radius_complete_6h'] == False]

        # Change the 24h stations corresponding with the 6h stations containing missing values by changing its value using the formula
        # (PR - Sum(PR6)) / (number of missing values), if and only if the resulting value is positive (>=0).
        missing_6h_slots_df = self.all_6h_df.loc[self.all_6h_df['count_6h_slots'] < self.num_6h_slots]
        missing_6h_slots_df.reset_index(drop=True, inplace=True)
        tree_missing_6h = cKDTree(missing_6h_slots_df[[self.COL_LON, self.COL_LAT]])
        for provider_id in self.provider_radius:
            radius = self.provider_radius[provider_id]
            self.df_24h_with_neighbors.loc[self.df_24h_with_neighbors[self.COL_PROVIDER_ID] == provider_id,
                                           self.COL_VALUE] = self.df_24h_with_neighbors.apply(self.get_decumulated_24h_value_for_missing_6h_values,
                                                                                              axis=1, tree=tree_missing_6h, provider_id=provider_id,
                                                                                              radius=radius, stations_6h_df=missing_6h_slots_df)
        self.df_24h_with_neighbors = self.df_24h_with_neighbors.loc[self.df_24h_with_neighbors[self.COL_VALUE] >= 0.0]

        # Clean the dataframes of internal columns before merging them
        self.df_24h_without_neighbors = self.df_24h_without_neighbors.drop(columns=self.INTERNAL_COLUMNS, errors='ignore')

        # Insert the decumulated stations in the respective 6h slots
        return_data_frames = [self.kiwis_24h_dataframe]
        for df in filtered_data_frames:
            df = df.drop(columns=self.INTERNAL_COLUMNS, errors='ignore')
            # Now we need to eliminate the stations that have neighbors on this 6h slot,
            # which means the slot of 6h have already a 6h value in the radius and no
            # decumulation is needed in this slot.
            df_decumulated_24h = self.df_24h_with_neighbors
            df_decumulated_24h['has_neighbor_within_radius'] = False
            tree = cKDTree(df[[self.COL_LON, self.COL_LAT]])
            df_decumulated_24h = self.update_column_if_provider_stations_are_in_radius(df=df_decumulated_24h, tree=tree)
            df_decumulated_24h = df_decumulated_24h[df_decumulated_24h['has_neighbor_within_radius'] == False]
            df_decumulated_24h.drop(columns=self.INTERNAL_COLUMNS, errors='ignore')
            # Append both decumulated 24h dataframes to the 6h slot
            df_filtered = pd.concat([df, self.df_24h_without_neighbors, df_decumulated_24h])
            return_data_frames.append(df_filtered)
        return return_data_frames

