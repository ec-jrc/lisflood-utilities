from pathlib import Path
import pandas as pd
import re
from datetime import datetime as dt

class KiwisFilter():
    
    def __init__(self, filter_columns: dict, filter_args: dict):
        self.args = filter_args
        self.filter_columns = filter_columns
        self.stati = {"Active": 1, "Inactive": 0, "yes": 0, "no": 1, "Closed": 0, "Under construction": 0}
        self.defaultReturn = 1
        self.cur_timestamp = ''
        print('ARGS: ', self.args)
        self.COL_LAT = self.__get_column_name('COL_LAT', 'station_latitude')
        self.COL_LON = self.__get_column_name('COL_LON', 'station_longitude')
        self.COL_VALUE = self.__get_column_name('COL_VALUE', 'ts_value')
        self.COL_QUALITY_CODE = self.__get_column_name('COL_QUALITY_CODE', 'q_code')
        self.COL_STATION_DIARY_STATUS = self.__get_column_name('COL_STATION_DIARY_STATUS', 'station_diary_status')
        self.COL_NO_GRIDDING = self.__get_column_name('COL_NO_GRIDDING', 'EFAS-ADDATTR-NOGRIDDING')
        self.COL_IS_IN_DOMAIN = self.__get_column_name('COL_IS_IN_DOMAIN', 'EFAS-ADDATTR-ISINARCMINDOMAIN')
        self.COL_EXCLUDE = self.__get_column_name('COL_EXCLUDE', 'EXCLUDE')
        self.COL_INACTIVE_HISTORIC = self.__get_column_name('COL_INACTIVE_HISTORIC', 'INACTIVE_histattr')
        
        self.OUTPUT_COLUMNS = [self.COL_LON, self.COL_LAT, self.COL_VALUE]

    
    def filter(self, kiwis_files: list, kiwis_timestamps: list, kiwis_data_frames: list) -> list:
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
            self.cur_timestamp = dt.strptime(f'{kiwis_timestamps[i]}00', "%Y%m%d%H%M%S")
            df_kiwis = self.__apply_filter(df_kiwis)
            filtered_data_frames.append(df_kiwis)
            i += 1
        return filtered_data_frames

    def __get_column_name(self, column_arg_key: str, column_default_name: str):
        return column_default_name if column_arg_key not in self.filter_columns else self.filter_columns[column_arg_key]
        
    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter kiwis leaving only the rows to be used for point file creation
        """
        # convert to string to make it easy to compare columns that have a mixture of string and number data
        df = df.astype(str)
        df = df.replace('nan', '')
        # Translate status columns
        df[f'{self.COL_STATION_DIARY_STATUS}_INTERNAL'] = df[self.COL_STATION_DIARY_STATUS].apply(self.rewritecol)
        df[f'{self.COL_INACTIVE_HISTORIC}_INTERNAL'] = df[self.COL_INACTIVE_HISTORIC].apply(self.rewritecol)
        # Apply filtering rules
        df = df.loc[((df[f'{self.COL_QUALITY_CODE}'] == '40') | (df[f'{self.COL_QUALITY_CODE}'] == '120')) & 
                    (df[f'{self.COL_NO_GRIDDING}'] == 'no') & (df[f'{self.COL_IS_IN_DOMAIN}'] == 'yes') &
                    (df[f'{self.COL_EXCLUDE}'] != 'yes') & (df[f'{self.COL_STATION_DIARY_STATUS}_INTERNAL'] == 1) &
                    (df[f'{self.COL_INACTIVE_HISTORIC}_INTERNAL'] == 1)]
        return df[self.OUTPUT_COLUMNS]
    
    def getStatus(self, status):
        try:
            curstate = int(self.stati[status])
            return curstate
        except:
            return 0

    def getNotStatus(self, status):
        try:
            curstate = int(not self.stati[status])
            return curstate
        except:
            return 0

    def rewritecol(self, cell_value: str):
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
                return_code = self.getNotStatus(status_list[0])
            elif to_eval_timestamp >= datetime_list[-1]:
                """Timestamp is after the last timestamp, using last timestamp as value"""
                return_code = self.getStatus(status_list[-1])
            else:
                for i in range(0, len(datetime_list)):
                    if i != len(datetime_list) - 1:
                        if datetime_list[i] <= to_eval_timestamp < datetime_list[i+1]:
                            return_code = self.getStatus(status_list[i])
        return return_code


# KIWIS_FILTER_PLUGIN_CLASSES = {'DowgradedObservationsKiwisFilter': {'1295': 1.0}, 'ObservationsKiwisFilter': {'1303': 100.0}}

class DowgradedObservationsKiwisFilter(KiwisFilter):
    
    def __init__(self, filter_columns: dict, filter_args: dict):
        super().__init__(filter_columns, filter_args)

    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().__apply_filter(df)
        return df

class ObservationsKiwisFilter(KiwisFilter):
    
    def __init__(self, filter_columns: dict, filter_args: dict):
        super().__init__(filter_columns, filter_args)

    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().__apply_filter(df)
        return df
    