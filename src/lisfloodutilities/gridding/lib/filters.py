from pathlib import Path
import pandas as pd

class KiwisFilter():
    
    def __init__(self, filter_args: dict):
        self.args = filter_args
    
    def filter(self, kiwis_files: array) -> array:
        filtered_data_frames = []
        for file_path in kiwis_files:
            df_kiwis = pd.read_csv(file_path, sep="\t")
            df_kiwis = self.__apply_filter(df_kiwis)
            filtered_data_frames.append(df_kiwis)
        return filtered_data_frames

    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        # Get the code to filter kiwis leaving only the rows to be used for point file creation
        return df

# KIWIS_FILTER_PLUGIN_CLASSES = {'DowgradedObservationsKiwisFilter': {'1295': 1.0}, 'ObservationsKiwisFilter': {'1303': 100.0}}

class DowgradedObservationsKiwisFilter(KiwisFilter):
    
    def __init__(self, filter_args: dict):
        super().__init__(filter_args)

    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().__apply_filter(df)
        return df

class ObservationsKiwisFilter(KiwisFilter):
    
    def __init__(self, filter_args: dict):
        super().__init__(filter_args)

    def __apply_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        df = super().__apply_filter(df)
        return df
    