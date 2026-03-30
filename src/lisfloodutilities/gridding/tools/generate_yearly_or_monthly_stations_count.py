import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
from dateutil import relativedelta
import sys
import argparse


TIMESTAMP_PATTERN = '%Y-%m-%d'
# For 6hourly data there are at max 4 values per day
MAX_SUB_DAILY_SAMPLES = 4


def days_of_month(current_date: str, datetime_pattern: str = TIMESTAMP_PATTERN) -> int:
    first_day_current_month = datetime.strptime(current_date, datetime_pattern).replace(day=1)
    first_day_next_month = first_day_current_month + relativedelta.relativedelta(months=1)
    last_day_current_month = first_day_next_month - timedelta(days=1)
    return int(last_day_current_month.day)

def generate_aggregated_stations(foldername: str, var_id: str, outputfolder: str, aggregate_monthly: bool = False, is_6hourly: bool = False):
    """
    Generate yearly stations count from the point files.
    In a given timestep, if there are more then one coordinate it counts as one.

    Parameters:
    - foldername (str): The name of the folder containing the point files for the year.
    - var_id (str): The current ID of the variable.
    - outputfolder (str): The output folder where the TSV results are going to be written.
    - aggregate_monthly (bool): Flag indicating if data should be aggregated monthly. By default will be aggregated yearly.
    - is_6hourly (bool): Flag indicating if data is 6hourly. By default it is daily data.
    """
    try:
        file_pattern = f'{var_id}*_*.txt'
        column_names = ['lon', 'lat', f'{var_id}']
        folderpath = Path(foldername)
        output_folder = Path(outputfolder, var_id)
        output_folder.mkdir(parents=True, exist_ok=True)
        if aggregate_monthly:
            mm = folderpath.name
            yyyy = folderpath.parent.name
            max_values_per_station = days_of_month(f'{yyyy}-{mm}-01')
            output_filename = Path(output_folder, f'stations_{var_id}_{yyyy}{mm}.tsv')
        else:
            yyyy = folderpath.name
            february_max_days = days_of_month(f'{yyyy}-02-01')
            max_values_per_station = 366 if february_max_days > 28 else 365
            output_filename = Path(output_folder, f'stations_{var_id}_{yyyy}.tsv')

        if is_6hourly:
            max_values_per_station *= MAX_SUB_DAILY_SAMPLES

        files = list(folderpath.rglob(file_pattern))
        if len(files) == 0:
            raise Exception('No input files were found.')
        
        # DEBUG: Log file count for diagnostics
        print(f"[DEBUG] Found {len(files)} files matching pattern '{file_pattern}'")
        
        agg_data = None
        for file in files:
            df = pd.read_csv(file, names=column_names, header=None, sep='\t')
            grouped = df.groupby(['lon', 'lat']).size().reset_index(name='count')
            
            # DEBUG: Log grouped DataFrame info
            print(f"[DEBUG] File: {file.name}, rows: {len(df)}, grouped rows: {len(grouped)}, columns: {list(grouped.columns)}")
            
            if agg_data is None:
                agg_data = grouped
            else:
                agg_data = pd.concat([agg_data, grouped]).groupby(['lon', 'lat'], as_index=False)['count'].sum()
        
        # DEBUG: Validate agg_data before operations
        print(f"[DEBUG] agg_data is None: {agg_data is None}")
        if agg_data is not None:
            print(f"[DEBUG] agg_data shape: {agg_data.shape}, columns: {list(agg_data.columns)}")
        
        # Ensure agg_data is not None and has 'count' column before operations
        if agg_data is None or agg_data.empty:
            raise ValueError("No data was aggregated from input files.")
        
        if 'count' not in agg_data.columns:
            raise ValueError("The 'count' column is missing from aggregated data. This may occur if input files are empty.")
        
        agg_data['count'] = agg_data['count'].clip(upper=max_values_per_station)
        agg_data.to_csv(output_filename, sep='\t', index=False)
        print(f"Wrote {output_filename}")
    except FileNotFoundError:
        print(f"Error: File '{foldername}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate yearly or monthly stations count.")
    # set defaults
    parser.set_defaults(aggregate_monthly=False, is_6hourly=False)
    # set arguments
    parser.add_argument("foldername", help="The name of the yearly or monthly folder.")
    parser.add_argument("-v", "--var-id", required=True, help="The current ID of the variable.")
    parser.add_argument("-6", "--6hourly", dest="is_6hourly", action="store_true",
    			help="Informs the data to be processed is 6hourly, otherwise it is daily. [default: %(default)s]")
    parser.add_argument("-o", "--outputfolder", required=True, help="The output folder for the resulting TSV files.")
    parser.add_argument("-a", "--monagg", dest="aggregate_monthly", action="store_true",
    			help="Aggregate monthly, otherwise aggregates yearly. [default: %(default)s]")
    args = parser.parse_args()

    generate_aggregated_stations(args.foldername, args.var_id, args.outputfolder, args.aggregate_monthly, args.is_6hourly)

if __name__ == "__main__":
    main()
