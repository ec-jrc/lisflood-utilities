import json
import os
import sys
import re
import math
import itertools
import time
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, least_squares
import multiprocessing as mp
from pyg2p.main.interpolation.latlong import LatLong
from pyg2p.main.readers.netcdf import NetCDFReader
from pyg2p.main.interpolation.scipy_interpolation_lib import ScipyInterpolation
from lisfloodutilities.writers import NetCDFWriter


DEBUG_CDDMAP = False
# DEBUG_MIN_LAT = 57.95-20
# DEBUG_MIN_LON = 16.80-0
# DEBUG_MAX_LAT = 57.95+0
# DEBUG_MAX_LON = 16.80+20
DEBUG_MIN_LAT = 45-5
DEBUG_MIN_LON = 8-5
DEBUG_MAX_LAT = 45+5
DEBUG_MAX_LON = 8+5



EARTH_RADIUS = 6371  # in km
MIN_COMMON_STATION_VALUES = 365*4  # no of common station values (one year at 6hourly)
CDD_INITIAL = 100  # CDD initial value for curve_fit

def read_files(base_dir):
    # Initialize empty dictionaries to hold the timeseries data and any conflicts
    timeseries = {}
    conflicted_timeseries = {}

    # Get the file paths for all txt files in the directory
    files = []
    for root, dirs, file_names in os.walk(base_dir):
        files += [os.path.join(root, file_name) for file_name in file_names if file_name.endswith(".txt")]

    # Loop through all the files and read the data
    num_files = len(files)
    progress_step = max(1,int(num_files/250))
    for i, file_path in enumerate(files):
        # Extract the timestamp from the file name
        file_name = os.path.basename(file_path)
        timestamp = re.search(r"pr6(\d{12})", file_name).group(1)
        year = timestamp[:4]
        month = timestamp[4:6]
        day = timestamp[6:8]
        hour = timestamp[8:10]
        minute = timestamp[10:12]

        # Open the txt file and read the data
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                # Split the line by tab and extract the x, y, and z values
                x, y, z = line.strip().split("\t")
                # Create the timeseries key using the combination of x and y
                key = f"{x}_{y}"
                # Construct the timestamp for the data using the timestamp from the file name
                timestamp_str = f"{year}-{month}-{day} {hour}:{minute}"
                # Check if the key is already present in the timeseries dict
                if key not in timeseries:
                    # Create a new empty dictionary for the key in the timeseries dict and add the data
                    timeseries[key] = {}
                    timeseries[key][timestamp_str] = float(z)
                else:
                    # Check if the timestamp already exists in the timeseries dict
                    if timestamp_str in timeseries[key]:
                        # Add the conflicting data to the conflicted_timeseries dict
                        if key not in conflicted_timeseries:
                            conflicted_timeseries[key] = {}
                        num_conflicts = len(conflicted_timeseries[key].get(timestamp_str, []))
                        conflicted_timeseries[key][timestamp_str] = {"new": (float(z), num_conflicts + 1), "old": timeseries[key][timestamp_str]}
                    else:
                        # Add the z value to the timeseries dict for the corresponding timestamp
                        timeseries[key][timestamp_str] = float(z)

        # Print progress message
        if i % progress_step == 0:
            print('\rReading Files: {:.2f}%'.format(i * 100. / num_files), end="")
    print("\n")

    return timeseries, conflicted_timeseries

def compute_distances(keys_subset, keys):
    # Calculate distance matrix (fast algorithm...)
    lon_lats_subset = np.array([key.split("_") for key in keys_subset], dtype=np.float64)
    lon_lats_subset = np.radians(lon_lats_subset)
    lon_lats = np.array([key.split("_") for key in keys], dtype=np.float64)
    lon_lats = np.radians(lon_lats)
    lon_diffs = lon_lats_subset[:, 0, None] - lon_lats[:, 0]
    lat_diffs = lon_lats_subset[:, 1, None] - lon_lats[:, 1]
    a = np.sin(lat_diffs / 2) ** 2 + np.cos(lon_lats_subset[:, 1, None]) * np.cos(lon_lats[:, 1]) * np.sin(lon_diffs / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distances = EARTH_RADIUS * c

    return distances

def compute_correlations(keys_subset):
    correlations = {}
    distances = {}
    i_keys_subset = range(len(keys_subset))
    i_shared_keys = range(len(shared_keys))
    maxdistance = global_maxdistance
    num_combinations = len(keys_subset) * len(shared_keys)
    progress_step = max(1, int(num_combinations/250))
    worker = mp.current_process()
    start = time.time()
    
    distances_fast = compute_distances(keys_subset,shared_keys)

    for i, (i_key1, i_key2) in enumerate(itertools.product(i_keys_subset, i_shared_keys)):
        key1 = keys_subset[i_key1]
        key2 = shared_keys[i_key2]

        if ((maxdistance is None) or (distances_fast[i_key1][i_key2]<=maxdistance)):

            # Get the time series data for both keys
            data1 = shared_timeseries[key1]
            data2 = shared_timeseries[key2]

            # Align the data on common index values
            common_index = list(set(data1.keys()) & set(data2.keys()))
            if len(common_index) >= MIN_COMMON_STATION_VALUES:
                aligned_data1 = [data1[key] for key in common_index]
                aligned_data2 = [data2[key] for key in common_index]

                # Calculate the Pearson correlation coefficient between the two columns using NumPy functions
                correlation = np.corrcoef(aligned_data1, aligned_data2)[0, 1]
                if not math.isnan(correlation):
                    # Add the correlation to the dictionary
                    if key1 not in correlations:
                        correlations[key1] = {}
                        distances[key1] = {}
                    correlations[key1][key2] = correlation
                    distances[key1][key2] = distances_fast[i_key1][i_key2]

        # Print progress message
        if i % progress_step == 0:
            print('\rCalculating Correlations and Distances: ({}):{:.2f}% {} out of {}'.format(worker.name, i * 100. / num_combinations, i, num_combinations), end="")
    # #test out of memory
    # if worker.name=='ForkPoolWorker-2':
    #     a=np.empty(1000*1000*1000*10)

    print("\n")
    print('Total time (sec): {} worker: {} keys: {}\n'.format(time.time() - start, worker.name, len(i_keys_subset)))
    correlation_file = 'temp_correlation_file_{}_{}.json'.format(worker.name,start)
    distance_file = 'temp_distance_file_{}_{}.json'.format(worker.name,start)
    write_correlations_and_distances(keys_subset, correlations, distances, correlation_file, distance_file)
    return correlations, distances

# initialize worker processes
def init_worker(timeseries, keys, maxdistance):
    # declare scope of a new global variable
    global shared_timeseries, shared_keys, global_maxdistance
    # store argument in the global variable for this process
    shared_timeseries = timeseries
    shared_keys = keys
    global_maxdistance = maxdistance
    worker = mp.current_process()
    print("\n")
    print('Worker {} started'.format(worker.name), end="")
    print("\n")

def map_with_exception(keys_subset):
    try:
        correlations, distances = compute_correlations(keys_subset)
        return correlations, distances, None
    except Exception as e:
        return None, None, e
    
def calc_correlations_and_distances_parallel(timeseries, start, end, maxdistance):
    # Create empty dictionaries to hold the correlation values for each key/timeseries and distances between each key
    correlations = {}
    distances = {}

    # Get keys
    #keys = list(timeseries.keys())[:2000]
    keys = list(timeseries.keys())

    # Compute correlations in parallel
    num_processes = min(mp.cpu_count(), len(keys))
    # limit process to 7
    num_processes = min(6, num_processes)
    print('Evaluating Correlation in parallel in range ({},{}), num_processes={}'.format(start,end,num_processes))
    pool = mp.Pool(processes=num_processes, initializer=init_worker, initargs=(timeseries,keys,maxdistance))
    key_sub=list(timeseries.keys())
    if (end is not None):
        key_sub=key_sub[:end]
    if (start is not None):
        key_sub=key_sub[start:]

    keys_batches = np.array_split(key_sub, num_processes)
    # results = []
    # for batch in keys_batches:
    #     async_results = [pool.apply_async(map_with_exception, (compute_correlations, [batch]))]
    # # wait for tasks to complete
    # results += [res.get() for res in async_results]
    #results = pool.map(map_with_exception, keys_batches)

    def error_callback(e):
        print("Exception in worker process:", e)
        pool.terminate()
        pool.join()
        raise e
    
    try:
        results = pool.map_async(map_with_exception, keys_batches, error_callback=error_callback).get()
    except Exception as e:
        error_callback(e)

    # Merge dictionaries from processes
    for result in results:
        exception = result[2]
        if exception is not None:
            # handle the exception here, e.g. raise an error or log it
            print('\nException reised {}'.format(exception))
            #raise exception
        if result[0]!=None:
            for key, values in result[0].items():
                key1=key[:]
                for key2, value in values.items():
                    if key1 not in correlations:
                        correlations[key1] = {}
                    correlations[key1][key2] = value
        if result[1]!=None:
            for key, values in result[1].items():
                key1=key[:]
                for key2, value in values.items():
                    if key1 not in distances:
                        distances[key1] = {}
                    distances[key1][key2] = value

    # Close and join the worker pool
    pool.close()
    pool.join()
    
    return correlations, distances

def calc_correlations_and_distances(timeseries,start,end,maxdistance):
    # Get up to 200 unique keys/timeseries in the dataset for debugging purposes
    #keys = list(timeseries.keys())[:2000]
    keys = list(timeseries.keys())
    # Create empty dictionaries to hold the correlation values for each key/timeseries and distances between each key
    correlations = {}
    distances = {}

    # Loop through all the combinations of keys/timeseries
    str_range = ''
    if (start is not None) or (end is not None):
        str_range='({},{})'.format(start,end)
        keys_subset=list(timeseries.keys())
        if (end is not None):
            keys_subset=keys_subset[:end]
        if (start is not None):
            keys_subset=keys_subset[start:]
        i_subset = range(len(keys_subset))
        i_keys = range(len(keys))
        loop=itertools.product(i_subset, i_keys)
        num_combinations = len(list(itertools.product(i_subset, i_keys)))
    else:
        i_subset = range(len(keys_subset))
        i_keys = range(len(shared_keys))
        loop=itertools.combinations(i_keys, 2)
        num_combinations = len(list(itertools.combinations(i_keys, 2)))

    distances_fast = compute_distances(keys_subset,keys)

    progress_step = max(1,int(num_combinations/250))
    for i, (i_key1, i_key2) in enumerate(loop):
        key1 = keys_subset[i_key1]
        key2 = keys[i_key2]

        if ((maxdistance is None) or (distances_fast[i_key1][i_key2]<=maxdistance)):
            # Get the time series data for both keys
            data1 = timeseries[key1]
            data2 = timeseries[key2]
            # Align the data on common index values
            common_index = list(set(data1.keys()) & set(data2.keys()))
            if len(common_index) >= MIN_COMMON_STATION_VALUES:
                aligned_data1 = [data1[key] for key in common_index]
                aligned_data2 = [data2[key] for key in common_index]
                # Calculate the Pearson correlation coefficient between the two columns using NumPy functions
                correlation = np.corrcoef(aligned_data1, aligned_data2)[0, 1]     
                if not math.isnan(correlation):
                    # Add the correlation to the dictionary
                    if key1 not in correlations:
                        correlations[key1] = {}
                        distances[key1] = {}
                    # if key2 not in correlations:
                    #     correlations[key2] = {}
                    #     distances[key2] = {}
                    correlations[key1][key2] = correlation
                    # correlations[key2][key1] = correlation
                    distances[key1][key2] = distances_fast[i_key1][i_key2]
                    # distances[key2][key1] = distances_fast[i_key1][i_key2]

        # Print progress message
        if i % progress_step == 0:
            print('\rCalculating Correlations and Distances {}: {:.2f}%'.format(str_range, i * 100. / num_combinations), end="")
    print("\n")

    return correlations, distances

def exponential(x, CDD):
    return np.exp(-x / CDD)

def exponential_lsq(params, x, y):
    CDD = params[0]
    return np.exp(-x / CDD) - y

def fit_exp_linear(x, y):
    y = np.log(y)
    K, _ = np.polyfit(x, y, 1)
    CDD = -1/K
    return CDD

def calculate_cdd_parameters(timeseries, correlations, distances, start, end, maxdistance):
    # Initialize the CDD dictionary to store the fitted parameters
    CDDs = {}
    str_range=''
    if (timeseries is not None) and ((start is not None) or (end is not None)):
        str_range='({},{})'.format(start,end)
        keys_subset=list(timeseries.keys())
        if (end is not None):
            keys_subset=keys_subset[:end]
        if (start is not None):
            keys_subset=keys_subset[start:]
        loop=keys_subset
        num_keys = len(keys_subset)
    else:
        num_keys = len(correlations.keys())
        loop=correlations.keys()
    progress_step = max(1,int(num_keys/250))

    # Loop over each key in the timeseries dictionary and fit the data to an exponential curve
    for i, key in enumerate(loop):
        # Create a list of tuples containing distances and correlations for the given timeseries key
        key_distances = []
        for other_key in correlations.get(key, {}).keys():
            if other_key != key:
                distance_value = distances.get(key, {}).get(other_key, None)
                correlation_value = correlations.get(key, {}).get(other_key, None)
                if distance_value is not None and correlation_value is not None and \
                    ((maxdistance is None) or (distance_value<=maxdistance)):
                        key_distances.append((distance_value, correlation_value))
                    
        # Sort the list of tuples by distance
        key_distances = sorted(key_distances, key=lambda x: x[0])

        # Convert the key_distances list to arrays
        xdata = np.array([t[0] for t in key_distances])
        ydata = np.array([t[1] for t in key_distances])

        # Fit the data to an exponential curve using least squares
        if len(xdata)>0 and len(ydata)>0:
            # with open('testxydata.txt', "w") as f: 
            #     for x,y in zip(xdata,ydata): 
            #         f.write(f"{x}\t{y}\n")
            try:
                CDD_curve, _ = curve_fit(exponential, xdata, ydata, p0=CDD_INITIAL)
                residuals = ydata - exponential(xdata, *CDD_curve)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((ydata-np.mean(ydata))**2)
                r_squared = 1 - (ss_res / ss_tot)
                # CDD_curve, _ = least_squares(exponential_lsq, CDD_INITIAL, args=(xdata, ydata))
                # CDD_linear, _ = fit_exp_linear(xdata, ydata)
                CDD=CDD_curve[0]

                # num of station<CDD
                num_stations = np.sum(xdata<=CDD)
                # num of station<CDD and CORR> 1/e
                num_stations_corr = np.sum(np.logical_and(xdata<=CDD,ydata>=1/np.e))
            
                # Add the fitted parameter to the CDD dictionary
                CDDs[key] = [CDD, r_squared, num_stations, num_stations_corr]
            except Exception as e:
                print(f'\nException in curve_fit for key {key}\n')
        # Print progress message
        if i % progress_step == 0:
            print('\rCalculating CDD values {}: {:.2f}%'.format(str_range, i * 100. / num_keys), end="")
    print("\n")

    return CDDs

def plot_correlation_distance(correlations, distances, CDDs, key, maxdistance):
    import matplotlib.pyplot as plt

    # Create a list of tuples containing distances and correlations for the given timeseries key
    key_distances = []
    for other_key in correlations.get(key, {}).keys():
        if other_key != key:
            distance_value = distances.get(key, {}).get(other_key, None)
            correlation_value = correlations.get(key, {}).get(other_key, None)
            if distance_value is not None and correlation_value is not None and \
                ((maxdistance is None) or (distance_value<=maxdistance)):
                key_distances.append((distance_value, correlation_value, other_key))

    # Sort the list of tuples by distance
    key_distances = sorted(key_distances, key=lambda x: x[0])

    # Convert the key_distances list to arrays
    xdata = np.array([t[0] for t in key_distances])
    ydata = np.array([t[1] for t in key_distances])
    labelsdata = np.array([t[2] for t in key_distances])

    if len(xdata)>0 and len(ydata)>0:
        # get previously fitted value
        CDD = CDDs[key][0]
        rsquared = CDDs[key][1]

        # Create the scatter plot
        fig, ax = plt.subplots(figsize=(12,10))
        ax.scatter(xdata, ydata)
        ax.set_xlabel("Distances")
        ax.set_ylabel("Correlations")
        ax.set_title(f"Timeseries {key} Correlations vs. Distances\nCDD: {CDD:.3g}, r2:{rsquared:.3g}")

        # Plot the exponential curve
        xfit = np.linspace(0, max(xdata), 10000)
        yfit = exponential(xfit, CDD)
        ax.plot(xfit, yfit, "-r", label=f"Fit: y = e^(-x / {CDD:.3f})")

        # Add vertical and horizontal lines
        ax.axvline(x=CDD, color='red', linestyle='--')
        ax.axhline(y=1/np.e, color='red', linestyle='--')

        # # Label points with correlation value >0.9 or <-0.9
        for i, xy in enumerate(zip(xdata, ydata)):
            if abs(ydata[i]) > 0.9:
                ax.annotate(labelsdata[i], xy=xy, textcoords="data", fontsize=0.5, rotation=45)

        # Show legend and save plot to file
        ax.legend(loc="best")

        str_maxdistance=''
        if maxdistance is not None:
            str_maxdistance='_maxdist{}'.format(maxdistance)
        output_file = 'plot_{}{}_p0_{}.png'.format(str(key),str_maxdistance,CDD_INITIAL)
        fig.savefig(output_file, dpi=1200)

        # Close the plot to free up memory
        plt.close()
 
def read_correlations_and_distances(correlation_file, distance_file):
    # Initialize the correlations and distances dictionaries
    correlations = {}
    distances = {}

    # Read the correlation values from the file
    with open(correlation_file, "r") as f:
        correlations = json.load(f)

    # Read the distance values from the file
    with open(distance_file, "r") as f:
        distances = json.load(f)

    return correlations, distances

def write_correlations_and_distances(keys_subset, correlations, distances, correlation_file, distance_file):
    if keys_subset is not None:
        filtered_correlations = {key: correlations.get(key, {}) for key in keys_subset}
        filtered_distances = {key: distances.get(key, {}) for key in keys_subset}
    else:
        filtered_correlations = correlations
        filtered_distances = distances
    # Write the correlation values to file
    with open(correlation_file, "w") as f:
        json.dump(filtered_correlations, f)
    # Write the distance values to file
    with open(distance_file, "w") as f:
        json.dump(filtered_distances, f)

def write_conflicted_timeseries(conflicted_timeseries, conflicted_timeseries_file):
    # Write the conflicted_timeseries values to file
    with open(conflicted_timeseries_file, "w") as f:
        json.dump(conflicted_timeseries, f)

def write_cdd_to_file(CDD, file_path):
    with open(file_path, "w") as f:
        f.write("x\ty\tCDD\trsquared\tnum <=CDD\tnum <=CDD AND >=1/e\n")
        for key, value in CDD.items():
            x, y = key.split("_")
            f.write(f"{x}\t{y}\t{value[0]}\t{value[1]}\t{value[2]}\t{value[3]}\n")

def read_cdd_from_file(file_path):
    CDD = {}
    with open(file_path, "r") as f:
        # Skip the header line
        next(f)
        # Read the data from the file
        for line in f:
            x, y, value, rsquared, num_stations, num_stations_corr = line.strip().split("\t")
            key = f"{x}_{y}"
            CDD[key] = [float(value), float(rsquared), int(num_stations), int(num_stations_corr)]
    return CDD

def read_station_list_from_file(file_path):
    StationKeys = {}
    with open(file_path, "r") as f:
        # Skip the header line
        next(f)
        # Read the data from the file
        for line in f:
            x, y, value = line.strip().split("\t")
            x = format(float(x), '.10f')
            y = format(float(y), '.10f')
            key = f"{x}_{y}"
            StationKeys[key] = float(value)
    return StationKeys

def read_timeseries_keys_from_file(filename):
    keys = []
    with open(filename, "r") as f:
        for line in f:
            key = line.strip()
            if len(key) > 0:
                x, y = key.split('_')
                x = format(float(x), '.10f')
                y = format(float(y), '.10f')
                new_key = f'{x}_{y}'
                keys.append(new_key)

    return keys

# def apply_precision_to_keys(dictionary):
#     new_dict = {}
#     for main_key, sub_dict in dictionary.items():
#         x, y = main_key.split('_')
#         x = format(float(x), '.10f')
#         y = format(float(y), '.10f')
#         new_main_key = f'{x}_{y}'
#         new_sub_dict = {}
#         for sub_key, value in sub_dict.items():
#             x, y = sub_key.split('_')
#             x = format(float(x), '.10f')
#             y = format(float(y), '.10f')
#             new_sub_key = f'{x}_{y}'
#             new_sub_dict[new_sub_key] = value
#         new_dict[new_main_key] = new_sub_dict
#     return new_dict

def apply_precision_to_keys(dictionary):
    new_dict = {}
    for main_key, sub_dict in dictionary.items():
        x, y = map(lambda val: format(float(val), '.10f'), main_key.split('_'))
        new_main_key = f'{x}_{y}'
        new_sub_dict = {f'{format(float(sub_key.split("_")[0]), ".10f")}_{format(float(sub_key.split("_")[1]), ".10f")}': value for sub_key, value in sub_dict.items()}
        new_dict[new_main_key] = new_sub_dict
    return new_dict

def apply_precision_to_keys_in_place(dictionary):
    keys_to_update = list(dictionary.keys())
    for main_key in keys_to_update:
        sub_dict = dictionary.pop(main_key)
        x, y = map(lambda val: format(float(val), '.10f'), main_key.split('_'))
        new_main_key = f'{x}_{y}'
        new_sub_dict = {f'{format(float(sub_key.split("_")[0]), ".10f")}_{format(float(sub_key.split("_")[1]), ".10f")}': value for sub_key, value in sub_dict.items()}
        dictionary[new_main_key] = new_sub_dict

def write_timeseries_to_files(keys, timeseries):
    for key in keys:
        if (key in timeseries.keys()):
            filename = key + "_timeserie.txt"
            with open(filename, "w") as f:
                for k, v in sorted(timeseries[key].items()):
                    f.write(str(k) + "\t" + str(v) + "\n")

def main():
    # Get the input directory from the command line argument
    if len(sys.argv) < 2:
        print("Usage: python cdd.py <directory>/<--analyze>/<--merge-and-filter-jsons>/<--generatemap>")
        print("[--start <start station> --end <end station> --parallel --only-extract-timeseries <timeseries_keys_file>] --maxdistance <max distance in km>")
        sys.exit(1)
    base_dir = sys.argv[1]
    # Parse the optional parameters
    start = None
    end = None
    parallel = False
    only_extract_timeseries_file = None
    maxdistance = None
    for i in range(2, len(sys.argv)):
        if sys.argv[i] == '--start':
            start = int(sys.argv[i+1])
        elif sys.argv[i] == '--end':
            end = int(sys.argv[i+1])
        elif sys.argv[i] == '--parallel':
            parallel = True
        elif sys.argv[i] == '--exclude-stations':
            exclude_stations_file = sys.argv[i+1]        
        elif sys.argv[i] == '--only-extract-timeseries':
            only_extract_timeseries_file = sys.argv[i+1]
        elif sys.argv[i] == '--maxdistance':
            maxdistance = float(sys.argv[i+1])

    # Read the txt files and construct the timeseries data and conflicted timeseries data
    timeseries, conflicted_timeseries = read_files(base_dir)
    correlations = {}
    distances = {}

    if (only_extract_timeseries_file is not None):
        keys_to_extract = read_timeseries_keys_from_file(only_extract_timeseries_file)
        write_timeseries_to_files(keys_to_extract,timeseries)
    else:
        correlation_file = "temp_correlation_file.json"
        distance_file = "temp_distance_file.json"
        conflicted_timeseries_file = "conflicted_timeseries_file.json"
        cdd_file = 'cdd_file_p0{}.txt'.format(CDD_INITIAL)
        str_parallel=''
        str_maxdistance=''
        if parallel:
            str_parallel='_parallel'
        if maxdistance is not None:
            str_maxdistance='_maxdist{}'.format(maxdistance)
        if (start is not None) or (end is not None) or (maxdistance is not None):
            cdd_file = 'cdd_file_{}_{}{}{}_p0{}.txt'.format(start,end,str_parallel,str_maxdistance,CDD_INITIAL)
            distance_file = 'temp_distance_file_{}_{}{}{}.json'.format(start,end,str_parallel,str_maxdistance)
            correlation_file = 'temp_correlation_file_{}_{}{}{}.json'.format(start,end,str_parallel,str_maxdistance)
            conflicted_timeseries_file = 'conflicted_timeseries_file.json_{}_{}{}{}.json'.format(start,end,str_parallel,str_maxdistance)
        keys_subset = None
        if base_dir == '--analyze':
            # use the read function to read correlations and distances instead of caltulating it again
            if os.path.exists(correlation_file) and os.path.exists(distance_file):
                correlations, distances = read_correlations_and_distances(correlation_file, distance_file)
            else:
                distance_file = 'temp_distance_file_{}_{}{}.json'.format(start,end,str_parallel)
                correlation_file = 'temp_correlation_file_{}_{}{}.json'.format(start,end,str_parallel)
                correlations, distances = read_correlations_and_distances(correlation_file, distance_file)

            if os.path.exists(cdd_file):
                CDDs = read_cdd_from_file(cdd_file)
            else:
                CDDs = calculate_cdd_parameters(None, correlations, distances, start, end, maxdistance)
                write_cdd_to_file(CDDs, cdd_file)
        elif base_dir == "--merge-and-filter-jsons":
            # Initialize the dictionary variable
            start = time.time()
            intermediate=start

            print(f"Starting merge...\n")

            # Iterate over all the files in the folder
            folder_path = '.'
            for file_name in os.listdir(folder_path):
                if file_name.startswith("temp_correlation_file") and file_name.endswith(".json"):
                    # Get the full correlation file path
                    correlation_file_path = os.path.join(folder_path, file_name)

                    # Get the corresponding distance file name
                    distance_file_name = file_name.replace("temp_correlation_file", "temp_distance_file")

                    # Get the full distance file path
                    distance_file_path = os.path.join(folder_path, distance_file_name)

                    print(f"Loading files {file_name} and {distance_file_name}...")
                    # Call the function and get the correlations and distances
                    file_correlations, file_distances = read_correlations_and_distances(correlation_file_path, distance_file_path)
                    print(f"Done (time={time.time()-intermediate}\n")
                    intermediate=time.time()

                    print(f"Adding data (len={len(file_correlations)})")
                    # Append the correlations to the correlations dictionary
                    correlations.update(file_correlations)
                    print(f"Done (time={time.time()-intermediate}\n")
                    intermediate=time.time()

                    # Append the distances to the distances dictionary
                    distances.update(file_distances)
                    print(f"Current key lenghts correlations={len(correlations)} distances={len(distances)}\n")

            if exclude_stations_file is not None:

                print(f"Loading keys to remove: {exclude_stations_file}")
                if os.path.exists(exclude_stations_file):
                    keys_to_remove = read_station_list_from_file(exclude_stations_file)
                print(f"Done (time={time.time()-intermediate}\n")
                intermediate=time.time()

                print(f"Keys to remove: {len(keys_to_remove)}\n")
                

                print(f"Appling same precision to correlations keys...")
                #apply_precision_to_keys_in_place(correlations)
                correlations = apply_precision_to_keys(correlations)
                print(f"Done (time={time.time()-intermediate}\n")
                intermediate=time.time()

                print(f"Appling same precision to distances keys...")
                #apply_precision_to_keys_in_place(correlations)
                distances = apply_precision_to_keys(distances)
                print(f"Done (time={time.time()-intermediate}\n")
                intermediate=time.time()
                
                # # Remove the keys from the correlations dictionary
                # for key_main in keys_to_remove.keys():
                #     if key_main in correlations:
                #         del correlations[key_main]
                #         del distances[key_main]

                # for key_main in correlations.keys():
                #     sub_dict = correlations[key_main]
                #     for key_sub in keys_to_remove:
                #         if key_sub in sub_dict:
                #             del sub_dict[key_sub]

                # for key_main in distances.keys():
                #     sub_dict = distances[key_main]
                #     for key_sub in keys_to_remove:
                #         if key_sub in sub_dict:
                #             del sub_dict[key_sub]

                print(f"Removing correlation keys...")
                # Remove the keys from the correlations dictionary
                correlations = {key_main: {key_sub: value for key_sub, value in sub_dict.items() if key_sub not in keys_to_remove}
                            for key_main, sub_dict in correlations.items() if key_main not in keys_to_remove}
                print(f"Done (time={time.time()-intermediate}\n")
                intermediate=time.time()

                print(f"Removing distances keys...")
                # Remove the keys from the distances dictionary
                distances = {key_main: {key_sub: value for key_sub, value in sub_dict.items() if key_sub not in keys_to_remove}
                            for key_main, sub_dict in distances.items() if key_main not in keys_to_remove}
                print(f"Done (time={time.time()-intermediate}\n")
                intermediate=time.time()
                print(f"Current key lenghts after filtering: correlations={len(correlations)} distances={len(distances)}\n")


            print(f"Calculating CDDs...")
            CDDs = calculate_cdd_parameters(None, correlations, distances, start, end, maxdistance)
            print(f"Done (time={time.time()-intermediate}\n")
            intermediate=time.time()
            print(f"Wrinting CDDs to file...")
            write_cdd_to_file(CDDs, cdd_file)
            print(f"Done (time={time.time()-intermediate}\n")
            print(f"Total time={time.time()-start}\n")
        elif base_dir == '--generatemap':
                if os.path.exists(cdd_file):
                    CDDs = read_cdd_from_file(cdd_file)
                    print(f"\r{len(CDDs)} keys loaded")
                    cdd_map_file = 'CDDmap_all_pr.nc'
                    print(f"\rGenerating CDD map File {cdd_map_file}...")

                    target_coords =  LatLong('template_European_01min.nc', 'template_European_01min.nc')
                    grid_x, grid_y = target_coords.lons, target_coords.lats
                    mv_target = target_coords.mv
                    if DEBUG_CDDMAP:
                        #write Debug Template Map
                        template_nc =  NetCDFReader('template_European_01min.nc')
                        template_nc_values = template_nc.values

                        # target_lats=target_lats[1800-(9*20):1800-(-16*20), 3600+(-15*20):3600+(16*20)]
                        # target_lons=target_lons[1800-(9*20):1800-(-16*20), 3600+(-15*20):3600+(16*20)]
                        if target_coords.lats.shape==(3600,7200):
                            # Global_3arcmin DEBUG
                            grid_y=target_coords.lats[1800-int(DEBUG_MAX_LAT*20):1800-int(DEBUG_MIN_LAT*20), 3600+int(DEBUG_MIN_LON*20):3600+int(DEBUG_MAX_LON*20)]
                            grid_x=target_coords.lons[1800-int(DEBUG_MAX_LAT*20):1800-int(DEBUG_MIN_LAT*20), 3600+int(DEBUG_MIN_LON*20):3600+int(DEBUG_MAX_LON*20)]
                            #latefas-=0.008369999999992217
                            # lonefas-=0.00851999999999431
                            #lonefas-=0.024519999999977227   
                            template_nc_values=template_nc_values[1800-int(DEBUG_MAX_LAT*20):1800-int(DEBUG_MIN_LAT*20), 3600+int(DEBUG_MIN_LON*20):3600+int(DEBUG_MAX_LON*20)]
                        else:
                            # European_1arcmin DEBUG
                            selection_lats = np.logical_and(target_coords.lats[:,0]>=DEBUG_MIN_LAT,target_coords.lats[:,0]<=DEBUG_MAX_LAT)
                            selection_lons = np.logical_and(target_coords.lons[0,:]>=DEBUG_MIN_LON,target_coords.lons[0,:]<=DEBUG_MAX_LON)
                            grid_y=target_coords.lats[selection_lats,:][:,selection_lons]
                            grid_x=target_coords.lons[selection_lats,:][:,selection_lons]
                            template_nc_values=template_nc_values[selection_lats,:][:,selection_lons]

                        template_nc_values = np.array(template_nc_values, dtype=np.int8)
                        metadata = {'variable': {'shortname': 'template', 'compression': 9, 'least_significant_digit': 5, 'mv': np.int8(template_nc.mv)},
                                    'geographical': {'variable_y_name': 'lat', 'variable_x_name': 'lon'},
                                    'dtype': template_nc_values.dtype,
                                    'rows': grid_y[:,0].size, 'cols': grid_x[0,:].size,
                                    'lats': grid_y[:,0], 'lons': grid_x[0,:],
                                    }

                        templace_writer = NetCDFWriter("template_test.nc",False,**metadata)
                        templace_writer.add_to_stack(template_nc_values, None)
                        templace_writer.finalize()
                        print('\tTemplace Map test Generation completed')


                    print('Starting interpolation...')
                    xp = []
                    yp = []
                    values = []
                    for key, value in CDDs.items():
                        x, y = key.split("_")
                        z = value[0]
                        xp.append(float(x))
                        yp.append(float(y))
                        values.append(float(z))
                    xp = np.array(xp,dtype='float64')
                    yp = np.array(yp,dtype='float64')
                    values = np.array(values,dtype='float64')
                    mv_source = 0
                    grid_details = {'gridType': 'NOT rotated', 'Nj': 1, 'radius': 1.0}
                    mode='adw'
                    nnear=11
                    # mode='nearest'
                    # nnear=1
                    scipy_interpolation = ScipyInterpolation(xp, yp, grid_details, values,
                                                            nnear, mv_target, mv_source,
                                                            target_is_rotated=False, parallel=True,
                                                            mode=mode)
                    results, weights, indexes = scipy_interpolation.interpolate(grid_x, grid_y)
                    print('Interpolation completed.')
                    grid_data = results.reshape(grid_x.shape)
                    metadata = {'variable': {'shortname': 'CDD', 'compression': 9, 'least_significant_digit': 5, 'mv': mv_target},
                                'geographical': {'variable_y_name': 'lat', 'variable_x_name': 'lon'},
                                'dtype': grid_data.dtype,
                                'rows': grid_y[:,0].size, 'cols': grid_x[0,:].size,
                                'lats': grid_y[:,0], 'lons': grid_x[0,:],
                                }

                    output_writer = NetCDFWriter(cdd_map_file,False,**metadata)
                    output_writer.add_to_stack(grid_data, None)
                    output_writer.finalize()
                    print('\nCDD Map Generation completed')
                else:
                    print(f"\rError: CDD File {cdd_file} not found.")
        else:
            if len(timeseries.keys())>0:
                if (start is not None) or (end is not None):
                    keys_subset=list(timeseries.keys())
                    print(f"\rTotal number of keys is: {len(keys_subset)}")
                    if (end is not None):
                        keys_subset=keys_subset[:end]
                    if (start is not None):
                        keys_subset=keys_subset[start:]
                    print(f"\rNow anaylzing keys from {start} to {end}: {len(keys_subset)} keys")

                # Calculate the correlation values and distances between each key
                if parallel:
                    correlations, distances = calc_correlations_and_distances_parallel(timeseries, start, end, maxdistance)
                else:
                    correlations, distances = calc_correlations_and_distances(timeseries, start, end, maxdistance)

                write_correlations_and_distances(keys_subset, correlations, distances, correlation_file, distance_file)
                write_conflicted_timeseries(conflicted_timeseries, conflicted_timeseries_file)

                CDDs = calculate_cdd_parameters(timeseries, correlations, distances, start, end, maxdistance)
                write_cdd_to_file(CDDs, cdd_file)
            else:
                print(f"\rNo timeserie found in folder {base_dir}")

        if base_dir != '--generatemap':
            # Print the resulting correlations dictionary and distances dictionary
            # print("Correlations:")
            # print(correlations)
            # print("\nDistances:")
            # print(distances)

            # Print the resulting conflicted timeseries dictionary
            # print("\nConflicted Timeseries:")
            # print(conflicted_timeseries)

            # Print the resulting CDDs
            # print("\nCDDs:")
            # print(CDDs)
            if keys_subset is not None:
                loop=keys_subset
            else:
                loop=correlations.keys()

            # Save plots for debugging  
            if DEBUG_CDDMAP:
                # plot the graph for the CDD
                for i, key in enumerate(loop):
                    if i>=10: #save only 10 plots
                        break
                    print(f"\rSaving plots: {i+1}", end="")
                    if key in correlations:
                        plot_correlation_distance(correlations, distances, CDDs, key, maxdistance)
            print("\n")

if __name__ == "__main__":
    main()