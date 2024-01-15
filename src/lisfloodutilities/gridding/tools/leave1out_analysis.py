__author__="Goncalo Gomes"
__date__="$Jul 12, 2023 12:01:00$"
__version__="0.1"
__updated__="$Jul 13, 2023 10:41:00$"

"""
Copyright 2019-2020 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import sys
import os
import csv
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError
from lisfloodutilities.gridding.lib.utils import FileUtils
import pandas as pd
import numpy as np
import rioxarray as rxr
from osgeo import gdal
import netCDF4 as nc
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.feature_selection import r_regression
from typing import List, Tuple


def mean_bias_error(y_true: np.array, y_pred: np.array) -> float:
    '''
    Parameters:
        y_true (array): Array of observed values
        y_pred (array): Array of prediction values

    Returns:
        mbe (float): Bias score
    '''
    mbe = np.mean(y_true - y_pred)
    return mbe


def critical_success_index(df: pd.DataFrame, limit_value: float) -> float:
    '''
    Critical Success Index is used only for Precipitation.

    Parameters:
        df (pandas dataframe): Dataframe containing both observed and predicted values
        limit_value (float): Limit in millimeters to be used in the formula

    csi = wc / (wc+wi+di) 
    Wc = total number of points where y_true and y_pred > 1mm
    Wi = total number of points where y_true  > 1mm and y_pred < 1mm
    Di = total number of points where y_true  < 1mm and y_pred > 1mm

    Returns:
        csi (float): Critical Success Index
    '''
    wc = len(df[(df['true_value'] > limit_value) & (df['interpolated_value'] > limit_value)])
    wi = len(df[(df['true_value'] > limit_value) & (df['interpolated_value'] < limit_value)])
    di = len(df[(df['true_value'] < limit_value) & (df['interpolated_value'] > limit_value)])
    csi = wc / (wc + wi + di) 
    return csi


def get_netcdf_meteo_variable_name(nc_file_obj):
    # Only one variable must be present in netcdf files
    num_dims = 3 if 'time' in nc_file_obj.variables else 2
    varname = [v for v in nc_file_obj.variables if len(nc_file_obj.variables[v].dimensions) == num_dims][0]
    return varname


def get_pixel_area(pixel_area_file: Path) -> np.ndarray:
    ncfile = nc.Dataset(pixel_area_file, 'r')
    var_name = get_netcdf_meteo_variable_name(ncfile)
    var = ncfile.variables[var_name]
    data = var[:]
    ncfile.close()
    return np.ascontiguousarray(data.flatten())


def write_results(outfile: Path, test_code: str, mae: float, mbe: float, pearson_r: float, values_sum: float,
                  count_pixels_1st_interval: float, count_pixels_2nd_interval: float, mse: float, csi: float):
    if outfile.is_file():
        df = pd.read_csv(outfile, delimiter='\t')
        new_data = [{
            'test': test_code,
            'mae': mae,
            'mbe': mbe,
            'mse': mse,
            'pearson_r': pearson_r,
            'csi': csi,
            'values_sum': values_sum,
            'pixels_with_values_0<X<=0.1': count_pixels_1st_interval,
            'pixels_with_values_0.1<X<1': count_pixels_2nd_interval,
        }]
        df = df.append(new_data, ignore_index=True)
    else:
        df = pd.DataFrame({
            'test': [test_code],
            'mae': [mae],
            'mbe': [mbe],
            'mse': [mse],
            'pearson_r': [pearson_r],
            'csi': [csi],
            'values_sum': [values_sum],
            'pixels_with_values_0<X<=0.1': [count_pixels_1st_interval],
            'pixels_with_values_0.1<X<1': [count_pixels_2nd_interval],
        })
    df.to_csv(outfile, sep='\t', index=False)


def unpack_data(data: np.ndarray, scale_factor: float = 1.0, offset: float = 0.0, no_data: float = -9999.0) -> np.ndarray:
    print(f'scale_factor: {scale_factor} offset: {offset} no_data: {no_data}')
    compressed_data = data.astype(float)
    compressed_data = np.where(compressed_data == no_data, np.nan, compressed_data)
    decompressed_data = compressed_data * scale_factor + offset
    decompressed_data = np.round(decompressed_data, 1)
    unpacked_no_data = no_data * scale_factor + offset
    decompressed_data = np.where(decompressed_data == unpacked_no_data, np.nan, decompressed_data)
    return decompressed_data


def format_metadata(metadata: dict, value: str, key: str, default: float) -> float:
    if value is None:
        return  float(metadata.get(key, default))
    return float(value)


def get_unpacking_metadata(file_interpolated_values: Path) -> Tuple[float, float, float]:
    ds=gdal.Open(file_interpolated_values.as_posix())
    gdal_nan = ds.GetRasterBand(1).GetNoDataValue()
    gdal_scale = ds.GetRasterBand(1).GetScale()
    gdal_offset = ds.GetRasterBand(1).GetOffset()
    print(f'gdal_nan: {gdal_nan} gdal_scale: {gdal_scale} gdal_offset: {gdal_offset}')
    metadata=ds.GetMetadata()
    print('metadata: ', metadata)
    scale_factor = format_metadata(metadata, gdal_scale, 'SCALE', 1.0)
    offset = format_metadata(metadata, gdal_offset, 'OFFSET', 0.0)
    no_data = format_metadata(metadata, gdal_nan, 'NODATA', -9999.0)
    ds = None
    return scale_factor, offset, no_data


def get_interpolated_values_dataframe(file_interpolated_values: Path) -> pd.DataFrame:
    scale_factor, offset, no_data = get_unpacking_metadata(file_interpolated_values)
    dataarray = rxr.open_rasterio(file_interpolated_values)
    band = dataarray[0]
    x, y, values = band.x.values, band.y.values, band.values
    x, y = np.meshgrid(x, y)
    x, y, values = x.flatten(), y.flatten(), unpack_data(values.flatten(), scale_factor, offset, no_data)
    df_interpolated_values = pd.DataFrame(
        {
            'x': x,
            'y': y,
            'value': values,
        }
    )
    # df_interpolated_values = df_interpolated_values.dropna()
    return df_interpolated_values


def get_interpolated_values(df_interpolated_values: pd.DataFrame, df_true_values: pd.DataFrame) -> np.ndarray:
    interpolated_values = np.ascontiguousarray(df_interpolated_values['value'].values)
    transformed_coordinates = (df_interpolated_values['x'].values, df_interpolated_values['y'].values)
    interpolated_values_query = cKDTree(data=np.vstack(transformed_coordinates).T, copy_data=True)
    stations = np.array([df_true_values['x'], df_true_values['y']]).T
    _, idx = interpolated_values_query.query(stations)
    return interpolated_values[idx]


def run(file_true_values: str, file_interpolated_values: str, file_pixel_area: str, outfile: str, run_csi: bool, limit_value: float):

    file_interpolated_values_path = Path(file_interpolated_values)
    outfile_path = Path(outfile)
    file_true_values_path = Path(file_true_values)
    # Load the CSV file into a pandas dataframe
    df_true_values = pd.read_csv(file_true_values_path, delimiter='\t', header=None, names=['x', 'y', 'true_value'])
    df_interpolated_values = get_interpolated_values_dataframe(file_interpolated_values_path)
    predicted_values_mm = np.ascontiguousarray(df_interpolated_values['value'].values)  # includes NaN
    df_interpolated_values = df_interpolated_values.dropna()
    interpolated_values = get_interpolated_values(df_interpolated_values, df_true_values)
    df_true_values['interpolated_value'] = interpolated_values
    # write intermediate results to a file
    df_true_values.to_csv(file_true_values_path.with_suffix('.tab'), sep="\t", index=False)

    csi = -9999.0
    if run_csi:
        csi = critical_success_index(df_true_values, 1.0)

    # process only values greater than 1.0 mm
    # limit_value = 1.0
    df_true_values = df_true_values[(df_true_values['true_value'] > limit_value) & (df_true_values['interpolated_value'] > limit_value)]

    true_values = df_true_values['true_value']
    predicted_values = df_true_values['interpolated_value']

    test_code = file_interpolated_values_path.parent.name
    mae = -9999.0
    mbe = -9999.0
    mse = -9999.0
    pearson_correlation_coeficient = -9999.0
    values_sum = -9999.0
    count_pixels_1st_interval = -9999.0
    count_pixels_2nd_interval = -9999.0
    
    if len(true_values) > 0 and len(predicted_values) > 0:
        mae = mean_absolute_error(true_values, predicted_values)
        mbe = mean_bias_error(true_values, predicted_values)
        mse = mean_squared_error(true_values, predicted_values)
        # The correlations needs at least 2 elements in each array
        if len(true_values) > 1 and len(predicted_values) > 1: 
            pearson_correlation_coeficient, p_value = pearsonr(true_values, predicted_values)
        # values_sum = np.round(np.nansum(np.ascontiguousarray(df_interpolated_values['value'].values)), 1)
        values_sum = np.round(np.nansum(np.ascontiguousarray(predicted_values)), 1)
    
        if len(file_pixel_area) > 0:
            pixel_area_m2 = get_pixel_area(Path(file_pixel_area))
            predicted_values_m = predicted_values_mm * 0.001
            predicted_values_m3 = predicted_values_m * pixel_area_m2
            values_sum = np.round(np.nansum(predicted_values_m3), 1)
            # number of pixels with values 0 < X <0.1 
            count_pixels_1st_interval = np.count_nonzero((~np.isnan(predicted_values_mm)) & (predicted_values_mm > 0) & (predicted_values_mm <= 0.1))
            # number of pixels with values 0.1 <= X < 1
            count_pixels_2nd_interval = np.count_nonzero((~np.isnan(predicted_values_mm)) & (predicted_values_mm > 0.1) & (predicted_values_mm < 1))

    # write the results to the output file
    write_results(outfile_path, test_code, mae, mbe, pearson_correlation_coeficient, values_sum,
                  count_pixels_1st_interval, count_pixels_2nd_interval, mse, csi)


def main(argv):
    '''Command line options.'''
    global quiet_mode

    program_name = os.path.basename(sys.argv[0])
    program_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    program_version = "v%s" % __version__
    program_build_date = "%s" % __updated__

    program_version_string = 'version %s (%s)\n' % (program_version, program_build_date)
    program_longdesc = '''
    This script uses RepeatedStratifiedKFold to generate several input files as training dataset where a small percentage
    of evenly distributed stations are removed to constitute the test cases. Later the training dataset will be interpolated
    with several algorithms: Spheremap, Inverse distance and Angular Distance Weighting. Then the error at the coordinates of
    the test cases is calculated by the difference between the interpolated value and the real value in the test case.
    '''
    program_license = """
    Copyright 2019-2020 European Union
    Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
    You may not use this work except in compliance with the Licence.
    You may obtain a copy of the Licence at:
    https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
    Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the Licence for the specific language governing permissions and limitations under the Licence.
    """

    try:
        # setup option parser
        parser = ArgumentParser(epilog=program_license, description=program_version_string+program_longdesc)

        # set defaults
        parser.set_defaults(run_csi=False, limit_value=1.0)

        parser.add_argument("-t", "--test", dest="file_true_values", required=True, type=FileUtils.file_type,
                            help="Set file path for true values form the test dataset in txt format.",
                            metavar="txt_file")
        parser.add_argument("-i", "--interpolated", dest="file_interpolated_values", required=True, type=FileUtils.file_type,
                            help="Set file path for interpolated values in tiff format.",
                            metavar="tiff_file")
        parser.add_argument("-p", "--pixelarea", dest="file_pixel_area", required=False, type=FileUtils.file_type,
                            help="Set file path for pixel area netCDF.",
                            metavar="pixel_area_file")
        parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_or_folder,
                            help="Set output file containing the statistics.",
                            metavar="output_file")
        parser.add_argument("-l", "--limit", dest="limit_value", required=False, type=float,
                            help="process only values greater than limit", metavar="limit_value")
        parser.add_argument("-c", "--csi", dest="run_csi", action="store_true",
                            help="Critical Success Index is used only for Precipitation. [default: %(default)s]")

        # process options
        args = parser.parse_args(argv)

        print(f"True values: {args.file_true_values}")
        print(f"Interpolated values: {args.file_interpolated_values}")
        print(f"Output File: {args.outfile}")
        print(f"Run CSI: {args.run_csi}")
        print(f"Limit Value: {args.limit_value}")
        file_pixel_area = ""
        if args.file_pixel_area is not None:
            file_pixel_area = args.file_pixel_area
            print(f"Pixel Area File: {file_pixel_area}")

        run(args.file_true_values, args.file_interpolated_values, file_pixel_area, args.outfile, args.run_csi, args.limit_value)
        print("Finished.")
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()