
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
import xarray as xr
from osgeo import gdal
import netCDF4 as nc


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


def get_unpacking_metadata(file_interpolated_values: Path) -> tuple[float, float, float]:
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


def get_interpolated_values(file_interpolated_values: Path, is_compressed_data: bool) -> np.ndarray:
    # values_ds = rxr.open_rasterio(file_interpolated_values, mask_and_scale=True)
    scale_factor, offset, no_data = get_unpacking_metadata(file_interpolated_values)
    values_ds = rxr.open_rasterio(file_interpolated_values)
    values_ds = values_ds.sel(band=1)
    values_ds = values_ds.rename({'band': 'values'})
    values_ds = values_ds.sortby('y', ascending=False)
    values = unpack_data(values_ds.values, scale_factor, offset, no_data)
    return values


def get_pixelarea_catchments_and_interpolated_values(pixelarea_file_path: Path, catchments_file_path: Path,
                                                     file_interpolated_values: Path, is_compressed_data: bool) -> xr.Dataset:
    ds_pixelarea = xr.open_dataset(pixelarea_file_path)
    ds_pixelarea = ds_pixelarea.sortby('lat', ascending=False)
    ds_pixelarea = ds_pixelarea.rename({"Band1": "pixarea"})
    pixelarea = ds_pixelarea.pixarea
    pixelarea_values = pixelarea.values

    ds_catchments = xr.open_dataset(catchments_file_path)
    ds_catchments = ds_catchments.sortby('lat', ascending=False)
    ds_catchments = ds_catchments.rename({"Band1": "catchments"})
    catchments = ds_catchments.catchments
    catchments_values = catchments.values

    values = get_interpolated_values(file_interpolated_values, is_compressed_data)

    combined_data = xr.Dataset({'pixarea': (['lat', 'lon'], pixelarea_values),
                               'catchments': (['lat', 'lon'], catchments_values),
                               'values': (['lat', 'lon'], values)},
                               coords={'lon': ds_pixelarea['lon'], 'lat': ds_pixelarea['lat']})
    combined_data = combined_data.where(combined_data.catchments != 0)
    return combined_data


def write_results(outfile: Path, df: pd.DataFrame):
    if outfile.is_file():
        df_from_file = pd.read_csv(outfile, delimiter='\t', index_col='catchments')
        df = df.drop('pixarea', axis=1)
        df = pd.merge(df_from_file, df, on=['catchments'])
    df.to_csv(outfile, sep='\t', index=True)


def run(file_catchments: str, file_interpolated_values: str, file_pixel_area: str, outfile: str, is_compressed_data: bool):
    file_catchments_path = Path(file_catchments)
    file_interpolated_values_path = Path(file_interpolated_values)
    file_pixel_area_path = Path(file_pixel_area)
    outfile_path = Path(outfile)
    ds = get_pixelarea_catchments_and_interpolated_values(file_pixel_area_path, file_catchments_path,
                                                          file_interpolated_values_path, is_compressed_data)

    # Convert the precipitation values from mm to m
    ds['values'] /= 1000.0
    # Convert precipitation from m to m3
    ds['values'] = ds['values'] * ds['pixarea']

    values_total_group = ds.groupby('catchments').sum()

    values_total_group_df = values_total_group.to_dataframe()
#     values_total_group_df = values_total_group_df.dropna()

    test_code = file_interpolated_values_path.parent.name
    # test_code = 'ALL'
    values_total_group_df = values_total_group_df.rename(columns={'values': f'{test_code}'})
    write_results(outfile_path, values_total_group_df)


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

#     try:
    # setup option parser
    parser = ArgumentParser(epilog=program_license, description=program_version_string+program_longdesc)

    # set defaults
    parser.set_defaults(is_compressed_data=False)

    parser.add_argument("-t", "--catchments", dest="file_catchments", required=True, type=FileUtils.file_type,
                        help="Set file path for catchments netCDF.",
                        metavar="txt_file")
    parser.add_argument("-i", "--interpolated", dest="file_interpolated_values", required=True, type=FileUtils.file_type,
                        help="Set file path for interpolated values in tiff format.",
                        metavar="tiff_file")
    parser.add_argument("-p", "--pixelarea", dest="file_pixel_area", required=True, type=FileUtils.file_type,
                        help="Set file path for pixel area netCDF.",
                        metavar="pixel_area_file")
    parser.add_argument("-o", "--out", dest="outfile", required=True, type=FileUtils.file_or_folder,
                        help="Set output file containing the statistics.",
                        metavar="output_file")
    parser.add_argument("-c", "--compressed", dest="is_compressed_data", action="store_true",
                        help="Indicates if data in the tiff file is compressed and decompresses it before using. [default: %(default)s]")

    # process options
    args = parser.parse_args(argv)

    print(f"Catchments File: {args.file_catchments}")
    print(f"Pixel Area File: {args.file_pixel_area}")
    print(f"Interpolated values: {args.file_interpolated_values}")
    print(f"Output File: {args.outfile}")
    print(f"Compressed Data: {args.is_compressed_data}")

    run(args.file_catchments, args.file_interpolated_values, args.file_pixel_area, args.outfile, args.is_compressed_data)
    print("Finished.")
#     except Exception as e:
#         indent = len(program_name) * " "
#         sys.stderr.write(program_name + ": " + repr(e) + "\n")
#         sys.stderr.write(indent + "  for help use --help")
#         return 2


def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_script()
