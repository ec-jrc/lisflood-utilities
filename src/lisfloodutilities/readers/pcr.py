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

import os
import glob
import re
import sys

import numpy as np

try:
    from osgeo import gdal
    from osgeo.gdalconst import GA_ReadOnly
except ImportError as e:
    print("""
    [!] GDAL is not installed. 
    [!] Please, install GDAL binaries and libraries for your system and then install the relative GDAL Python interface.
    [!] If you are using a conda env, you can install with `conda install gdal`.
    [!] For a system installation in Debian/Ubuntu:
    [!] sudo apt-get install libgdal-dev libgdal
    [!] export CPLUS_INCLUDE_PATH=/usr/include/gdal
    [!] export C_INCLUDE_PATH=/usr/include/gdal
    [!] To know the gdal library version you have installed on your machine, run `gdal-config --version` 
    [!] e.g. 3.1.0 ==> pip install GDAL==3.1.0
    """)
    sys.exit(1)


class PCRasterMap:
    """
    A class representing a PCRaster map file.
    """

    def __init__(self, pcr_map_filename):
        dataset = gdal.Open(pcr_map_filename.encode('utf-8'), GA_ReadOnly)
        self.filename = pcr_map_filename
        self.dataset = dataset
        self.geo_transform = dataset.GetGeoTransform()
        self.cols = self.dataset.RasterXSize
        self.rows = self.dataset.RasterYSize
        self.band = self.dataset.GetRasterBand(1)

    def __eq__(self, other):
        res = (self.cols == other.cols and self.rows == other.rows
               and self.min == other.min and self.max == other.max
               and self.geo_transform == other.geo_transform)
        if not res:
            return False
        diff = np.abs(self.data - other.data)
        return np.allclose(diff, np.zeros(diff.shape), rtol=1.e-2, atol=1.e-3, equal_nan=True)

    def __ne__(self, other):
        return not self.__eq__(other)

    def close(self):
        """
        Release resources
        """
        self.band = None
        self.dataset = None

    @property
    def data(self):
        """
        :return: Values as 2D numpy array with shape (rows, cols)
        """
        data = self.band.ReadAsArray(0, 0, self.cols, self.rows)
        return data

    @property
    def mv(self):
        """
        :return: Value for missing data
        """
        return self.band.GetNoDataValue()

    @property
    def xul(self):
        """
        :return: top left x
        """
        return self.geo_transform[0]

    @property
    def yul(self):
        """
        :return: top left y
        """
        return self.geo_transform[3]

    @property
    def cell_length(self):
        """
        :return: w-e pixel resolution
        """
        return self.geo_transform[1]

    @property
    def cell_height(self):
        """
        :return: n-s pixel resolution (negative value)
        """
        return self.geo_transform[5]

    @property
    def min(self):
        """
        :return: Minimum value
        """
        return self.band.GetMinimum()

    @property
    def max(self):
        """
        :return: Maximum value
        """
        return self.band.GetMaximum()

    @property
    def lats(self):
        """
        :return: numpy 1D array representing lats
        """
        step = abs(self.cell_height)
        start = self.yul - step / 2
        end = self.yul - step * self.rows
        return np.arange(start, end, -step)

    @property
    def lons(self):
        """
        :return: numpy 1D array representing lons
        """
        step = self.cell_length
        start = self.xul + step / 2
        end = self.xul + step * self.cols
        return np.arange(start, end, step)

    @classmethod
    def build(cls, pcr_map_filename):
        """
        Class method to build a :class:`.PCRasterMap` object
        :param pcr_map_filename: path to PCRaster map file
        :return: PCRasterMap object
        """
        return cls(pcr_map_filename)


class PCRasterReader:
    """
    A class to read PCRaster maps dataset
    """

    FORMAT = 'PCRaster'
    digits = re.compile(r'\d+')

    def __init__(self, input_set):
        self.input_set = input_set
        self._driver = gdal.GetDriverByName(self.FORMAT)
        self._driver.Register()

    def get_metadata_from_set(self):
        """
        Get metadata from first map of a set.
        :return: A dictionary with keys ['rows', 'cols', 'lats', 'lons', 'dtype']
        """
        first_map, _ = next(self.fileset)
        return {'rows': first_map.rows, 'cols': first_map.cols,
                'lats': first_map.lats, 'lons': first_map.lons,
                'dtype': first_map.data.dtype}

    @property
    def fileset(self):
        """
        A generator that yield a tuple (PCRasterMap object, timestep int) at each iteration.

        :return: A generator object yielding a tuple of PCRasterMap object and an integer timestep
        """
        if self.input_is_single_file():
            yield PCRasterMap.build(self.input_set), None
        elif self.input_is_wildcard():
            filelist = sorted(glob.glob(self.input_set), key=self._extract_timestep)
            for f in filelist:
                yield PCRasterMap.build(f), self._extract_timestep(f)
        elif self.input_is_dir():
            filelist = sorted((f for f in os.listdir(self.input_set) if not f.endswith('.nc')),
                              key=self._extract_timestep)
            for f in filelist:
                f = os.path.join(self.input_set, f)
                yield PCRasterMap.build(f), self._extract_timestep(f)
        else:
            raise ValueError('Cannot guess input PCRaster files')

    def input_is_single_file(self):
        """
        Check if input_set is a single file
        :return: bool. True if input is a single file
        """
        return os.path.isfile(self.input_set)

    def input_is_wildcard(self):
        """
        Check if input_set is a Unix-like path expansion (e.g. /path/to/file*.0*)
        :return: bool. True if input is a path expansion
        """
        return '*' in self.input_set

    def input_is_dir(self):
        """
        Check if input_set is a directory
        :return: bool. True if input is a directory
        """
        return os.path.isdir(self.input_set)

    @classmethod
    def _extract_timestep(cls, f):
        filename = os.path.basename(f).replace('.', '')
        try:
            step = cls.digits.findall(filename)[0]
        except IndexError:
            return None
        else:
            return int(step)
