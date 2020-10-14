import os
import uuid

from osgeo import gdal

import numpy as np
import numpy.ma as ma


class PCRasterWriter:
    FORMAT = 'PCRaster'
    dtype_to_valuescale = {
        'uint8': (gdal.GDT_Byte, 'VS_BOOLEAN', np.int8),
        'int8': (gdal.GDT_Int32, 'VS_NOMINAL', np.int32),
        'int16': (gdal.GDT_Int32, 'VS_ORDINAL', np.int32),
        'float32': (gdal.GDT_Float32, 'VS_SCALAR', np.float32),
        'float64': (gdal.GDT_Float32, 'VS_SCALAR', np.float32),
    }

    def __init__(self, *args, **kwargs):
        """
        Constructor accept only keywords arguments:
        Must have either a clonemap or coordinates in input
        clonemap: path to a clone PCRaster map
        coordinates: dict of numpy 1D arrays representing coordinates, with 'x' and 'y' keys
        mv: Missing value (optional)
        """
        self._drv = gdal.GetDriverByName(self.FORMAT)
        self._drv.Register()
        self._clone_map = kwargs.get('clonemap')
        self._coordinates = kwargs.get('coordinates')
        mv = kwargs.get('mv')
        self.flipped_x = False
        self.flipped_y = False
        if not (self._clone_map or self._coordinates):
            raise ValueError('Need either a PCRaster clonemap or coordinates={"x": [...], "y": [...]}')

        if self._clone_map:
            src_ds = gdal.Open(self._clone_map.encode('utf-8'), gdal.GA_ReadOnly)
            src_band = src_ds.GetRasterBand(1)

            # Producing mask array
            self.cols = src_ds.RasterXSize
            self.rows = src_ds.RasterYSize
            self.src_geo_transform = src_ds.GetGeoTransform()

            rs = src_band.ReadAsArray()
            no_data_value = src_band.GetNoDataValue()
            src_ds = None
            del src_ds
            src_band = None
            del src_band
            rs = ma.masked_values(rs, no_data_value, copy=False)
            self._mask = ma.getmask(rs)
            rs = None
            self.mv = no_data_value if not mv else float(mv)
        else:
            self.cols = self._coordinates['x'].shape[0]
            self.rows = self._coordinates['y'].shape[0]
            self.src_geo_transform = self._get_geo_transform_from_coords()
            self.mv = np.nan if not mv else float(mv)
            self.flipped_x = self._coordinates['x'][0] > self._coordinates['x'][1]
            self.flipped_y = self._coordinates['y'][0] < self._coordinates['y'][1]

    def write(self, output_map_name, values, is_ldd=False):
        if is_ldd:
            etype = gdal.GDT_Byte
            valuescale = 'VS_LDD'
            as_type = np.int8
        else:
            etype, valuescale, as_type = self.get_valuescale(values)

        if not self._clone_map:
            # if no clone map was provided, we can calculate mask once we get values
            masked_values = ma.masked_values(values, self.mv, copy=False)
            self._mask = ma.getmask(masked_values)

        filled_values = self._mask_values(values, self.mv)
        filled_values = filled_values.astype(as_type)
        tmpfilename = 'tmp{}.map'.format(uuid.uuid4())
        dest_ds = self._drv.Create(tmpfilename, xsize=self.cols, ysize=self.rows, bands=1,
                                   eType=etype,
                                   options=["PCRASTER_VALUESCALE={}".format(valuescale)])

        # =============================================================================
        # Create a MEM clone of the dest file as Create returns readonly dataset
        # =============================================================================
        mem_drive = gdal.GetDriverByName('MEM')
        mem_drive.Register()
        mem_ds = mem_drive.CreateCopy('mem', dest_ds)
        dest_ds = None
        del dest_ds
        mem_ds.SetGeoTransform(self.src_geo_transform)
        rs = mem_ds.GetRasterBand(1)
        if self.flipped_y:
            # when we do not use a clone map, we need to flip up down values
            # not clear..it should be done always since PCRaster y coordinates are upside down.
            filled_values = np.flipud(filled_values)
        if self.flipped_x:
            # when we do not use a clone map, we need to flip up down values
            # not clear..it should be done always since PCRaster y coordinates are upside down.
            filled_values = np.fliplr(filled_values)
        rs.WriteArray(filled_values)
        # FIXME this set missing value statement is ignored once file is written on disk
        rs.SetNoDataValue(self.mv)
        rs.FlushCache()
        # =============== Write on disk
        out_ds = self._drv.CreateCopy(output_map_name.encode('utf-8'), mem_ds)
        out_ds = None
        del out_ds
        mem_ds = None
        del mem_ds
        os.unlink(tmpfilename)

    def _mask_values(self, values, mv=None):
        if isinstance(values, ma.core.MaskedArray):
            masked = ma.masked_where((self._mask | values.mask), values.data, copy=False)
        else:
            masked = ma.masked_where(self._mask, values, copy=False)
        mv = self.mv if not mv else mv
        masked = ma.filled(masked, mv)
        return masked

    def close(self):
        pass

    def get_valuescale(self, values):
        dtype = str(values.dtype)
        return self.dtype_to_valuescale[dtype]

    def _get_geo_transform_from_coords(self):
        # geotransform is just a tuple e.g.  (2500000.0, 5000.0, 0.0, 5500000.0, 0.0, -5000.0)
        # (-25.25, 0.016666666666666666, 0.0, 72.25, 0.0, -0.016666666666666666) geotransform with 1 arcmin europe maps
        # 0 - xul
        # 1 - w-e pixel resolution
        # 2 - 0.0 (?)
        # 3 - yul
        # 4 - 0.0
        # 5 - n-s pixel resolution (negative value)
        #
        top_left_x = np.round(self._coordinates['x'][0], 2)
        w_e_resolution = self._coordinates['x'][1] - self._coordinates['x'][0]

        top_left_y = np.round(self._coordinates['y'][0], 2)
        n_s_resolution = - (self._coordinates['y'][1] - self._coordinates['y'][0])
        return top_left_x, w_e_resolution, 0.0, top_left_y, 0.0, n_s_resolution
