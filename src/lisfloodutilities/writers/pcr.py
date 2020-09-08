from osgeo import gdal

import numpy.ma as ma


class PCRasterWriter:
    FORMAT = 'PCRaster'

    def __init__(self, *args):
        self._clone_map = args[0]
        # =============================================================================
        # Create a MEM clone of the source file.
        # =============================================================================

        self._src_drv = gdal.GetDriverByName(self.FORMAT)
        self._src_drv.Register()
        self._src_ds = gdal.Open(self._clone_map.encode('utf-8'))
        self._src_band = self._src_ds.GetRasterBand(1)

        self._mem_ds = gdal.GetDriverByName('MEM').CreateCopy('mem', self._src_ds)

        # Producing mask array
        cols = self._src_ds.RasterXSize
        rows = self._src_ds.RasterYSize
        rs = self._src_band.ReadAsArray(0, 0, cols, rows)
        self.mv = self._src_band.GetNoDataValue()
        rs = ma.masked_values(rs, self.mv)
        self._mask = ma.getmask(rs)

    def write(self, output_map_name, values):
        drv = gdal.GetDriverByName(self.FORMAT)
        masked_values = self._mask_values(values)
        n = ma.count_masked(masked_values)
        self._mem_ds.GetRasterBand(1).SetNoDataValue(self.mv)
        self._mem_ds.GetRasterBand(1).WriteArray(masked_values)
        out_ds = drv.CreateCopy(output_map_name.encode('utf-8'), self._mem_ds)
        out_ds = None
        del out_ds

    def _mask_values(self, values):
        if isinstance(values, ma.core.MaskedArray):
            masked = ma.masked_where((self._mask | values.mask), values.data, copy=False)
        else:
            masked = ma.masked_where(self._mask, values, copy=False)
        masked = ma.filled(masked, self.mv)
        return masked

    def close(self):
        self._mem_ds = None
        self._src_ds = None
        self._src_band = None
