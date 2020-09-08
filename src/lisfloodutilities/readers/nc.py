import xarray as xr


class NetCDFMap:
    """
    A class representing a netCDF file.
    """
    def __init__(self, netcdf_filename):
        self.ds = xr.open_dataset(netcdf_filename)
