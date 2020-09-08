import xarray as xr


class NetCDFMap:
    """
    A class representing a netCDF map.
    Current version does not support mapstacks.
    """
    def __init__(self, netcdf_filename):
        self.ds = xr.open_dataset(netcdf_filename)

    @property
    def data(self):
        for varname, variable in self.ds.variables:
            if len(variable.dims) < 2:
                continue
            yield varname, variable.values

    def close(self):
        self.ds.close()
