import xarray as xr

coordinates_names = ('lat', 'lon', 'x', 'y', 'longitude', 'latitude')


class NetCDFMap:
    """
    A class representing a netCDF map.
    Current version does not support mapstacks.
    e.g. To convert ldd.nc
    """
    def __init__(self, netcdf_filename):
        self.ds = xr.open_dataset(netcdf_filename)

    @property
    def data(self):
        for varname, variable in self.ds.variables.items():
            if len(variable.dims) < 2:
                continue
            yield varname, variable.values

    @property
    def coordinates(self):
        for varname, variable in self.ds.variables.items():
            if varname.lower() not in coordinates_names:
                continue
            yield varname, variable.values

    def close(self):
        self.ds.close()
