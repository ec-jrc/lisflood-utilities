import xarray as xr

x_coordinates_names = ('lon', 'x', 'longitude')
y_coordinates_names = ('lat', 'y', 'latitude')
coordinates_names = x_coordinates_names + y_coordinates_names


class NetCDFMap:
    """
    A class representing a netCDF map.
    Current version does not support mapstacks.
    e.g. To convert ldd.nc
    """
    def __init__(self, netcdf_filename):
        self.ds = xr.open_dataset(netcdf_filename, decode_cf=False)

    @property
    def mv(self):
        for variable in self.ds.variables.values():
            if len(variable.dims) < 2:
                continue
            return variable.attrs['_FillValue'] if '_FillValue' in variable.attrs else variable.attrs['missing_value']

    @property
    def data(self):
        for varname, variable in self.ds.variables.items():
            if len(variable.dims) < 2:
                continue
            yield varname, variable

    @property
    def coordinates(self):
        res = {}
        for varname, variable in self.ds.variables.items():
            if varname.lower() not in coordinates_names:
                continue
            key = 'x' if varname.lower() in x_coordinates_names else 'y'
            res[key] = variable.values
        return res

    def close(self):
        self.ds.close()
