import sys

import numpy as np
from netCDF4 import Dataset


def main():
    ndim = 128  # Size of the matrix ndim*ndim
    xdimension = 0.75
    ydimension = 0.75

    file_a = Dataset('data/folder_d/a.nc', 'w', format='NETCDF4')
    file_a.description = 'Test data A'
    file_b = Dataset('data/folder_d/b.nc', 'w', format='NETCDF4')
    file_b.description = 'Test data B'
    file_mask = Dataset('data/folder_d/mask.nc', 'w', format='NETCDF4')
    file_mask.description = 'Test data mask'

    # dimensions
    file_a.createDimension('time', None)
    file_a.createDimension('x', ndim)
    file_a.createDimension('y', ndim)
    file_b.createDimension('time', None)
    file_b.createDimension('x', ndim)
    file_b.createDimension('y', ndim)
    file_mask.createDimension('x', ndim)
    file_mask.createDimension('y', ndim)

    # variables
    time_a = file_a.createVariable('time', 'f8', ('time',))
    x_a = file_a.createVariable('x', 'f4', ('x',))
    y_a = file_a.createVariable('y', 'f4', ('y',))
    field_a = file_a.createVariable('field', 'f8', ('time', 'x', 'y',))
    time_b = file_b.createVariable('time', 'f8', ('time',))
    x_b = file_b.createVariable('x', 'f4', ('x',))
    y_b = file_b.createVariable('y', 'f4', ('y',))
    field_b = file_b.createVariable('field', 'f8', ('time', 'x', 'y',))
    x_c = file_mask.createVariable('x', 'f4', ('x',))
    y_c = file_mask.createVariable('y', 'f4', ('y',))
    field_mask = file_mask.createVariable('field', 'f8', ('x', 'y',))

    # data
    x_range = np.linspace(0, xdimension, ndim)
    y_range = np.linspace(0, ydimension, ndim)
    x_a[:] = x_range
    y_a[:] = y_range
    x_b[:] = x_range
    y_b[:] = y_range
    x_c[:] = x_range
    y_c[:] = y_range
    field_mask[:, :] = np.bool8(np.ones(shape=(len(x_range), len(y_range))))
    for i in range(5):
        time_a[i] = i * 50.0
        time_b[i] = i * 50.0
        values = np.random.uniform(size=(len(x_range), len(y_range)))
        field_a[i, :, :] = values
        field_b[i, :, :] = values + 0.001

    file_a.close()
    file_b.close()
    file_mask.close()


if __name__ == '__main__':
    sys.exit(main())
