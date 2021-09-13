#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from skimage.transform import resize
from skimage.transform import downscale_local_mean
from datetime import datetime, timedelta

def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()

def get_row_compressor(old_dimension, new_dimension):
    dim_compressor = np.zeros((new_dimension, old_dimension))
    bin_size = float(old_dimension) / new_dimension
    next_bin_break = bin_size
    which_row = 0
    which_column = 0
    while which_row < dim_compressor.shape[0] and which_column < dim_compressor.shape[1]:
        if round(next_bin_break - which_column, 10) >= 1:
            dim_compressor[which_row, which_column] = 1
            which_column += 1
        elif next_bin_break == which_column:

            which_row += 1
            next_bin_break += bin_size
        else:
            partial_credit = next_bin_break - which_column
            dim_compressor[which_row, which_column] = partial_credit
            which_row += 1
            dim_compressor[which_row, which_column] = 1 - partial_credit
            which_column += 1
            next_bin_break += bin_size
    dim_compressor /= bin_size
    return dim_compressor

def get_column_compressor(old_dimension, new_dimension):
    return get_row_compressor(old_dimension, new_dimension).transpose()
    
def compress_and_average(array, new_shape):
    # Note: new shape should be smaller in both dimensions than old shape
    return np.mat(get_row_compressor(array.shape[0], new_shape[0])) * \
           np.mat(array) * \
           np.mat(get_column_compressor(array.shape[1], new_shape[1]))
           
def imresize_majority(oldarray,newshape):
    # Resize array using majority filter
    
    cls = np.unique(oldarray)
    
    if len(cls)>100:
        raise ValueError('More than 100 classes, are you sure this is a thematic map?')    
    
    threshold = 1/len(cls)
    newarray = np.zeros(newshape,dtype=type(oldarray))
    for cl in cls:        
        temp = resize((oldarray==cl).astype(np.single),newshape,order=1,mode='constant',anti_aliasing=False)
        newarray[temp>=threshold] = cl

    return newarray
    
def imresize_mean(oldarray,newshape):
    # Resize array using averaging
    
    oldshape = oldarray.shape
    
    factor = oldshape[0]/newshape[0]
    
    if factor==int(factor):
        factor = int(factor)
        newarray = downscale_local_mean(oldarray,(factor,factor))
        
    else:
        factor = 1
        while newshape[0]*factor<oldshape[0]:
            factor = factor+1        
        intarray = resize(oldarray,(int(newshape[0]*factor),int(newshape[1]*factor)),order=0,mode='constant',anti_aliasing=False)
        newarray = downscale_local_mean(intarray,(factor,factor))

    return newarray
    
def save_netcdf_3d(file, varname, index, data, varunits, timeunits, ts, least_sig_dig, lat, lon):

    if os.path.isfile(file)==False:
    
        if os.path.isdir(os.path.dirname(file))==False:
            os.mkdir(os.path.dirname(file))
        
        ncfile = Dataset(file, 'w', format='NETCDF4')
        ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

        ncfile.createDimension('lon', len(lon))
        ncfile.createDimension('lat', len(lat))
        ncfile.createDimension('time', None)

        ncfile.createVariable('lon', 'f4', ('lon',))
        ncfile.variables['lon'][:] = lon
        ncfile.variables['lon'].units = 'degrees_east'
        ncfile.variables['lon'].long_name = 'longitude'

        ncfile.createVariable('lat', 'f8', ('lat',))
        ncfile.variables['lat'][:] = lat
        ncfile.variables['lat'].units = 'degrees_north'
        ncfile.variables['lat'].long_name = 'latitude'

        ncfile.createVariable('time', 'f8', 'time')
        ncfile.variables['time'].units = timeunits
        ncfile.variables['time'].long_name = 'time'
    
    else:
        ncfile = Dataset(file, 'r+', format='NETCDF4')   
    
    if varname not in ncfile.variables.keys():
        ncfile.createVariable(varname, data.dtype, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999, least_significant_digit=least_sig_dig)
    
    ncfile.variables['time'][index] = ts
    ncfile.variables[varname][index,:,:] = data
    ncfile.variables[varname].units = varunits
    
    ncfile.close()

def save_netcdf(file, varname, data, least_sig_dig, lat, lon):

    if os.path.isfile(file)==True: 
        os.remove(file)

    ncfile = Dataset(file, 'w', format='NETCDF4')

    ncfile.createDimension('lon', len(lon))
    ncfile.createDimension('lat', len(lat))

    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'

    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'][:] = lat
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'].long_name = 'latitude'
    
    ncfile.createVariable(varname, data.dtype, ('lat', 'lon'), zlib=True, chunksizes=(32,32,), fill_value=-9999, least_significant_digit=least_sig_dig)

    ncfile.variables[varname][:,:] = data
    
    ncfile.close()
    
def fill(data, invalid=None):
    # Nearest neighbor interpolation gap fill by Juh_
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]  