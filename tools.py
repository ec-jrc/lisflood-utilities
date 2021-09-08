#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
import pcraster as pcr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import subprocess
import rasterio
from skimage.transform import resize
from scipy import ndimage as nd

def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()
    
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
    
def save_netcdf_3d(file, varname, data, varunits, ts, least_sig_dig):

    if os.path.isfile(file)==False:
        
        mapsize = data.shape
        res = 360/float(mapsize[1])
        lon = np.arange(mapsize[1], dtype=np.single)*res - 180.0 + res/2
        lat = 90-np.arange(mapsize[0], dtype=np.single)*res - res/2
                
        mkdir(os.path.dirname(file))
        
        ncfile = Dataset(file, 'w', format='NETCDF4')
        ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

        ncfile.createDimension('lon', len(lon))
        ncfile.createDimension('lat', len(lat))
        ncfile.createDimension('time', None)

        ncfile.createVariable('lon', 'f4', ('lon',))
        ncfile.variables['lon'][:] = lon
        ncfile.variables['lon'].units = 'degrees_east'
        ncfile.variables['lon'].long_name = 'longitude'

        ncfile.createVariable('lat', 'f4', ('lat',))
        ncfile.variables['lat'][:] = lat
        ncfile.variables['lat'].units = 'degrees_north'
        ncfile.variables['lat'].long_name = 'latitude'

        ncfile.createVariable('time', 'i4', 'time')
        ncfile.variables['time'][:] = (pd.to_datetime(ts)-pd.to_datetime(datetime(1900, 1, 1))).total_seconds()/86400
        ncfile.variables['time'].units = 'days since 1900-1-1 00:00:00'
        ncfile.variables['time'].long_name = 'time'
    
    else:
        ncfile = Dataset(file, 'r+', format='NETCDF4')   
    
    if varname not in ncfile.variables.keys():
        ncfile.createVariable(varname, data.dtype, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999, least_significant_digit=least_sig_dig) #'f4'

    ncfile.variables[varname][0,:,:] = data
    ncfile.variables[varname].units = varunits
    
    ncfile.close()
    
def fill(data, invalid=None):
    # Nearest neighbor interpolation gap fill by Juh_
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]  