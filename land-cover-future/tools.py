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
from scipy import ndimage as nd

def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()
            
def imresize_mean(oldarray,newshape):
    '''Resample an array using simple averaging'''
    
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

def fill(data, invalid=None):
    '''Nearest neighbor interpolation gap fill by Juh_'''
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]    
    
def load_config(filepath):
    '''Load configuration file into dict'''
    
    df = pd.read_csv(filepath,header=None,index_col=False)
    config = {}
    for ii in np.arange(len(df)): 
        string = df.iloc[ii,0].replace(" ","")
        varname = string.rpartition('=')[0]
        varcontents = string.rpartition('=')[2]
        try:
            varcontents = float(varcontents)
        except:
            pass
        config[varname] = varcontents
    return config
    
def load_data_interp(year_target,file_pattern):
    '''Linear interplation between maps for different years.'''
    
    files = glob.glob(file_pattern)
    files = sorted(files)
    files_years = np.array([int(os.path.splitext(os.path.basename(file))[0][-4:]) for file in files])

    # Backward nearest neighbor extrapolation
    if year_target<=files_years[0]:
        return np.load(files[0])['data']        
        
    # Forward nearest neighbor extrapolation
    elif year_target>=files_years[-1]:
        return np.load(files[-1])['data']
    
    # Linear interpolation
    else:
        ind_sort = np.argsort(np.abs(np.array(files_years)-year_target))[:2]
        ind_sort = np.sort(ind_sort)
        year_1 = files_years[ind_sort[0]]
        year_2 = files_years[ind_sort[1]]
        period = year_2-year_1
        contribution_map_2 = (year_target-year_1)/period
        contribution_map_1 = 1-contribution_map_2
        map_1 = np.load(files[ind_sort[0]])['data']
        map_2 = np.load(files[ind_sort[1]])['data']
        return map_1*contribution_map_1+map_2*contribution_map_2
    