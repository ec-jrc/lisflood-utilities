#!/usr/bin/env python  
# -*- coding: utf-8 -*-
"""
Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""
__author__ = "Hylke E. Beck"
__date__ = "January 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from skimage.transform import resize
from skimage.transform import downscale_local_mean
from datetime import datetime, timedelta
from scipy import ndimage as nd
import rasterio
        
def load_country_code_map(filepath,mapsize):
    
    src = rasterio.open(filepath)
    country_code_map = src.read(1).astype(np.single)
    country_code_map_refmatrix = src.get_transform()
    src.close()

    # Check if country border raster covers entire globe
    assert (country_code_map_refmatrix[0]==-180) & (country_code_map_refmatrix[3]==90)

    country_code_map = resize(country_code_map,mapsize,order=0,mode='edge',anti_aliasing=False)
    # apply the fill algorithm to fill the gaps of the Disputed Territories (code>1000), so that the fill algorithm will set each pixel value to their nearest neighbour
    country_code_map = fill(country_code_map, country_code_map>1000) # CR: check here: the fill algorithm includes also zero as a neighbour country code, so that isles and some land area will be replaced by ocean
    country_code_map[country_code_map==0] = np.NaN
    
    return country_code_map
    
def load_us_state_code_map(filepath,mapsize):
    
    src = rasterio.open(filepath)
    state_code_map = src.read(1).astype(np.single)
    state_code_map = resize(state_code_map,mapsize,order=0,mode='edge',anti_aliasing=False)
    state_code_map_refmatrix = src.get_transform()
    src.close()

    # Check if state border raster covers entire globe
    assert (state_code_map_refmatrix[0]==-180) & (state_code_map_refmatrix[3]==90)

    return state_code_map
    
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