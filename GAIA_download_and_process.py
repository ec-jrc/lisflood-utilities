#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import geopandas as gpd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import subprocess
import rasterio
from tools import *
    
# Load configuration file
config = pd.read_csv(sys.argv[1],header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")    
    try:
        exec(string.replace("=","=r")) 
    except:
        exec(string) 

if os.path.isdir(os.path.join(gaia_folder,'geotiffs'))==False:
    os.makedirs(os.path.join(gaia_folder,'geotiffs'))
    
if os.path.isdir(os.path.join(gaia_folder,'netCDF_1km'))==False:
    os.makedirs(os.path.join(gaia_folder,'netCDF_1km'))

if use_proxy==True:
    proxy_credentials = pd.read_csv('proxy_credentials.cfg',header=None,index_col=False)
    
# Load GAIA shapefile
GAIA_shape = gpd.read_file(os.path.join(gaia_folder,'GAIA_shape','GAIA_nameID_1deg.shp'))


############################################################################
#   Download and untar all MERIT Hydro upstream data
############################################################################

url_pre = 'http://data.ess.tsinghua.edu.cn/data/GAIA/'

for ii in np.arange(GAIA_shape.shape[0]):    
    string = GAIA_shape.loc[ii,'fName_ID']
    
    # Fix string padding...
    underscore = string.find('_')
    lonstr = str(string[underscore+1:]).zfill(4)
    if lonstr[0]=='0': lonstr = lonstr[1:]
    latstr = str(string[:underscore]).zfill(3)
    if latstr[0]=='0': latstr = latstr[1:]
    string = latstr+'_'+lonstr
    
    filename = 'GAIA_1985_2018_'+string+'.tif'
    if os.path.isfile(os.path.join(gaia_folder,'geotiffs',filename)): continue
    if use_proxy==True:
        proxy_string = '-e use_proxy=yes -e http_proxy=http://'+proxy_credentials.iloc[0,0]+':'+proxy_credentials.iloc[1,0]+'@'+proxy_credentials.iloc[2,0]+' '
    else:
        proxy_string = ''
    print('===============================================================================')
    command = 'wget '+proxy_string+url_pre+filename+' --no-clobber --directory-prefix='+os.path.join(gaia_folder,'geotiffs')
    subprocess.call(command,shell=True)


############################################################################
#   Loop over years and make impervious maps
############################################################################

years = np.arange(1985,2019)
values = np.arange(34,0,-1)
if len(years)!=len(values):
    raise ValueError('Something wrong')
    
geotiff_files = glob.glob(os.path.join(gaia_folder,'geotiffs','*.tif'))
for jj in np.arange(len(years)):
    t0 = time.time()
    year = years[jj]
    value = values[jj]
    print('Processing '+str(year))
    
    output_res = 0.01
    global_impervious = np.zeros((int(180/output_res),int(360/output_res)),dtype=np.single)
    for ii in np.arange(len(geotiff_files)):
        
        # Load and resample tile
        raw = rasterio.open(geotiff_files[ii]).read(1)
        impervious = raw>=value
        impervious = imresize_mean(np.single(impervious),(int(1/output_res),int(1/output_res)))
           
        # Find top row and left column of tile
        filename = os.path.basename(geotiff_files[ii])
        try:
            tile_lat_top = float(filename[15:18])
        except:
            tile_lat_top = float(filename[15:17])
        try:
            tile_lon_left = float(filename[-8:-4])
        except:
            tile_lon_left = float(filename[-7:-4])
        tile_row_top, tile_col_left = latlon2rowcol(tile_lat_top-output_res/2,tile_lon_left+output_res/2,output_res,90,-180)
        
        # Insert tile into global map
        global_impervious[tile_row_top:tile_row_top+impervious.shape[0],tile_col_left:tile_col_left+impervious.shape[1]] = impervious
        
    # Save global map to netCDF
    save_netcdf(
        file=os.path.join(gaia_folder,'netCDF_1km',str(year)+'.nc'), 
        varname='impervious_fraction', 
        data=global_impervious, 
        least_sig_dig=3, 
        lat=np.arange(90-output_res/2,-90-output_res/2,-output_res), 
        lon=np.arange(-180+output_res/2,180+output_res/2,output_res)
        )
        
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
pdb.set_trace()