#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "August 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
import pcraster as pcr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import subprocess
import rasterio

def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()
    
def imresize_max(oldarray,newshape):
    # Resize array using max filter. Inefficient, but works.
    
    # Determine resize factor
    oldshape = oldarray.shape
    factor = oldshape[0]/newshape[0]    
    if factor!=np.round(factor): 
        raise ValueError('Resize factor of '+str(factor)+' not integer, needs to be integer')    
    factor = factor.astype(int)
    
    # Loop over new grid and compute max
    newarray = np.zeros(newshape,dtype=np.single)*np.NaN
    for ii in np.arange(newshape[0]):
        for jj in np.arange(newshape[1]):
            newarray[ii,jj] = np.nanmax(oldarray[ii*factor:(ii+1)*factor,jj*factor:(jj+1)*factor])
    
    return newarray

def save_netcdf_2d(file, varname, data, lat, lon, varunits, least_sig_dig):

    if os.path.isfile(file)==True: os.remove(file)
    
    mkdir(os.path.dirname(file))
    
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
    ncfile.variables[varname].units = varunits
    
    ncfile.close()
    
# Load configuration file
config = pd.read_csv('config.cfg',header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")    
    try:
        exec(string) 
    except:
        string = string.replace("=","=r")
        exec(string) 

if os.path.isdir(merit_folder)==False:
    os.mkdir(merit_folder)
    

############################################################################
#   Download and untar all MERIT Hydro upstream data
############################################################################

url_pre = 'http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/distribute/v1.0/'
username = 'hydrography'
password = 'rivernetwork'

lats = ['n60','n30','n00','s30','s60']
lons = ['w180','w150','w120','w090','w060','w030','e000','e030','e060','e090','e120','e150']
for lat in lats:
    for lon in lons:
        if os.path.isdir(os.path.join(merit_folder,'upa_'+lat+lon)): continue
        filename = 'upa_'+lat+lon+'.tar'
        command = 'wget '+url_pre+filename+' --no-clobber --user='+username+' --password='+password+' --directory-prefix='+merit_folder
        subprocess.call(command,shell=True)
        command = 'tar -xvf '+os.path.join(merit_folder,filename)+' -C '+merit_folder
        subprocess.call(command,shell=True)
        if os.path.isfile(os.path.join(merit_folder,filename)):
            os.remove(os.path.join(merit_folder,filename))
        

############################################################################
#   Make global upstream map at specified resolution 
#   based on MERIT Hydro data
############################################################################

global_shape = (int(180/res),int(360/res))
upstream_area_global = np.zeros(global_shape,dtype=np.single)*np.NaN
'''
for subdir, dirs, files in os.walk(merit_folder):
    for file in files:        
        
        print('--------------------------------------------------------------------------------')
        print('Processing '+file)
        t1 = time.time()
        
        # Resize using maximum filter
        oldarray = rasterio.open(os.path.join(subdir, file)).read(1)
        factor = res/(5/6000)
        factor = np.round(factor*1000000000000)/1000000000000
        if factor!=np.round(factor): 
            raise ValueError('Resize factor of '+str(factor)+' not integer, needs to be integer')
        factor = factor.astype(int)
        newshape = (oldarray.shape[0]//factor,oldarray.shape[1]//factor)
        newarray = imresize_max(oldarray,newshape)
        
        # Insert into global map
        tile_lat_bottom = float(file[:3].replace("n","").replace("s","-"))
        tile_lon_left = float(file[3:7].replace("e","").replace("w","-"))
        tile_row_bottom, tile_col_left = latlon2rowcol(tile_lat_bottom+res/2,tile_lon_left+res/2,res,90,-180)
        upstream_area_global[tile_row_bottom-newshape[0]+1:tile_row_bottom+1,tile_col_left:tile_col_left+newshape[1]] = newarray
        
        print('Time elapsed is ' + str(time.time() - t1) + ' sec')
'''

############################################################################
#   Create ldd based on global upstream area map
############################################################################

# Use upstream area as elevation to derive ldd
dset = Dataset(clonemap)
clone_lat = dset.variables['lat'][:]
clone_lon = dset.variables['lon'][:]
clone_res = clone_lon[1]-clone_lon[0]
varname = list(dset.variables.keys())[-1]
clone_np = np.array(dset.variables[varname][:])
pcr.setclone(clone_np.shape[0],clone_np.shape[1],clone_res,clone_lon[0]-clone_res/2,clone_lat[0]-clone_res/2)

# Subset global ldd to region defined by clonemap
if np.round(clone_res*1000000)!=np.round(res*1000000):
    raise ValueError('Clone resolution of '+str(clone_res)+' does not match target resolution of '+str(res))
row_upper,col_left = latlon2rowcol(clone_lat[0],clone_lon[0],res,90,-180)
upstream_area = upstream_area_global[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]

upstream_area_pcr = pcr.numpy2pcr(pcr.Scalar,upstream_area,mv=-999999999)
pdb.set_trace()
ldd_global = pcr.lddcreate(-upstream_area_pcr,9999999,9999999,9999999,9999999)



# Produce netcdf
if os.path.isdir()==False:
    mkdir(output_folder)
    
ldd_subset

pdb.set_trace()

plt.figure(0)
plt.imshow(upstream_area_global)
plt.colorbar()
plt.show()
        

pdb.set_trace()