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
from skimage.transform import resize

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

def save_netcdf(file, varname, data, lat, lon):

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
    
    ncfile.createVariable(varname, data.dtype, ('lat', 'lon'), zlib=True, chunksizes=(32,32,), fill_value=-9999)

    ncfile.variables[varname][:,:] = data
    
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

if os.path.isdir(temp_folder)==False:
    os.mkdir(temp_folder)

if os.path.isdir(output_folder)==False:
    os.mkdir(output_folder)
   
# Load clone map
dset = Dataset(clonemap_path)
clone_lat = dset.variables['lat'][:]
clone_lon = dset.variables['lon'][:]
clone_res = clone_lon[1]-clone_lon[0]
varname = list(dset.variables.keys())[-1]
clone_np = np.array(dset.variables[varname][:])
#pcr.setclone(clone_np.shape[0],clone_np.shape[1],clone_res,clone_lon[0]-clone_res/2,clone_lat[0]-clone_res/2)

mapsize_global = (np.round(180/clone_res).astype(int),np.round(360/clone_res).astype(int))
mapsize_area = clone_np.shape


############################################################################
#   Load data
############################################################################

dset = Dataset(os.path.join(hildaplus_folder,'hildaplus_vGLOB-1.0-f_states.nc'))
latitude = np.array(dset.variables['latitude'][:])
longitude = np.array(dset.variables['longitude'][:])
for year in np.arange(year_start,year_end+1):
    print(year)
    t0 = time.time()
    
    # Read raw data (global)
    ind = year-1899
    raw = np.array(dset.variables['LULC_states'][ind,:,:])
    
    # Resample raw data (nearest neighbour)
    data = np.round(resize(raw,mapsize_global,order=0,mode='constant',anti_aliasing=False)*255).astype(int)
    
    # HILDA+ legend
    # 11 Urban 
    # 22 Cropland 
    # 33 Pasture 
    # 40: Forest (Unknown/Other) 
    # 41: Forest (Evergreen, needle leaf) 
    # 42: Forest ( Evergreen, broad leaf) 
    # 43: Forest (Deciduous, needle leaf) 
    # 44: Forest (Deciduous, broad leaf) 
    # 45: Forest (Mixed) 
    # 55 Grass/shrubland 
    # 66 Other land 
    # 77 Water 
    fracwater = data==77
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    plt.imshow(data)
    plt.colorbar()
    plt.show()
    pdb.set_trace()

pdb.set_trace()


############################################################################
#   Make global upstream map at specified resolution based on MERIT Hydro 
#   upstream data
############################################################################

if os.path.isfile(os.path.join(output_folder,'upstream_area_global.npy'))==False:

    global_shape = (int(180/res),int(360/res))
    upstream_area_global = np.zeros(global_shape)*np.NaN
    
    for subdir, dirs, files in os.walk(merit_folder):
        for file in files:        
            
            print('--------------------------------------------------------------------------------')
            print('Processing '+file)
            t1 = time.time()
            
            # Resize using maximum filter
            oldarray = rasterio.open(os.path.join(subdir, file)).read(1)
            oldarray[oldarray<0] = 9999999 # Necessary to ensure that all rivers flow into the ocean
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

    with open(os.path.join(output_folder,'upstream_area_global.npy'), 'wb') as f:
        np.save(f, upstream_area_global)

# Load resampled global upstream area map from disk
with open(os.path.join(output_folder,'upstream_area_global.npy'), 'rb') as f:
    upstream_area_global = np.load(f)

    
############################################################################
#   Create ldd based on global upstream area map
############################################################################

# Load clone map
dset = Dataset(clonemap_path)
clone_lat = dset.variables['lat'][:]
clone_lon = dset.variables['lon'][:]
clone_res = clone_lon[1]-clone_lon[0]
varname = list(dset.variables.keys())[-1]
clone_np = np.array(dset.variables[varname][:])
pcr.setclone(clone_np.shape[0],clone_np.shape[1],clone_res,clone_lon[0]-clone_res/2,clone_lat[0]-clone_res/2)

# Create grid-cell area map
xi, yi = np.meshgrid(clone_lon, clone_lat)
area_np = (40075*res/360)**2*np.cos(np.deg2rad(yi))
area_pcr = pcr.numpy2pcr(pcr.Scalar,area_np,mv=-9999)

# Subset global upstream map to clone map region
if np.round(clone_res*1000000)!=np.round(res*1000000):
    raise ValueError('Clone resolution of '+str(clone_res)+' does not match target resolution of '+str(res))
row_upper,col_left = latlon2rowcol(clone_lat[0],clone_lon[0],res,90,-180)
upstream_area = upstream_area_global[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]

# Produce ldd by inverting upstream map
fake_elev_np = 9999999-upstream_area
fake_elev_np[np.isnan(fake_elev_np)] = 0
fake_elev_pcr = pcr.numpy2pcr(pcr.Scalar,fake_elev_np,mv=9999999)
ldd_pcr = pcr.lddcreate(fake_elev_pcr,0,0,0,0)
ldd_np = pcr.pcr2numpy(ldd_pcr,mv=-9999)
upstreamarea_pcr = pcr.accuflux(ldd_pcr,area_pcr)
upstreamarea_np = pcr.pcr2numpy(upstreamarea_pcr,mv=9999999)

# Save results
save_netcdf(os.path.join(output_folder,'ldd.nc'), 'ldd', ldd_np, clone_lat, clone_lon)
save_netcdf(os.path.join(output_folder,'ups.nc'), 'ups', upstreamarea_np, clone_lat, clone_lon)


############################################################################
#   Plot some maps
############################################################################

plt.figure(1)
plt.imshow(np.log10(upstreamarea_np),vmin=0,vmax=6)
plt.title('New upstream area (upstreamarea_np)')

plt.figure(2)
plt.imshow(np.log10(upstream_area_global),vmin=0,vmax=6)
plt.title('MERIT Hydro upstream area (upstream_area_global)')

plt.show(block=False)

pdb.set_trace()