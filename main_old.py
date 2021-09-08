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
    
def fill(data, invalid=None):
    # Nearest neighbor interpolation gap fill by Juh_
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]    
        
# Load configuration file
config = pd.read_csv('config.cfg',header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")    
    try:
        exec(string.replace("=","=r"))     
    except:
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

# Determine map sizes
mapsize_global = (np.round(180/clone_res).astype(int),np.round(360/clone_res).astype(int))
mapsize_clone = clone_np.shape
row_upper,col_left = latlon2rowcol(clone_lat[0],clone_lon[0],clone_res,90,-180)

# List of years with HYDE data
hyde_files = glob.glob(os.path.join(hyde_folder,'baseline','zip','cropland*'))
hyde_years = [int(os.path.basename(hyde_file)[8:12]) for hyde_file in hyde_files]

# List of years with GSWE data
gswe_files = glob.glob(os.path.join(gswe_folder,'*'))
gswe_years = [int(os.path.basename(gswe_file)[:4]) for gswe_file in gswe_files]
gswe_years = np.unique(gswe_years)


############################################################################
#   Load data
############################################################################

dset = Dataset(os.path.join(hildaplus_folder,'hildaplus_vGLOB-1.0-f_states.nc'))
latitude = np.array(dset.variables['latitude'][:])
longitude = np.array(dset.variables['longitude'][:])
for year in np.arange(year_start,year_end+1):
    print('===============================================================================')
    print('Year: '+str(year))
    t0 = time.time()
    
    print('Loading HILDA+ data')
    ind = year-1899
    hilda_raw = np.array(dset.variables['LULC_states'][ind,:,:])
    
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
    
    print('Resampling HILDA+ data')
    fracwater_hilda = resize(np.single((hilda_raw==0) | (hilda_raw==77)),mapsize_global,order=1,mode='constant',anti_aliasing=False)
    fracforest_hilda = resize(np.single((hilda_raw>=40) & (hilda_raw<=45)),mapsize_global,order=1,mode='constant',anti_aliasing=False)    
    fracsealed_hilda = resize(np.single(hilda_raw==11),mapsize_global,order=1,mode='constant',anti_aliasing=False)
    fraccrop_hilda = resize(np.single(hilda_raw==22),mapsize_global,order=1,mode='constant',anti_aliasing=False)
    del hilda_raw
    
    # The HILDA+ cropland fraction will be replaced with HYDE and the HILDA+ 
    # water fraction will be replaced with GSWE and the other fractions
    # (forest and sealed) will be adjusted accordingly. However, if the 
    # HILDA+ cropland or water fractions are 1, the other fraction cannot be 
    # adjusted, as they will be 0. As a workaround, we reduce the HILDA+ 
    # cropland and water fractions by a tiny amount while increasing the 
    # other fractions by a tiny amount using interpolated non-zero values.
    
    print('Fixing fully water-covered HILDA+ grid-cells')
    mask = fracwater_hilda==1
    fracwater_hilda[mask] = 1-0.000001
    fracforest_hilda = fracforest_hilda+0.000001*fill(fracforest_hilda,invalid=mask)
    fracsealed_hilda = fracsealed_hilda+0.000001*fill(fracsealed_hilda,invalid=mask)
    fraccrop_hilda = fraccrop_hilda+0.000001*fill(fraccrop_hilda,invalid=mask)
    
    print('Fixing fully crop-covered HILDA+ grid-cells')
    mask = fraccrop_hilda==1
    fraccrop_hilda[mask] = 1-0.000001
    fracforest_hilda = fracforest_hilda+0.000001*fill(fracforest_hilda,invalid=mask)
    fracsealed_hilda = fracsealed_hilda+0.000001*fill(fracsealed_hilda,invalid=mask)
    fracwater_hilda = fracwater_hilda+0.000001*fill(fracwater_hilda,invalid=mask)

    print('Loading HYDE data')
    idx = (np.abs(np.array(hyde_years)-year)).argmin()
    hyde_year = hyde_years[idx]
    hyde_garea_cr = np.array(pd.read_csv(os.path.join(hyde_folder,'general_files','garea_cr.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.single) # total gridcell area in km2
    hyde_garea_cr = hyde_garea_cr[:,:-1] # Rogue last column due to spaces after last value in asc file
    hyde_cropland = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','cropland'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.single)/hyde_garea_cr # total cropland area
    hyde_cropland[np.isnan(hyde_cropland)] = 0
    hyde_cropland = resize(hyde_cropland,mapsize_global,order=0,mode='constant',anti_aliasing=False)    
    hyde_ir_norice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_norice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.single)/hyde_garea_cr # irrigated other crops area (no rice) area
    hyde_ir_norice[np.isnan(hyde_ir_norice)] = 0
    hyde_ir_norice = resize(hyde_ir_norice,mapsize_global,order=0,mode='constant',anti_aliasing=False)    
    hyde_ir_rice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_rice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.single)/hyde_garea_cr # irrigated rice area (no rice)
    hyde_ir_rice[np.isnan(hyde_ir_rice)] = 0
    hyde_ir_rice = resize(hyde_ir_rice,mapsize_global,order=0,mode='constant',anti_aliasing=False)
    hyde_rf_rice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','rf_rice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.single)/hyde_garea_cr # rainfed rice area (no rice)
    hyde_rf_rice[np.isnan(hyde_rf_rice)] = 0
    hyde_rf_rice = resize(hyde_rf_rice,mapsize_global,order=0,mode='constant',anti_aliasing=False)    
 
    print('Adjusting HYDE fractions based on HILDA+ sealed')
    total = fracsealed_hilda+hyde_cropland
    corr = 1-total
    corr[corr>1
    hyde_cropland = hyde_cropland-corr
    hyde_ir_norice = hyde_ir_norice*corr_factor
    hyde_ir_rice = hyde_ir_rice*corr_factor
    hyde_rf_rice = hyde_rf_rice*corr_factor

    print('Adjusting HILDA+ fractions based on HYDE cropland')
    fraccrop_init = hyde_cropland.copy()
    fracirrigation_init = hyde_ir_norice.copy()
    fracrice_init = hyde_ir_rice+hyde_rf_rice
    fracsealed_init = fracsealed_hilda.copy()
    corr_factor = (1-fraccrop_init-fracsealed_init)/(1-fraccrop_hilda-fracsealed_init)
    fracforest_init = fracforest_hilda*corr_factor
    fracwater_init = fracwater_hilda*corr_factor
    fracother_init = 1-fracirrigation_init-fracrice_init-fracsealed_init-fracforest_init-fracwater_init
    
    plt.figure(0)
    plt.imshow(fracother_init)
    plt.colorbar()
    plt.title('fracother')
    
    plt.show(block=False)

    pdb.set_trace()

    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    for month in np.arange(1,13):
        print('-------------------------------------------------------------------------------')
        print('Month: '+str(month))
        t0 = time.time()

        print('Loading GSWE data')
        idx = (np.abs(np.array(gswe_years)-year)).argmin()
        gswe_year = gswe_years[idx]
        dset = Dataset(os.path.join(gswe_folder,str(gswe_year)+'_'+str(month).zfill(2)+'_B.nc'))
        gswe_lats = np.array(dset.variables['lat'][:])
        gswe_lons = np.array(dset.variables['lon'][:])
        gswe_res = gswe_lats[0]-gswe_lats[1]
        add_top = int(np.round((90-gswe_lats[0])/gswe_res))
        add_bottom = int(np.round((90+gswe_lats[-1])/gswe_res))
        gswe_shape = (int(len(gswe_lons)/2),len(gswe_lons))
        gswe_raw = np.zeros(gswe_shape,dtype=np.single)*np.NaN
        gswe_raw[add_top+1:gswe_shape[0]-add_bottom,:] = dset.variables['GSWD_'+str(gswe_year)+' B'][:]
        
        print('Resampling GSWE data')
        gswe_resized = resize(gswe_raw,mapsize_global,order=1,mode='constant',anti_aliasing=False)
        del gswe_raw
            
        print('Inserting GSWE data')
        fracwater = fracwater_init.copy()
        fracwater[np.isnan(gswe_resized)==False] = gswe_resized[np.isnan(gswe_resized)==False]
        
        print('Adjusting other fractions (forest, sealed, crop) according to GSWE')
        corr_factor = (1-fracwater)/(1-fracwater_init)
        fracforest = fracforest_init*corr_factor
        fracsealed = fracsealed_init*corr_factor
        fraccrop = fraccrop_init*corr_factor
        
        pdb.set_trace()
        
        plt.figure(0)
        plt.imshow(hyde_cropland)
                
        plt.figure(1)
        plt.imshow(fraccrop)
        plt.show()
        
    
        
        
        print('Subset maps to area')
        fracforest = fracforest[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        
        print("Time elapsed is "+str(time.time()-t0)+" sec")

pdb.set_trace()

plt.show(block=False)

pdb.set_trace()