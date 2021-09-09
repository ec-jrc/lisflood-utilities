#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from skimage.transform import resize
from tools import *

print('===============================================================================')
        
# Load configuration file
config = pd.read_csv(sys.argv[1],header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    print(string)
    string = string.replace(" ","")    
    try:
        exec(string.replace("=","=r"))     
    except:
        exec(string)         

# Create output folder
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
#   Loop over years and months
############################################################################

for year in np.arange(year_start,year_end+1):
    print('===============================================================================')
    print('Year: '+str(year))
    t0 = time.time()
    
    print('Loading HILDA+ data')
    ind = year-1899
    dset = Dataset(os.path.join(hildaplus_folder,'hildaplus_vGLOB-1.0-f_states.nc'))
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
    
    print('Resampling HILDA+ data to clone map resolution')
    fracwater_init = resize(np.single((hilda_raw==0) | (hilda_raw==77)),mapsize_global,order=1,mode='constant',anti_aliasing=False).astype(np.double)
    fracforest_init = resize(np.single((hilda_raw>=40) & (hilda_raw<=45)),mapsize_global,order=1,mode='constant',anti_aliasing=False).astype(np.double)    
    fracsealed_init = 0.75*resize(np.single(hilda_raw==11),mapsize_global,order=1,mode='constant',anti_aliasing=False).astype(np.double)
    fracother_hilda = 1-fracwater_init-fracforest_init-fracsealed_init    
    del hilda_raw
    
    idx = (np.abs(np.array(hyde_years)-year)).argmin()
    hyde_year = hyde_years[idx]
    print('Loading HYDE data for '+str(hyde_year)+' and resampling to clone map resolution')
    hyde_garea_cr = np.array(pd.read_csv(os.path.join(hyde_folder,'general_files','garea_cr.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double) # total gridcell area in km2
    hyde_garea_cr = hyde_garea_cr[:,:-1] # Rogue last column due to spaces after last value in asc file...
    hyde_cropland = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','cropland'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # total cropland area
    hyde_cropland[np.isnan(hyde_cropland)] = 0
    hyde_cropland = resize(hyde_cropland,mapsize_global,order=1,mode='constant',anti_aliasing=False)    
    hyde_ir_norice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_norice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # irrigated other crops area (no rice) area
    hyde_ir_norice[np.isnan(hyde_ir_norice)] = 0
    hyde_ir_norice = resize(hyde_ir_norice,mapsize_global,order=1,mode='constant',anti_aliasing=False)    
    hyde_ir_rice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_rice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # irrigated rice area (no rice)
    hyde_ir_rice[np.isnan(hyde_ir_rice)] = 0
    hyde_ir_rice = resize(hyde_ir_rice,mapsize_global,order=1,mode='constant',anti_aliasing=False)
    hyde_rf_rice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','rf_rice'+str(hyde_year)+'AD.asc'), 
        header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # rainfed rice area (no rice)
    hyde_rf_rice[np.isnan(hyde_rf_rice)] = 0
    hyde_rf_rice = resize(hyde_rf_rice,mapsize_global,order=1,mode='constant',anti_aliasing=False)    
  
    # hyde_ir_norice exceeds 1 for some grid-cells, which shouldn't be possible
    print('Fixing HYDE fractions')
    totals = hyde_ir_rice+hyde_rf_rice+hyde_ir_norice
    corr_factor = (hyde_cropland+0.000001)/(totals+0.000001)
    corr_factor[corr_factor>1] = 1
    hyde_ir_rice = hyde_ir_rice*corr_factor
    hyde_rf_rice = hyde_rf_rice*corr_factor
    hyde_ir_norice = hyde_ir_norice*corr_factor
    
    print('Making sure HYDE cropland fraction does not exceed HILDA+ other fraction')
    fracrice_hyde = hyde_ir_rice+hyde_rf_rice
    fracirrigation_hyde = hyde_ir_norice.copy()
    corr_factor = (fracother_hilda+0.000001)/(hyde_cropland+0.000001)
    corr_factor[corr_factor>1] = 1
    fracrice_init = fracrice_hyde*corr_factor
    fracirrigation_init = fracirrigation_hyde*corr_factor
    
    print('Make sure sum of fractions (without other) is <1')
    total = fracwater_init+fracforest_init+fracsealed_init+fracrice_init+fracirrigation_init
    total[total<1] = 1
    fracwater_init = (fracwater_init/total).clip(0,1)
    fracforest_init = (fracforest_init/total).clip(0,1)
    fracsealed_init = (fracsealed_init/total).clip(0,1)
    fracrice_init = (fracrice_init/total).clip(0,1)
    fracirrigation_init = (fracirrigation_init/total).clip(0,1)
    
    print('Calculating other fraction using HYDE rice and irrigation (no rice)')
    fracother_init = 1-fracwater_init-fracforest_init-fracsealed_init-fracrice_init-fracirrigation_init

    # The initial water fraction will be replaced with GSWE and the other
    # five fractions will be rescaled accordingly. However, if the the initial
    # water fraction is 1, the other fraction cannot be adjusted, as they will
    # be 0. As a workaround, we reduce the initial water fraction by a tiny
    # amount while increasing the other fractions by a tiny amount using 
    # interpolated (non-zero) values.    
    print('Fixing fully water-covered grid-cells')
    mask = fracwater_init==1
    fracforest_init = fracforest_init+0.000001*fill(fracforest_init,invalid=mask)
    fracsealed_init = fracsealed_init+0.000001*fill(fracsealed_init,invalid=mask)
    fracrice_init = fracrice_init+0.000001*fill(fracrice_init,invalid=mask)
    fracirrigation_init = fracirrigation_init+0.000001*fill(fracirrigation_init,invalid=mask)
    fracother_init = fracother_init+0.000001*fill(fracother_init,invalid=mask)
   
    print('Rescale fractions to sum to 1')
    total = fracwater_init+fracforest_init+fracsealed_init+fracrice_init+fracirrigation_init+fracother_init
    fracwater_init = (fracwater_init/total).clip(0,1)
    fracforest_init = (fracforest_init/total).clip(0,1)
    fracsealed_init = (fracsealed_init/total).clip(0,1)
    fracrice_init = (fracrice_init/total).clip(0,1)
    fracirrigation_init = (fracirrigation_init/total).clip(0,1)
    fracother_init = (fracother_init/total).clip(0,1)

    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    for month in np.arange(1,13):
        print('-------------------------------------------------------------------------------')
        print()
        print('Year: '+str(year)+' Month: '+str(month))
        t0 = time.time()

        idx = (np.abs(np.array(gswe_years)-year)).argmin()
        gswe_year = gswe_years[idx]
        print('Loading GSWE data ('+os.path.join(gswe_folder,str(gswe_year)+'_'+str(month).zfill(2)+'_B.nc')+')')
        dset = Dataset(os.path.join(gswe_folder,str(gswe_year)+'_'+str(month).zfill(2)+'_B.nc'))
        gswe_lats = np.array(dset.variables['lat'][:])
        gswe_lons = np.array(dset.variables['lon'][:])
        gswe_res = gswe_lats[0]-gswe_lats[1]
        add_top = int(np.round((90-gswe_lats[0])/gswe_res))
        add_bottom = int(np.round((90+gswe_lats[-1])/gswe_res))
        gswe_shape = (int(len(gswe_lons)/2),len(gswe_lons))
        gswe_raw = np.zeros(gswe_shape,dtype=np.single)*np.NaN
        varname = list(dset.variables.keys())[-1] # Variable name differs for some years...
        gswe_raw[add_top+1:gswe_shape[0]-add_bottom,:] = dset.variables[varname][:]
        
        print('Resampling GSWE data')
        gswe_resized = resize(gswe_raw,mapsize_global,order=1,mode='constant',anti_aliasing=False).astype(np.double)
        del gswe_raw
        
        print('Inserting GSWE data')
        fracwater = fracwater_init.copy()
        fracwater[np.isnan(gswe_resized)==False] = gswe_resized[np.isnan(gswe_resized)==False]
        
        print('Adjusting other fractions according to GSWE')
        corr_factor = (1-fracwater+0.000001)/(1-fracwater_init+0.000001)
        fracforest = fracforest_init*corr_factor
        fracsealed = fracsealed_init*corr_factor
        fracrice = fracrice_init*corr_factor
        fracirrigation = fracirrigation_init*corr_factor
        fracother = fracother_init*corr_factor        
         
        print('Rescale fractions to sum to 1')
        total = fracwater+fracforest+fracsealed+fracrice+fracirrigation+fracother
        fracwater = (fracwater/total).clip(0,1)
        fracforest = (fracforest/total).clip(0,1)
        fracsealed = (fracsealed/total).clip(0,1)
        fracrice = (fracrice/total).clip(0,1)
        fracirrigation = (fracirrigation/total).clip(0,1)
        fracother = (fracother/total).clip(0,1)
           
        print('Subsetting data to clone map area')
        fracwater = fracwater[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        fracforest = fracforest[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        fracsealed = fracsealed[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        fracrice = fracrice[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        fracirrigation = fracirrigation[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        fracother = fracother[row_upper:row_upper+len(clone_lat),col_left:col_left+len(clone_lon)]
        
        print('Saving data to netCDF')
        vars = ['fracwater','fracforest','fracsealed','fracrice','fracirrigation','fracother']
        for vv in np.arange(len(vars)):
            save_netcdf_3d(
                file = os.path.join(output_folder,vars[vv]+'.nc'),
                varname = vars[vv], 
                index = (year-year_start)*12+month-1,
                data = eval(vars[vv]).astype(np.single),
                varunits = 'fraction',
                timeunits = 'days since 1979-01-02 00:00:00',
                ts = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400,
                least_sig_dig = 3,
                lat = clone_lat,
                lon = clone_lon)
        
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
pdb.set_trace()