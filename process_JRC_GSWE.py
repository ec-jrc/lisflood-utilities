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
from tools import *

# Load configuration file
config = load_config(sys.argv[1])

def main():
    
    # Load template map
    dset = Dataset(config['templatemap_path'])
    template_lat = dset.variables['lat'][:]
    template_lon = dset.variables['lon'][:]
    template_res = template_lon[1]-template_lon[0]
    varname = list(dset.variables.keys())[-1]
    template_np = np.array(dset.variables[varname][:])

    # Determine map sizes
    mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
    mapsize_template = template_np.shape
    row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)
    
    # Create output folder
    if os.path.isdir(os.path.join(config['output_folder'],'JRC_GSWE'))==False:
        os.makedirs(os.path.join(config['output_folder'],'JRC_GSWE'))
        
    # Climatology years
    years = np.arange(2010,2020)    
    
    # Loop over months and years
    months = np.arange(1,13)
    for month in months:
        print('-------------------------------------------------------------------------------')
        print('Computing fracwater climatology for month '+str(month))
        t0 = time.time()
        fracwater_all = np.zeros((mapsize_global[0],mapsize_global[1],len(years)),dtype=np.single)
        for jj in np.arange(len(years)):
            year = years[jj]
        
            # Load HILDA+ data            
            ind = year-1899
            dset = Dataset(os.path.join(config['hildaplus_folder'],'hildaplus_vGLOB-1.0-f_states.nc'))
            hilda_raw = np.array(dset.variables['LULC_states'][ind,:,:])
            dset.close()
            fracwater_hilda = imresize_mean(np.single((hilda_raw==0) | (hilda_raw==77)),mapsize_global).astype(np.single)
            del hilda_raw
        
            # Load GSWE data
            month_gswe = month.copy()
            if (year==2018) & (month==10): month_gswe = 9 # October 2018 GSWE data are erroneous
            dset = Dataset(os.path.join(config['gswe_folder'],str(year)+'_'+str(month_gswe).zfill(2)+'_B.nc'))
            gswe_lats = np.array(dset.variables['lat'][:])
            gswe_lons = np.array(dset.variables['lon'][:])
            gswe_res = gswe_lats[0]-gswe_lats[1]
            add_top = int(np.round((90-gswe_lats[0])/gswe_res))
            add_bottom = int(np.round((90+gswe_lats[-1])/gswe_res))
            gswe_shape = (int(len(gswe_lons)/2),len(gswe_lons))
            gswe_raw = np.zeros(gswe_shape,dtype=np.single)*np.NaN
            varname = list(dset.variables.keys())[-1] # Variable name differs for some years...
            gswe_raw[add_top+1:gswe_shape[0]-add_bottom,:] = dset.variables[varname][:].clip(0,1)
            fracwater = imresize_mean(gswe_raw,mapsize_global).astype(np.single)
            del gswe_raw
            dset.close()
            
            # Fill gaps in GSWE fracwater (at high latitudes) with HILDA+ fracwater
            fracwater[np.isnan(fracwater)] = fracwater_hilda[np.isnan(fracwater)]            
            fracwater_all[:,:,jj] = fracwater
            
        # Compute and save climatology
        fracwater_clim = np.mean(fracwater_all,axis=2)
        np.savez_compressed(os.path.join(config['output_folder'],'JRC_GSWE',str(month).zfill(2)+'.npz'),data=fracwater_clim)
        print("Time elapsed is "+str(time.time()-t0)+" sec")                
        
    pdb.set_trace()
    
if __name__ == '__main__':
    main()