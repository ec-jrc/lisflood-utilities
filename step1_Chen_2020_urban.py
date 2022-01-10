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
import rasterio

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
    
    # Process the following scenarios and years
    scenarios = ['ssp1', 'ssp2', 'ssp3', 'ssp4',  'ssp5']
    years = np.arange(2020,2110,10)
    
    # Loop through scenarios
    for scenario in scenarios:
        
        # Create output folder
        if os.path.isdir(os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario))==False:
            os.makedirs(os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario))
            
        # Loop through years
        for year in years:
            print('Processing '+scenario+' '+str(year))
            t0 = time.time()
            src = rasterio.open(os.path.join(config['chen_2020_urban_folder'],scenario.upper(),'global_'+scenario.upper()+'_'+str(year)+'.tif'))
            data = src.read(1).astype(np.single)
            src.close()
            data = data-1
            data[data<0] = 0            
            data = imresize_mean(data,mapsize_global)            
            np.savez_compressed(os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario,'fracsealed_'+str(year)+'.npz'),data=data)
            print("Time elapsed is "+str(time.time()-t0)+" sec")
            
    pdb.set_trace()
       
    
if __name__ == '__main__':
    main()