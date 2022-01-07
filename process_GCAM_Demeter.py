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
    
    # Process the following scenarios, years, and classes
    # See table 2 of Chen et al. (2020; https://doi.org/10.1038/s41597-020-00669-x) for the legend
    scenarios = ['ssp1_rcp26', 'ssp1_rcp45', 'ssp1_rcp60', 'ssp2_rcp26',  'ssp2_rcp45',  'ssp2_rcp60',  'ssp3_rcp45',  'ssp3_rcp60',  'ssp4_rcp26',  'ssp4_rcp45',  'ssp4_rcp60',  'ssp5_rcp26',  'ssp5_rcp45',  'ssp5_rcp60', 'ssp5_rcp85']
    years = np.arange(2015,2105,5)
    klassen = {}
    klassen['forest'] = np.arange(1,9)                      # Forest
    klassen['rice'] = np.array([24])                        # Irrigated rice
    klassen['irrigation'] = np.array([16,18,20,22,26,28])   # Irrigated crops (no rice)   
    
    # Loop through scenarios
    for scenario in scenarios:
        
        # Create output folder
        if os.path.isdir(os.path.join(config['output_folder'],'GCAM_Demeter',scenario))==False:
            os.mkdir(os.path.join(config['output_folder'],'GCAM_Demeter',scenario))
            
        # Loop through years and classes (forest, rice, and irrigation)
        for year in years:
            print('Processing '+scenario+' '+str(year))
            t0 = time.time()
            dset = Dataset(os.path.join(config['gcam_demeter_folder'],'GCAM_Demeter_LU_'+scenario+'_modelmean_'+str(year)+'.nc'))
            for klas in klassen.keys():
                data_sum = np.zeros(mapsize_global,dtype=np.single)
                for ii in np.arange(len(klassen[klas])):
                    type = klassen[klas][ii]
                    data = np.array(dset.variables['PFT'+str(type)][:]).transpose()
                    data = imresize_mean(data,mapsize_global)
                    data_sum += data/100
                np.savez_compressed(os.path.join(config['output_folder'],'GCAM_Demeter',scenario,'frac'+klas+'_'+str(year)+'.npz'),data=data_sum)
            dset.close()
            print("Time elapsed is "+str(time.time()-t0)+" sec")
            
    pdb.set_trace()
       
    
if __name__ == '__main__':
    main()