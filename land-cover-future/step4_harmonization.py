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
import gc

# Load configuration file
config = load_config(sys.argv[1])

def main():
    print('===============================================================================')

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
    
    # List of years
    years = np.arange(config['year_start'],config['year_end']+1).astype(int)
    
    # List of years with GCAM-Demeter data
    gcam_demeter_scenarios = glob.glob(os.path.join(config['output_folder'],'step2_GCAM_Demeter','*'))
    gcam_demeter_files = glob.glob(os.path.join(gcam_demeter_scenarios[0],'fracforest*'))
    gcam_demeter_years = np.array([int(os.path.basename(gcam_demeter_file)[-8:-4]) for gcam_demeter_file in gcam_demeter_files])
    
    # List of years with Chen et al. (2020) urban data
    chen_2020_scenarios = glob.glob(os.path.join(config['output_folder'],'step1_Chen_2020_urban','*'))
    chen_2020_files = glob.glob(os.path.join(chen_2020_scenarios[0],'fracsealed*'))
    chen_2020_years = np.array([int(os.path.basename(chen_2020_file)[-8:-4]) for chen_2020_file in chen_2020_files])
    
    
    ############################################################################
    #   Maps are only available every decade (every five years for GCAM-Demeter).
    #   Here we loop through decades to calculate fractions (no temporal 
    #   interpolation).
    ############################################################################
    
    scenarios = glob.glob(os.path.join(config['output_folder'],'step2_GCAM_Demeter','*'))
    scenarios = [os.path.basename(scenario) for scenario in scenarios]
    
    for scenario in scenarios:
        scenario_ssp = scenario[:4]
        scenario_rcp = scenario[5:]
        
        # Create output folder
        if os.path.isdir(os.path.join(config['output_folder'],'step4_harmonization',scenario))==False:
            os.makedirs(os.path.join(config['output_folder'],'step4_harmonization',scenario))
            
        # Loop through years and months
        for year in years:       
            for month in np.arange(1,13):
                print('-------------------------------------------------------------------------------')
                print('Scenario: '+scenario+' Year: '+str(year)+' Month: '+str(month))
                t0 = time.time()
                frac = {}
                
                print('Loading Chen et al. (2020) fracsealed data (interpolating between years)')
                frac['sealed'] = load_data_interp(year,os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario_ssp,'fracsealed_*')).clip(0,1)
                
                print('Loading GSWE fracwater data (monthly climatology)')
                frac['water'] = np.load(os.path.join(config['output_folder'],'step3_JRC_GSWE',str(month).zfill(2)+'.npz'))['data'].clip(0,1)
                
                print('Reduce fracwater if fracsealed+fracwater >1 (fracsealed overrides fracwater)')
                totals = frac['sealed']+frac['water']
                mask = totals>1
                frac['water'][mask] = frac['water'][mask]-(totals[mask]-1)
                frac['water'] = frac['water'].clip(0,1)
                
                print('Computing fracother as residual of fracsealed and fracwater')
                fracother_init = 1-frac['sealed']-frac['water']
                fracother_init = fracother_init.clip(0,1)
                
                print('Loading GCAM-Demeter fracforest, fracirrigation, and fracrice data (interpolating between years)')
                frac['forest'] = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracforest_*'))).clip(0,1)
                frac['irrigation'] = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracirrigation_*'))).clip(0,1)
                frac['rice'] = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracrice_*'))).clip(0,1)
                
                print('Reducing fracforest, fracirrigation, and fracrice if sum exceeds fracother')
                totals = frac['forest']+frac['rice']+frac['irrigation']
                mask = totals>fracother_init
                excess = totals-fracother_init
                frac['forest'][mask] = frac['forest'][mask]-excess[mask]*frac['forest'][mask]/totals[mask]
                frac['rice'][mask] = frac['rice'][mask]-excess[mask]*frac['rice'][mask]/totals[mask]
                frac['irrigation'][mask] = frac['irrigation'][mask]-excess[mask]*frac['irrigation'][mask]/totals[mask]
                
                print('Reducing precision to save space (rounding down to avoid sum >1)')
                precision = 2
                factor = 10**2
                for key in frac.keys(): 
                    frac[key] = np.floor(frac[key]*factor)/factor
                    frac[key] = frac[key].clip(0,1)
                
                print('Recalculating fracother as residual')
                frac['other'] = 1-frac['water']-frac['sealed']-frac['forest']-frac['rice']-frac['irrigation']
                
                print('Subsetting to template map area')
                for key in frac.keys():
                    frac[key] = frac[key][row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]

                print('Saving in Numpy format')
                for key in frac.keys():
                    np.savez_compressed(os.path.join(config['output_folder'],'step4_harmonization',scenario,str(year)+str(month).zfill(2)+'_frac'+key+'.npz'),data=frac[key])
                    
                print("Time elapsed is "+str(time.time()-t0)+" sec")


        ############################################################################
        #   Convert to netCDF. least_significant_digit option (useful to conserve
        #   disk space) not used as it caused issues when running LISFLOOD. 
        ############################################################################

        for key in frac.keys(): 
            varname = 'frac'+key

            file = os.path.join(config['output_folder'],'step4_harmonization',scenario,varname+'.nc')            
            
            if os.path.isfile(file):
                os.remove(file)
                
            ncfile = Dataset(file, 'w', format='NETCDF4')
            ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

            ncfile.createDimension('lon', len(template_lon))
            ncfile.createDimension('lat', len(template_lat))
            ncfile.createDimension('time', None)

            ncfile.createVariable('lon', 'f8', ('lon',))
            ncfile.variables['lon'][:] = template_lon
            ncfile.variables['lon'].units = 'degrees_east'
            ncfile.variables['lon'].long_name = 'longitude'

            ncfile.createVariable('lat', 'f8', ('lat',))
            ncfile.variables['lat'][:] = template_lat
            ncfile.variables['lat'].units = 'degrees_north'
            ncfile.variables['lat'].long_name = 'latitude'

            ncfile.createVariable('time', 'f8', 'time')
            ncfile.variables['time'].units = 'days since 1979-01-02 00:00:00'
            ncfile.variables['time'].long_name = 'time'
            ncfile.variables['time'].calendar = 'proleptic_gregorian'

            ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999)
            ncfile.variables[varname].units = 'fraction'

            for year in years:
                print('-------------------------------------------------------------------------------')
                print('Saving to netCDF Scenario: '+scenario+' Year: '+str(year)+' Var: '+varname)
                t0 = time.time()
                
                for month in np.arange(1,13):
                    data = np.load(os.path.join(config['output_folder'],'step4_harmonization',scenario,str(year)+str(month).zfill(2)+'_'+varname+'.npz'))['data']
                    index = (year-config['year_start'])*12+month-1
                    ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
                    ncfile.variables[varname][index,:,:] = data                    
                
                print("Time elapsed is "+str(time.time()-t0)+" sec")
                    
            ncfile.close()


        ############################################################################
        #   Verify that sum of fractions is 1
        ############################################################################

        print('-------------------------------------------------------------------------------')
        print('Sum verification')
        
        sum = np.zeros(mapsize_template)
        for key in frac.keys(): 
            varname = 'frac'+key
            file = os.path.join(config['output_folder'],'step4_harmonization',scenario,varname+'.nc')
            ncfile = Dataset(file)
            data = np.array(ncfile.variables[varname][0,:,:])
            sum += data

        print('Max sum of fractions: '+str(np.nanmax(sum)))
        print('Min sum of fractions: '+str(np.nanmin(sum)))
            
    pdb.set_trace()
        
        
if __name__ == '__main__':
    main()