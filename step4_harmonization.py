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
    
    vars = ['fracwater','fracforest','fracsealed','fracrice','fracirrigation','fracother']
    
    
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
                
                print('Loading Chen et al. (2020) fracsealed data (interpolating between years)')
                fracsealed = load_data_interp(year,os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario_ssp,'fracsealed_*')).clip(0,1)
                
                print('Loading GSWE fracwater data (monthly climatology)')
                fracwater = np.load(os.path.join(config['output_folder'],'step3_JRC_GSWE',str(month).zfill(2)+'.npz'))['data'].clip(0,1)
                
                print('Reduce fracwater if fracsealed+fracwater >1 (fracsealed overrides fracwater)')
                totals = fracsealed+fracwater
                mask = totals>1
                fracwater[mask] = fracwater[mask]-(totals[mask]-1)
                fracwater = fracwater.clip(0,1)
                
                print('Computing fracother as residual of fracsealed and fracwater')
                fracother_init = 1-fracsealed-fracwater
                fracother_init = fracother_init.clip(0,1)
                
                print('Loading GCAM-Demeter fracforest, fracirrigation, and fracrice data (interpolating between years)')
                fracforest = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracforest_*'))).clip(0,1)
                fracirrigation = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracirrigation_*'))).clip(0,1)
                fracrice = fill(load_data_interp(year,os.path.join(config['output_folder'],'step2_GCAM_Demeter',scenario,'fracrice_*'))).clip(0,1)
                
                print('Reducing fracforest, fracirrigation, and fracrice if sum exceeds fracother')
                totals = fracforest+fracrice+fracirrigation
                mask = totals>fracother_init
                excess = totals-fracother_init
                fracforest[mask] = fracforest[mask]-excess[mask]*fracforest[mask]/totals[mask]
                fracrice[mask] = fracrice[mask]-excess[mask]*fracrice[mask]/totals[mask]
                fracirrigation[mask] = fracirrigation[mask]-excess[mask]*fracirrigation[mask]/totals[mask]
                
                print('Reducing precision to save space (round down to avoid sum >1)')
                precision = 2
                factor = 10**2
                fracwater = np.floor(fracwater*factor)/factor
                fracsealed = np.floor(fracsealed*factor)/factor
                fracforest = np.floor(fracforest*factor)/factor
                fracrice = np.floor(fracrice*factor)/factor
                fracirrigation = np.floor(fracirrigation*factor)/factor
                
                print('Recalculating fracother as residual')
                fracother = 1-fracwater-fracsealed-fracforest-fracrice-fracirrigation
                
                print('Subsetting data to template map area')
                fracwater = fracwater[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                fracforest = fracforest[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                fracsealed = fracsealed[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                fracrice = fracrice[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                fracirrigation = fracirrigation[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                fracother = fracother[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]

                print('Saving data in Numpy format')
                t0 = time.time()
                
                sum = np.zeros(mapsize_template)
                
                for vv in np.arange(len(vars)):            
                    data = eval(vars[vv])
                    data = np.round(data,precision)
                    np.savez_compressed(os.path.join(config['output_folder'],'step4_harmonization',scenario,str(year)+str(month).zfill(2)+'_'+vars[vv]),data=data)

                    sum += data
                
                print('Max sum of fractions: '+str(np.nanmax(sum)))
                print('Min sum of fractions: '+str(np.nanmin(sum)))

                print("Time elapsed is "+str(time.time()-t0)+" sec")


        ############################################################################
        #   Convert to netCDF. least_significant_digit option (useful to conserve
        #   disk space) not used as it caused issues when running LISFLOOD. 
        ############################################################################

        for vv in np.arange(len(vars)):

            file = os.path.join(config['output_folder'],'step4_harmonization',scenario,vars[vv]+'.nc')
            varname = vars[vv]
            
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

            ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999) #, least_significant_digit=9
            ncfile.variables[varname].units = 'fraction'

            for year in years:
                for month in np.arange(1,13):
                    print('-------------------------------------------------------------------------------')
                    print('Saving to netCDF scenario: '+scenario+' var: '+varname+' year: '+str(year)+' month: '+str(month))
                    t0 = time.time()
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
        for vv in np.arange(len(vars)):
            varname = vars[vv]        
            file = os.path.join(config['output_folder'],'step4_harmonization',scenario,varname+'.nc')    
            ncfile = Dataset(file)
            data = np.array(ncfile.variables[varname][0,:,:])
            sum += data

        print('Max sum of fractions: '+str(np.nanmax(sum)))
        print('Min sum of fractions: '+str(np.nanmin(sum)))
            
    pdb.set_trace()
        
        
if __name__ == '__main__':
    main()