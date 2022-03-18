#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "March 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
import netCDF4
from tools import *
from skimage.transform import resize
import matplotlib.pyplot as plt
import rasterio
import shutil
from datetime import timedelta

# Load configuration file 
config = load_config(sys.argv[1])

def main():

    # Output dates
    year_start,year_end = 1979,1979
    out_dates_dly = pd.date_range(start=datetime(year_start,1,1), end=datetime(year_end+1,1,1)-pd.Timedelta(days=1), freq='D')
    
    # Load template map
    dset = netCDF4.Dataset(config['templatemap_path'])
    template_lat = np.array(dset.variables['lat'][:])
    template_lon = np.array(dset.variables['lon'][:])
    template_res = template_lon[1]-template_lon[0]
    varname = list(dset.variables.keys())[-1]
    template_np = np.array(dset.variables[varname][:])

    # Determine map sizes
    mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
    mapsize_template = template_np.shape
    row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)
    
    # Load elevation data, append zeros to top and bottom to make global, and resample to template resolution
    elev = np.zeros((21600,43200),dtype=np.single)
    src = rasterio.open(os.path.join(config['gmted2010_folder'],'elevation_1KMmn_GMTEDmn.tif'))
    elev[720:17520,:] = src.read(1)
    src.close()
    elev_global = imresize_mean(elev,mapsize_global)
    elev_template = elev_global[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]

    # Prepare temperature and air pres downscaling
    tmp = imresize_mean(elev_global,(360,720))
    tmp = resize(tmp,mapsize_global,order=1,mode='edge',anti_aliasing=False)
    elev_delta = elev_global-tmp
    temp_delta = -6.5*elev_delta/1000 # Simple 6.5 degrees C/km lapse rate
    pres_delta = 1013*((((293-0.0065*elev_global)/293)**5.26)-(((293-0.0065*tmp)/293)**5.26)) # Allen et al. (1994) equation 7
    lat_global = np.repeat(np.resize(np.arange(90-template_res/2,-90-template_res/2,-template_res),(1800,1)),mapsize_global[1],axis=1)
    lat_template = lat_global[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
    
    # Check if output already exists in scratch folder or output folder
    scratchoutdir = os.path.join(config['scratch_folder'],'MSWX_MSWEP_reanalysis')
    finaloutdir = os.path.join(config['output_folder'],'MSWX_MSWEP_reanalysis')            
    if (config['delete_existing']==False) & (os.path.isfile(os.path.join(finaloutdir,'ta.nc'))==True):
        print('Already processed, skipping this model')
        quit()
        
    # Initialize output files
    if os.path.isdir(scratchoutdir)==False:
        os.makedirs(scratchoutdir)
    ncfile = {}
    ncfile['ta'] = initialize_netcdf(os.path.join(scratchoutdir,'ta.nc'),template_lat,template_lon,'ta','degree_Celsius',1)
    ncfile['pr'] = initialize_netcdf(os.path.join(scratchoutdir,'pr.nc'),template_lat,template_lon,'pr','mm d-1',1)
    ncfile['et'] = initialize_netcdf(os.path.join(scratchoutdir,'et.nc'),template_lat,template_lon,'et','mm d-1',1)
    ncfile['ew'] = initialize_netcdf(os.path.join(scratchoutdir,'ew.nc'),template_lat,template_lon,'ew','mm d-1',1)
    ncfile['es'] = initialize_netcdf(os.path.join(scratchoutdir,'es.nc'),template_lat,template_lon,'es','mm d-1',1)

    # Loop over days of input file    
    for ii in np.arange(len(out_dates_dly)):
        current_date = out_dates_dly[ii]
        
        # Data missing on some days
        if (current_date==datetime(1979,1,1)) | (current_date==datetime(1979,1,31)):
            current_date = current_date+timedelta(days=1)
            
        # Open input files        
        print('Processing '+datetime.strftime(current_date,'%Y%j'))
        t0 = time.time()
        invars = {'Temp','Tmin','Tmax','RelHum','Wind','Pres','SWd','LWd','P'}
        dset = {}
        for invar in invars:
            file = os.path.join(config['mswx_folder'],'Past',invar,'Daily',datetime.strftime(current_date,'%Y%j')+'.nc')
            dset[invar] = netCDF4.Dataset(file.replace('Temp',invar))
        file_pr1 = os.path.join(config['mswep_folder'],'Past','Daily',datetime.strftime(current_date,'%Y%j')+'.nc')
        file_pr2 = os.path.join(config['mswep_folder'],'NRT','Daily',datetime.strftime(current_date,'%Y%j')+'.nc')
        if os.path.isfile(file_pr1):
            dset['P'] = netCDF4.Dataset(file_pr1) # mm/d
        else:
            dset['P'] = netCDF4.Dataset(file_pr2) # mm/d            
            
        # Read data from input files
        data = {}
        data['tmean'] = np.squeeze(np.array(dset['Temp'].variables['air_temperature'])) # degrees C
        data['tmin'] = np.squeeze(np.array(dset['Tmin'].variables['air_temperature'])) # degrees C
        data['tmax'] = np.squeeze(np.array(dset['Tmax'].variables['air_temperature'])) # degrees C
        data['relhum'] = np.squeeze(np.array(dset['RelHum'].variables['relative_humidity'])) # %
        data['wind'] = np.squeeze(np.array(dset['Wind'].variables['wind_speed'])) # m/s
        data['pres'] = np.squeeze(np.array(dset['Pres'].variables['surface_pressure']))/100 # mbar
        data['swd'] = np.squeeze(np.array(dset['SWd'].variables['downward_shortwave_radiation'])) # W/m2
        data['lwd'] = np.squeeze(np.array(dset['LWd'].variables['downward_longwave_radiation'])) # W/m2
        data['pr'] = np.squeeze(np.array(dset['P'].variables['precipitation'])) # mm/d
   
        # Simple lapse rate downscaling of temperature and air pressure, nearest-neighbor resampling of other vars
        for key in data.keys():
            if (key=='tmean') | (key=='tmin') | (key=='tmax'):
                data[key] = resize(data[key],mapsize_global,order=1,mode='edge',anti_aliasing=False)
                data[key] = data[key]+temp_delta
            if key=='pres':
                data[key] = resize(data[key],mapsize_global,order=1,mode='edge',anti_aliasing=False)
                data[key] = data[key]+pres_delta
            else:
                data[key] = resize(data[key],mapsize_global,order=0,mode='edge',anti_aliasing=False)

        # Subset data to template region
        for key in data.keys():
            data[key] = data[key][row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        
        # Compute potential evaporation
        albedo = {'et':0.23,'ew':0.05,'es':0.15}
        factor = {'et':1,'ew':0.5,'es':0.75}
        doy = int(datetime.strftime(current_date,'%j'))
        pet = potential_evaporation(data,albedo,factor,doy,lat_template,elev_template)
        
        # Write data to output netCDFs
        time_value = (current_date-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400                    
        index = np.where(out_dates_dly==current_date)[0][0]
        ncfile['pr'].variables['time'][index] = time_value
        ncfile['pr'].variables['pr'][index,:,:] = data['pr']
        ncfile['ta'].variables['time'][index] = time_value
        ncfile['ta'].variables['ta'][index,:,:] = data['tmean']
        ncfile['et'].variables['time'][index] = time_value
        ncfile['et'].variables['et'][index,:,:] = pet['et']
        ncfile['ew'].variables['time'][index] = time_value
        ncfile['ew'].variables['ew'][index,:,:] = pet['ew']
        ncfile['es'].variables['time'][index] = time_value
        ncfile['es'].variables['es'][index,:,:] = pet['es']
            
        # Generate figures to verify output        
        if ii==0:
            makefig('figures','et',pet['et'],0,12)
            makefig('figures','ew',pet['ew'],0,12)
            makefig('figures','es',pet['es'],0,12)
            for key in data.keys():
                makefig('figures',key,data[key],np.min(data[key]),np.max(data[key]))
            makefig('figures','elev_template',elev_template,0,6000)
        
        # Close input files
        for key in dset.keys():
            dset[key].close()
        
        print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    # Close output files
    for key in ncfile.keys():
        ncfile[key].close()

    # Move output from scratch folder to output folder
    print('-------------------------------------------------------------------------------')
    if os.path.isdir(finaloutdir)==False:
        os.makedirs(finaloutdir)
    for file in glob.glob(os.path.join(scratchoutdir,'*')):
        t0 = time.time()
        print('Moving '+os.path.basename(file)+' ('+str(round(os.path.getsize(file)/10**9))+' GB) to '+finaloutdir)
        shutil.copy(file, finaloutdir)
        print("Time elapsed is "+str(time.time()-t0)+" sec")

if __name__ == '__main__':
    main()
