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
import traceback

# Load configuration file 
config = load_config(sys.argv[1])

def main():
    
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
    
    # Loop through scenarios and models
    files = glob.glob(os.path.join(config['meteo_data_folder'],'*.nc'))
    scenarios = np.unique([os.path.basename(x).split('_')[3] for x in files])
    for scenario in scenarios:
        files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+scenario+'*.nc'))
        files = sorted(files)
        models = np.unique([os.path.basename(x).split('_')[0] for x in files])
        for model in models:           
            print('===============================================================================')            
            print('scenario: '+scenario+' model: '+model)
            print('-------------------------------------------------------------------------------')
            
            # What is time span of output?
            if scenario=='historical': 
                year_start,year_end = 1979,2014
            else: 
                year_start,year_end = 2015,2100
            out_dates_dly = pd.date_range(start=datetime(year_start,1,1), end=datetime(year_end+1,1,1)-pd.Timedelta(days=1), freq='D')

            # Check if output already exists in scratch folder or output folder
            scratchoutdir = os.path.join(config['scratch_folder'],scenario,model)
            finaloutdir = os.path.join(config['output_folder'],scenario,model)            
            if (config['delete_existing']==False) & ((os.path.isfile(os.path.join(scratchoutdir,'ta.nc'))==True) | (os.path.isfile(os.path.join(finaloutdir,'ta.nc'))==True)):
                print('Already processed, skipping this model')
                continue
                
            # Check if all input variables are present
            varnames = ['tas','tasmin','tasmax','hurs','sfcwind','ps','rsds','rlds','pr']
            nfiles = {}
            for varname in varnames:
                files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+model+'*'+scenario+'*_'+varname+'_*.nc'))    
                nfiles[varname] = len(files)            
            nfiles_arr = np.array(list(nfiles.values()))
            if (np.max(nfiles_arr)==0) | (any(nfiles_arr<np.max(nfiles_arr))):
                print('Not all input variables present, skipping this model')
                print(nfiles)
                continue
                
            # Initialize output files
            if os.path.isdir(scratchoutdir)==False:
                os.makedirs(scratchoutdir)
            ncfile_ta = initialize_netcdf(os.path.join(scratchoutdir,'ta.nc'),template_lat,template_lon,'ta','degree_Celsius',1)
            ncfile_pr = initialize_netcdf(os.path.join(scratchoutdir,'pr.nc'),template_lat,template_lon,'pr','mm d-1',1)
            ncfile_et = initialize_netcdf(os.path.join(scratchoutdir,'et.nc'),template_lat,template_lon,'et','mm d-1',1)
            ncfile_ew = initialize_netcdf(os.path.join(scratchoutdir,'ew.nc'),template_lat,template_lon,'ew','mm d-1',1)
            ncfile_es = initialize_netcdf(os.path.join(scratchoutdir,'es.nc'),template_lat,template_lon,'es','mm d-1',1)

            # Loop over input files (MFDataset doesn't work properly)
            files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+model+'*'+scenario+'*_tas_*.nc'))
            for file in files:
                file_year_start = int(os.path.basename(file).split('_')[7])
                file_year_end = int(os.path.basename(file).split('_')[8][:-3])
                file_dates_dly = pd.date_range(start=datetime(file_year_start,1,1), end=datetime(file_year_end+1,1,1)-pd.Timedelta(days=1), freq='D')
                hits = np.sum((out_dates_dly.year>=file_year_start) & (out_dates_dly.year<=file_year_end))
                if hits==0:
                    continue
                
                # Open input files
                print('Processing '+os.path.basename(file))
                t0 = time.time()
                dset_tmean = netCDF4.Dataset(file,diskless=True) # K
                dset_tmin = netCDF4.Dataset(file.replace('tas_','tasmin_'),diskless=True) # K
                dset_tmax = netCDF4.Dataset(file.replace('tas_','tasmax_'),diskless=True) # K
                dset_relhum = netCDF4.Dataset(file.replace('tas_','hurs_'),diskless=True) # kg/kg
                dset_wind = netCDF4.Dataset(file.replace('tas_','sfcwind_'),diskless=True) # m/s
                dset_pres = netCDF4.Dataset(file.replace('tas_','ps_'),diskless=True) # Pa
                dset_swd = netCDF4.Dataset(file.replace('tas_','rsds_'),diskless=True) # W/m2
                dset_lwd = netCDF4.Dataset(file.replace('tas_','rlds_'),diskless=True) # W/m2
                dset_pr = netCDF4.Dataset(file.replace('tas_','pr_'),diskless=True) # mm/d
                    
                # Loop over days of input file
                for ii in np.arange(len(file_dates_dly)):
                    if file_dates_dly[ii] not in out_dates_dly:
                        continue
                    #print('Processing '+os.path.basename(file)+', ii: '+str(ii)+', time stamp: '+datetime.now().strftime("%d/%m/%Y, %H:%M:%S")+')')                    
                        
                    # Read data from input files
                    data = {}
                    index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]                            
                    data['tmean'] = np.array(dset_tmean.variables['tas'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['tmin'] = np.array(dset_tmin.variables['tasmin'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['tmax'] = np.array(dset_tmax.variables['tasmax'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['relhum'] = np.array(dset_relhum.variables['hurs'][ii,:,:],dtype=np.single) # %
                    data['wind'] = np.array(dset_wind.variables['sfcwind'][ii,:,:],dtype=np.single)*0.75 # m/s (factor 0.75 to translate from 10-m to 2-m height)
                    data['pres'] = np.array(dset_pres.variables['ps'][ii,:,:],dtype=np.single)/100 # mbar
                    data['swd'] = np.array(dset_swd.variables['rsds'][ii,:,:],dtype=np.single) # W/m2
                    data['lwd'] = np.array(dset_lwd.variables['rlds'][ii,:,:],dtype=np.single) # W/m2
                    data['pr'] = np.array(dset_pr.variables['pr'][ii,:,:],dtype=np.single)*86400 # mm/d
                    
                    # Switch eastern and western hemispheres
                    #for key in data.keys():
                    #    data[key] = np.roll(data[key],int(data[key].shape[1]/2),axis=1)
                
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
                    doy = int(datetime.strftime(file_dates_dly[ii],'%j'))
                    pet = potential_evaporation(data,albedo,factor,doy,lat_template,elev_template)
                    
                    # Write data to output netCDFs
                    time_value = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400                    
                    index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]
                    ncfile_pr.variables['time'][index] = time_value
                    ncfile_pr.variables['pr'][index,:,:] = data['pr']
                    ncfile_ta.variables['time'][index] = time_value
                    ncfile_ta.variables['ta'][index,:,:] = data['tmean']
                    ncfile_et.variables['time'][index] = time_value
                    ncfile_et.variables['et'][index,:,:] = pet['et']
                    ncfile_ew.variables['time'][index] = time_value
                    ncfile_ew.variables['ew'][index,:,:] = pet['ew']
                    ncfile_es.variables['time'][index] = time_value
                    ncfile_es.variables['es'][index,:,:] = pet['es']
                        
                    # Generate figures to verify output
                    if ii==0:
                        makefig('figures','et',pet['et'],0,12)
                        makefig('figures','ew',pet['ew'],0,12)
                        makefig('figures','es',pet['es'],0,12)
                        for key in data.keys():
                            makefig('figures',key,data[key],np.min(data[key]),np.max(data[key]))
                        makefig('figures','elev_template',elev_template,0,6000)
                
                # Close input files
                dset_tmean.close()
                dset_tmin.close()
                dset_tmax.close()
                dset_relhum.close()
                dset_wind.close()
                dset_pres.close()
                dset_swd.close()
                dset_lwd.close()
                dset_pr.close()
                
                print("Time elapsed is "+str(time.time()-t0)+" sec")
            
            # Close output files
            ncfile_pr.close()
            ncfile_ta.close()
            ncfile_et.close()
            ncfile_ew.close()
            ncfile_es.close()

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