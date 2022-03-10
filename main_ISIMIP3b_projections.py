#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
import netCDF4
from tools import *
from skimage.transform import resize
import matplotlib.pyplot as plt
import rasterio
import shutil

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
    temp_delta = -6.5*elev_delta/1000
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

            # What is time span of output?
            if scenario=='historical':
                year_start = 1979
                year_end = 2014
            else:
                year_start = 2015
                year_end = 2100
            out_dates_dly = pd.date_range(start=datetime(year_start,1,1), end=datetime(year_end+1,1,1)-pd.Timedelta(days=1), freq='D')

            '''
                
            # Loop over variables
            vars = np.array([\
                ['ta','tas_',1,-273.15,'degree_Celsius',1],\
                ['pr','pr_',86400,0,'mm d-1',1],\
                ['et','tas_',1,0,'mm d-1',1],\
                ['ew','tas_',1,0,'mm d-1',1],\
                ['es','tas_',1,0,'mm d-1',1],\
                ])
                
                
            for vv in np.arange(len(vars)):
                
                # Initialize netCDF
                outdir = os.path.join(config['output_folder'],scenario,model)
                if os.path.isdir(outdir)==False:
                    os.makedirs(outdir)
                outfile = os.path.join(outdir,vars[vv,0]+'.nc')
                if os.path.isfile(outfile):
                    os.remove(outfile)                
                ncfile = initialize_netcdf(outfile,template_lat,template_lon,vars[vv,0],vars[vv,4],int(vars[vv,5]))
                

                #------------------------------------------------------------------------------
                #   Precipitation and temperature
                #------------------------------------------------------------------------------                            
                
                if (vars[vv,0]=='pr') | (vars[vv,0]=='ta'):
                
                    # Loop over input files (MFDataset unfortunately doesn't work properly)
                    files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+model+'*'+scenario+'*'+vars[vv,1]+'*.nc'))
                    for file in files:
                        file_year_start = int(os.path.basename(file).split('_')[7])
                        file_year_end = int(os.path.basename(file).split('_')[8][:-3])
                        file_dates_dly = pd.date_range(start=datetime(file_year_start,1,1), end=datetime(file_year_end+1,1,1)-pd.Timedelta(days=1), freq='D')                        
                        hits = np.sum((out_dates_dly.year>=file_year_start) & (out_dates_dly.year<=file_year_end))
                        if hits==0:
                            continue
                        
                        # Open input file
                        print('Processing '+file)
                        t0 = time.time() 
                        dset = netCDF4.Dataset(file)
                            
                        # Loop over days of input file
                        for ii in np.arange(10): #np.arange(len(file_dates_dly)):
                            if file_dates_dly[ii] not in out_dates_dly:
                                continue
                            
                            # Read data from input file, fix units, and shift map
                            index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]
                            data = np.array(dset.variables[vars[vv,1][:-1]][ii,:,:],dtype=np.single)
                            data = float(vars[vv,2])*data+float(vars[vv,3])
                            data = np.roll(data,int(data.shape[1]/2),axis=1)
                        
                            # Simple lapse rate downscaling of temperature, nearest-neighbor resampling of other vars
                            if vars[vv,0]=='ta':
                                data = resize(data,mapsize_global,order=1,mode='edge',anti_aliasing=False)
                                data = data+temp_delta
                            else:
                                data = resize(data,mapsize_global,order=0,mode='edge',anti_aliasing=False)
                                                    
                            # Subset data and write to netCDF
                            data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                            ncfile.variables['time'][index] = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400
                            ncfile.variables[vars[vv,0]][index,:,:] = data                    
                        
                        print("Time elapsed is "+str(time.time()-t0)+" sec")
                    
                    
                #------------------------------------------------------------------------------
                #   Potential evaporation
                #------------------------------------------------------------------------------                            
                
                if (vars[vv,0]=='et') | (vars[vv,0]=='ew') | (vars[vv,0]=='es'):
            '''
            
            # Initialize netCDFs
            outdir = os.path.join(config['output_folder'],scenario,model)
            if os.path.isdir(outdir)==False:
                os.makedirs(outdir)
            else:
                shutil.rmtree(outdir,ignore_errors=True)
            ncfile_ta = initialize_netcdf(os.path.join(outdir,'ta.nc'),template_lat,template_lon,'ta','degree_Celsius',1)
            ncfile_pr = initialize_netcdf(os.path.join(outdir,'pr.nc'),template_lat,template_lon,'pr','mm d-1',1)
            ncfile_et = initialize_netcdf(os.path.join(outdir,'et.nc'),template_lat,template_lon,'et','mm d-1',1)
            ncfile_ew = initialize_netcdf(os.path.join(outdir,'ew.nc'),template_lat,template_lon,'ew','mm d-1',1)
            ncfile_es = initialize_netcdf(os.path.join(outdir,'es.nc'),template_lat,template_lon,'es','mm d-1',1)

            # Loop over input files (MFDataset unfortunately doesn't work properly)
            files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+model+'*'+scenario+'*_tas_*.nc'))
            for file in files:
                file_year_start = int(os.path.basename(file).split('_')[7])
                file_year_end = int(os.path.basename(file).split('_')[8][:-3])
                file_dates_dly = pd.date_range(start=datetime(file_year_start,1,1), end=datetime(file_year_end+1,1,1)-pd.Timedelta(days=1), freq='D')
                hits = np.sum((out_dates_dly.year>=file_year_start) & (out_dates_dly.year<=file_year_end))
                if hits==0:
                    continue
                
                # Open input files
                print('Processing '+file)
                t0 = time.time() 
                dset_tmean = netCDF4.Dataset(file) # K
                dset_tmin = netCDF4.Dataset(file.replace('tas_','tasmin_')) # K
                dset_tmax = netCDF4.Dataset(file.replace('tas_','tasmax_')) # K
                dset_relhum = netCDF4.Dataset(file.replace('tas_','hurs_')) # kg/kg
                dset_wind = netCDF4.Dataset(file.replace('tas_','sfcwind_')) # m/s
                dset_pres = netCDF4.Dataset(file.replace('tas_','ps_')) # Pa
                dset_swd = netCDF4.Dataset(file.replace('tas_','rsds_')) # W/m2
                dset_lwd = netCDF4.Dataset(file.replace('tas_','rlds_')) # W/m2
                dset_pr = netCDF4.Dataset(file.replace('tas_','pr_')) # mm/d
                
                # Loop over days of input file
                #for ii in np.arange(10):
                for ii in np.arange(len(file_dates_dly)):
                    if file_dates_dly[ii] not in out_dates_dly:
                        continue
                    
                    # Read data from input files
                    data = {}
                    index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]                            
                    data['tmean'] = np.array(dset_tmean.variables['tas'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['tmin'] = np.array(dset_tmin.variables['tasmin'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['tmax'] = np.array(dset_tmax.variables['tasmax'][ii,:,:],dtype=np.single)-273.15 # degrees C
                    data['relhum'] = np.array(dset_relhum.variables['hurs'][ii,:,:],dtype=np.single) # %
                    data['wind'] = np.array(dset_wind.variables['sfcwind'][ii,:,:],dtype=np.single)*0.75 # m/s (factor 0.75 to convert from 10-m to 2-m height)
                    data['pres'] = np.array(dset_pres.variables['ps'][ii,:,:],dtype=np.single)/100 # mbar
                    data['swd'] = np.array(dset_swd.variables['rsds'][ii,:,:],dtype=np.single) # W/m2
                    data['lwd'] = np.array(dset_lwd.variables['rlds'][ii,:,:],dtype=np.single) # W/m2
                    data['pr'] = np.array(dset_pr.variables['pr'][ii,:,:],dtype=np.single)*86400 # mm/d
                    #for key in data.keys():
                    #    data[key] = np.roll(data[key],int(data[key].shape[1]/2),axis=1)
                
                    # Simple lapse rate downscaling of temperature and air pressure, nearest-neighbor resampling of other vars
                    for key in data.keys():
                        if (key=='tmean') | (key=='tmin') | (key=='tmax'):
                            data[key] = resize(data[key],mapsize_global,order=1,mode='edge',anti_aliasing=False)
                            #plt.imshow(data[key][750:1300,5000:6000],vmin=-30,vmax=30)
                            #plt.savefig('air_temp_mbar.png',dpi=300)
                            #plt.close()          
                            data[key] = data[key]+temp_delta
                            #plt.imshow(data[key][750:1300,5000:6000],vmin=-30,vmax=30)
                            #plt.savefig('air_temp_down_mbar.png',dpi=300)
                            #plt.close()
                        if key=='pres':
                            data[key] = resize(data[key],mapsize_global,order=1,mode='edge',anti_aliasing=False)
                            #plt.imshow(data[key][750:1300,5000:6000],vmin=500,vmax=1100)
                            #plt.savefig('surface_air_pressure_mbar.png',dpi=300)
                            #plt.close()
                            data[key] = data[key]+pres_delta
                            #plt.imshow(data[key][750:1300,5000:6000],vmin=500,vmax=1100)
                            #plt.savefig('surface_air_pressure_down_mbar.png',dpi=300)
                            #plt.close()
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
                    
                    # Write data to netCDFs
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
                    
                    '''
                    # Compute potential evaporation
                    
                    if vars[vv,0]=='et': # Vegetated
                        albedo,factor = 0.23,1
                    if vars[vv,0]=='ew': # Water
                        albedo,factor = 0.05,0.5
                    if vars[vv,0]=='es': # Soil
                        albedo,factor = 0.15,0.75
                    
                    
                    # Subset data and write to netCDF
                    ncfile.variables['time'][index] = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400
                    ncfile.variables[vars[vv,0]][index,:,:] = pet[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                    '''
                print("Time elapsed is "+str(time.time()-t0)+" sec")
        
            ncfile_pr.close()
            ncfile_ta.close()
            ncfile_et.close()
            ncfile_ew.close()
            ncfile_es.close()

            pdb.set_trace()
            
if __name__ == '__main__':
    main()