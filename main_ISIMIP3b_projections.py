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

# Load configuration file 
config = load_config(sys.argv[1])

def main():
    
    # Load template map
    dset = netCDF4.Dataset(config['templatemap_path'])
    template_lat = dset.variables['lat'][:]
    template_lon = dset.variables['lon'][:]
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
    elev = imresize_mean(elev,mapsize_global)

    # Prepare temperature lapse rate downscaling assuming lapse rate of 6.5°C/km
    tmp = imresize_mean(elev,(360,720))
    tmp = resize(tmp,mapsize_global,order=1,mode='edge',anti_aliasing=False)
    elev_delta = elev-tmp
    temp_delta = -6.5*elev_delta/1000
    
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
            
            # Loop over variables
            vars = np.array([\
                ['et','tas_',1,0,'mm d-1',1],\
                ['e','tas_',1,0,'mm d-1',1],\
                ['es','tas_',1,0,'mm d-1',1],\
                ['ta','tas_',1,-273.15,'degree_Celsius',1],\
                ['pr','pr_',86400,0,'mm d-1',1],\
                ])
                
            for vv in np.arange(len(vars)):
                
                # Initialize netCDF
                outdir = os.path.join(config['output_folder'],scenario,model)
                if os.path.isdir(outdir)==False:
                    os.makedirs(outdir)
                outfile = os.path.join(outdir,vars[vv,0]+'.nc')
                if os.path.isfile(outfile):
                    os.remove(outfile)                
                ncfile = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
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
                ncfile.createVariable(vars[vv,0], np.single, ('time', 'lat', 'lon'), zlib=True,\
                    chunksizes=(1,450,450,), fill_value=-9999,least_significant_digit=int(vars[vv,5]))
                ncfile.variables[vars[vv,0]].units = vars[vv,4]
                
                
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
                
                if (vars[vv,0]=='et') | (vars[vv,0]=='e') | (vars[vv,0]=='es'):
                
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
                        dset_tmean = netCDF4.Dataset(file) # K
                        dset_tmin = netCDF4.Dataset(file.replace('tas_','tasmin_')) # K
                        dset_tmax = netCDF4.Dataset(file.replace('tas_','tasmax_')) # K
                        dset_sh = netCDF4.Dataset(file.replace('tas_','huss_')) # kg kg-1
                        dset_wind = netCDF4.Dataset(file.replace('tas_','sfcwind_')) # m s-1
                        dset_swd = netCDF4.Dataset(file.replace('tas_','rsds_')) # W m-2
                        dset_lwd = netCDF4.Dataset(file.replace('tas_','rlds_')) # W m-2
                        
                        # Loop over days of input file
                        for ii in np.arange(10): #np.arange(len(file_dates_dly)):
                            if file_dates_dly[ii] not in out_dates_dly:
                                continue
                            
                            # Read data from input files
                            data = {}
                            index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]                            
                            data['tmean'] = np.array(dset_tmean.variables['tas'][ii,:,:],dtype=np.single)-273.15 # degrees C
                            data['tmin'] = np.array(dset_tmin.variables['tasmin'][ii,:,:],dtype=np.single)-273.15 # degrees C
                            data['tmax'] = np.array(dset_tmax.variables['tasmax'][ii,:,:],dtype=np.single)-273.15 # degrees C
                            data['sh'] = np.array(dset_sh.variables['huss'][ii,:,:],dtype=np.single) # kg kg-1
                            data['wind'] = np.array(dset_wind.variables['sfcwind'][ii,:,:],dtype=np.single) # m s-1 #### IS THIS AT 2-M HEIGHT????!?
                            data['swd'] = np.array(dset_swd.variables['rsds'][ii,:,:],dtype=np.single) # W m-2
                            data['lwd'] = np.array(dset_lwd.variables['rlds'][ii,:,:],dtype=np.single) # W m-2
                            for key in data.keys():
                                data[key] = np.roll(data[key],int(data[key].shape[1]/2),axis=1)
                        
                            # Simple lapse rate downscaling of temperature, nearest-neighbor resampling of other vars
                            for key in data.keys():
                                if (key=='tmean') | (key=='tmin') | (key=='tmax'):
                                    data[key] = resize(data[key],mapsize_global,order=1,mode='edge',anti_aliasing=False)
                                    data[key] = data[key]+temp_delta
                                else:
                                    data[key] = resize(data[key],mapsize_global,order=0,mode='edge',anti_aliasing=False)
                            
                            pdb.set_trace()
                            
                            data[varname[vv,0]] = potential_evaporation(data,albedo,cropf)
                            
                            # Subset data and write to netCDF
                            ncfile.variables['time'][index] = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400
                            ncfile.variables[vars[vv,0]][index,:,:] = data[varname[vv,0]][row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                        
                        print("Time elapsed is "+str(time.time()-t0)+" sec")
                    
                    
                ncfile.close()

            pdb.set_trace()
            
if __name__ == '__main__':
    main()