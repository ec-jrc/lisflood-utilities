#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from tools import *
from skimage.transform import resize

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
    
    # Loop through scenarios and models
    files = glob.glob(os.path.join(config['meteo_data_folder'],'*.nc'))
    scenarios = np.unique([os.path.basename(x).split('_')[3] for x in files])
    for scenario in scenarios:
        files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+scenario+'*.nc'))
        sort(files)
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
            
            # Initialize output netCDF
            varname = 'pr'
            varname_isimip3b = 'pr_' # 'tas_'
            factor = 86400
            offset = 0
            units = 'mm/d'
            if os.path.isdir(os.path.join(config['output_folder'],scenario,model))==False:
                os.makedirs(os.path.join(config['output_folder'],scenario,model))
            outfile = os.path.join(config['output_folder'],scenario,model,varname+'.nc')
            precision = 1
            if os.path.isfile(outfile):
                os.remove(outfile)                
            ncfile = Dataset(outfile, 'w', format='NETCDF4')
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
            ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,450,450,), fill_value=-9999,least_significant_digit=precision)
            ncfile.variables[varname].units = units

            # Loop over input files
            files = glob.glob(os.path.join(config['meteo_data_folder'],'*'+model+'*'+scenario+'*'+varname_isimip3b+'*.nc'))
            for file in files:
                file_year_start = int(os.path.basename(file).split('_')[7])
                file_year_end = int(os.path.basename(file).split('_')[8][:-3])
                count = np.sum((out_dates_dly.year>=file_year_start) & (out_dates_dly.year<=file_year_end))
                if count==0:
                    continue
                
                # Loop over days of input file
                print('Processing '+file)
                t0 = time.time()
                file_dates_dly = pd.date_range(start=datetime(file_year_start,1,1), end=datetime(file_year_end+1,1,1)-pd.Timedelta(days=1), freq='D')                
                dset = Dataset(file)
                for ii in np.arange(len(file_dates_dly)):
                    if file_dates_dly[ii] not in out_dates_dly:
                        continue
                    
                    index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]
                    data = np.array(dset.variables['pr'][ii,:,:])
                    data = factor*data+offset
                    data = np.roll(data,360,axis=1)
                    data = resize(data,mapsize_global,order=0,mode='edge',anti_aliasing=False)
                    data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                    
                    # Write data to netCDF
                    ncfile.variables['time'][index] = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400
                    ncfile.variables[varname][index,:,:] = data                    
                
                print("Time elapsed is "+str(time.time()-t0)+" sec")
                   
            ncfile.close()

            pdb.set_trace()
        
            
            
            
    # Loop through scenarios
    for scenario in scenarios:
        
        # Create output folder
        if os.path.isdir(os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario))==False:
            os.makedirs(os.path.join(config['output_folder'],'step1_Chen_2020_urban',scenario))
            
        # Loop through years
        for year in years:
            print('Processing '+scenario+' '+str(year))
            t0 = time.time()
            src = rasterio.open(os.path.join(config['chen_2020_urban_folder'],scenario.upper(),'global_'+scenario.upper()+'_'+str(year)+'_reprojected.tif'))
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