#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "January 2022"

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
    print('===============================================================================')
    
    # Create output folder
    if os.path.isdir(os.path.join(config['output_folder'],'step1_population_density'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step1_population_density'))
       
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

    # Compute area for each grid-cell (includes correction because lat/lon 
    # values in templates are often not rounded...
    xi, yi = np.meshgrid(np.arange(-180+template_res/2,180+template_res/2,template_res), np.arange(90-template_res/2,-90-template_res/2,-template_res))
    if yi.shape[0]>np.round(180/template_res):
        yi, xi = yi[:-1,:], xi[:-1,:]
    if yi.shape[1]>np.round(360/template_res):
        yi, xi = yi[:,:-1], xi[:,:-1]
    area_map = (40075*template_res/360)**2*np.cos(np.deg2rad(yi))
    
    # List of years with population data
    pop_folders = sorted(glob.glob(os.path.join(config['ghsl_folder'],'*')))
    pop_years = np.array([int(os.path.basename(pop_folder)[9:13]) for pop_folder in pop_folders])

    # List of years
    years = np.arange(config['year_start'],config['year_end']+1).astype(int)


    ############################################################################
    #   Load GHSL population data and upscale to template map resolution
    ############################################################################

    print('-------------------------------------------------------------------------------')
    t0 = time.time()

    pop_upscaled = np.zeros((mapsize_global[0],mapsize_global[1],len(pop_years)),dtype=np.single)    
    for ii in np.arange(len(pop_years)):
        print('Loading and upscaling '+str(pop_years[ii])+' population data')

        # Load raw data (total population per grid-cell)
        # Factor 16 to account for resampling from 0.0025 to 0.01 degrees
        src = rasterio.open(os.path.join(config['ghsl_folder'],'GHS_POP_E'+str(pop_years[ii])+'_GLOBE_R2019A_4326_9ss_V1_0','GHS_POP_E'+str(pop_years[ii])+'_GLOBE_R2019A_4326_9ss_V1_0_reprojected.tif'))
        pop_data = src.read(1).astype(np.single)*16
        pop_data[pop_data<0] = 0
        refmatrix = src.get_transform()
        src.close()
            
        # Make fully global
        add_top = int(np.round((90-refmatrix[3])/-refmatrix[5]))
        add_bottom = int(180/-refmatrix[5]-pop_data.shape[0]-add_top)
        pop_shape = (int(pop_data.shape[1]/2),pop_data.shape[1])
        pop_raw = np.zeros(pop_shape,dtype=np.single)
        pop_raw[add_top:pop_shape[0]-add_bottom,:] = pop_data    
        
        # Area map
        res = 180/pop_raw.shape[0]
        _, yi = np.meshgrid(np.arange(-180+res/2,180+res/2,res), np.arange(90-res/2,-90-res/2,-res))
        yi = yi.astype(np.single)
        ghsl_area_map = (40075*res/360)**2*np.cos(np.deg2rad(yi))
        
        # Convert to density and upscale
        pop_raw = pop_raw/ghsl_area_map
        pop_map = imresize_mean(pop_raw,mapsize_global)
        pop_upscaled[:,:,ii] = pop_map
        
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    

    ############################################################################
    #   Temporally inter- and extrapolate population data
    ############################################################################

    print('-------------------------------------------------------------------------------')
    t0 = time.time()

    pop_big = np.zeros((mapsize_global[0],mapsize_global[1],len(years)),dtype=np.single)    
    for ii in np.arange(len(years)):
        print('Calculating '+str(years[ii])+' population')
        
        # Extrapolate into past
        if years[ii]<np.min(pop_years): 
            annual_change = (pop_upscaled[:,:,1]-pop_upscaled[:,:,0])/(pop_years[1]-pop_years[0])
            number_of_years = years[ii]-pop_years[0]
            pop_big[:,:,ii] = pop_upscaled[:,:,-1]+annual_change*number_of_years
            
        # Extrapolate into future
        elif years[ii]>np.max(pop_years): 
            annual_change = (pop_upscaled[:,:,-1]-pop_upscaled[:,:,-2])/(pop_years[-1]-pop_years[-2])
            number_of_years = years[ii]-pop_years[-1]
            pop_big[:,:,ii] = pop_upscaled[:,:,-1]+annual_change*number_of_years
            
        # Year is included in GHSL, no inter- or extrapolation necessary  
        elif years[ii] in pop_years:
            index = np.argsort(np.abs(pop_years-years[ii]))[0]
            pop_big[:,:,ii] = pop_upscaled[:,:,index]
            
        # Interpolate
        else:
            index = np.argwhere(pop_years<years[ii])[-1][0]
            annual_change = (pop_upscaled[:,:,index+1]-pop_upscaled[:,:,index])/(pop_years[index+1]-pop_years[index])
            number_of_years = years[ii]-pop_years[index]
            pop_big[:,:,ii] = pop_upscaled[:,:,index]+annual_change*number_of_years

    pop_big[pop_big<0] = 0

    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Save population data
    ############################################################################

    print('-------------------------------------------------------------------------------')
    t0 = time.time()
        
    for ii in np.arange(len(years)):
        print('Saving '+str(years[ii])+' population to npz')
        data = pop_big[:,:,ii]
        np.savez_compressed(os.path.join(config['output_folder'],'step1_population_density',str(years[ii])+'.npz'),data=data)

    print("Time elapsed is "+str(time.time()-t0)+" sec")

    
    ############################################################################
    #   Convert to netCDF
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Saving as netCDF')
    t0 = time.time()
    
    varname = 'pop'
    file = os.path.join(config['output_folder'],'step1_population_density',varname+'.nc')

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

    ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,200,200,), fill_value=-9999, least_significant_digit=2)
    ncfile.variables[varname].units = 'total population per grid-cell'

    for year in years:
        print('Saving '+str(year)+' population to netCDF')
        for month in np.arange(1,13):
            
            data = np.load(os.path.join(config['output_folder'],'step1_population_density',str(year)+'.npz'))['data']
            data = np.round(data/10)*10
            data[np.isnan(data)] = 0            
            data = data*area_map # Convert from population density per km2 to total population per grid-cell
            data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)] # Subset to template map area        
            
            index = (year-config['year_start'])*12+month-1
             
            ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
            ncfile.variables[varname][index,:,:] = data
            
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    ncfile.close()
    
if __name__ == '__main__':
    main()