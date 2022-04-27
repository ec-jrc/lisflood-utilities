# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:09:58 2022

@author: tilloal
"""

#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "March 2022"


import rioxarray # for the extension to load
import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
import xarray as xr


dir_path = os.path.dirname(os.path.realpath(__file__))
# Change working directory
os.chdir(dir_path)

from config_comp import *
import netCDF4 as nc
from tools import *
from skimage.transform import resize
from skimage.io import imread
import matplotlib.pyplot as plt
#import rasterio
import shutil
from datetime import timedelta

# Load configuration file 
config = load_config(sys.argv[1])
#config = load_config("D:/tilloal/Documents/Lisflood_Meteo/config.cfg")
namefiles=config['namefiles']
petc = config['petc']
compression = config['compression']

eu_area = [72.25, -25.25, 22.25, 50.25,]
e5land_res= config['e5land_res']
cover = config['cover']
start=  int(sys.argv[2])
end = int(sys.argv[3])

#start=1987
#end=1988
print("working on " + cover)

if compression=='1':
    print('netcdf files will be compressed using add_offset and scale_factor')

def main():

    # Output dates
    year_start,year_end = start, end
    out_dates_dly = pd.date_range(start=datetime(year_start,1,1), end=datetime(year_end+1,1,1)-pd.Timedelta(days=1), freq='D')
    
    # Load template map
    dset = nc.Dataset(config['templatemap_path'])
    template_lat = np.array(dset.variables['lat'][:])
    template_lon = np.array(dset.variables['lon'][:])
    template_res = template_lon[1]-template_lon[0]
    varname = list(dset.variables.keys())[-1]
    template_np = np.array(dset.variables[varname][:])
    
    #load of template domain with xarray to set a mask on the new data
    Source = xr.open_dataset(config['templatemap_path'])
    obj=Source['area']
    #land mask at template resolution
    condition = obj.isnull()

    # Determine map sizes

    mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
    
    # this is the map size for pan-european hydrological anaysis
    mapsize_europe = (np.round((eu_area[0]-eu_area[2])/template_res).astype(int),np.round((eu_area[3]-eu_area[1])/template_res).astype(int))
    row_ue,col_lue = latlon2rowcol(template_lat[0],template_lon[0],template_res,eu_area[0],eu_area[1])

    mapsize_template = template_np.shape

    # locate the european and template domain in the global domain
    row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)
    row_uppeu,col_lefeu = latlon2rowcol(eu_area[0],eu_area[1],template_res,90,-180)

    # Load elevation data, append zeros to top and bottom to make global, and resample to template resolution
    elev = np.zeros((21600,43200),dtype=np.single)

    srcz=imread(os.path.join(config['dem_folder'],'elevation_1KMmn_GMTEDmn.tif'), plugin="tifffile")
    #src = rasterio.open(os.path.join(config['dem_folder'],'elevation_1KMmn_GMTEDmn.tif'))

    elev[720:17520,:] = srcz
    #elev[720:17520,:] = src.read(1)
    #src.close()

    elev_global = imresize_mean(elev,mapsize_global)
    elev_europe = elev_global[row_uppeu:row_uppeu+mapsize_europe[0],col_lefeu:col_lefeu+mapsize_europe[1]]
    elev_template = elev_global[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
    

    # Prepare temperature downscaling
    tmp = imresize_mean(elev_global,(1800,3600)) # Resample to dimensions of input data - global
    
    tmp = resize(tmp,mapsize_global,order=1,mode='edge',anti_aliasing=False)
 
    #Add option: "Europe" vs "Global"
    if cover=="global":
        elev_delta = elev_global-tmp
    if cover=="europe":
        tmp=tmp[row_uppeu:row_uppeu+mapsize_europe[0],col_lefeu:col_lefeu+mapsize_europe[1]]
        elev_delta = elev_europe-tmp
    #plt.imshow(temp_delta,vmin=-0.5,vmax=0.5)
    #plt.show()
    temp_delta = -6.5*elev_delta/1000 # Simple 6.5 degrees C/km lapse rate
    lat_global = np.repeat(np.resize(np.arange(90-template_res/2,-90-template_res/2,-template_res),(10800,1)),mapsize_global[1],axis=1)
    lat_template = lat_global[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
    # Check if output already exists in scratch folder or output folder
    scratchoutdir = os.path.join(config['scratch_folder'],'e5land_reanalysis')
    finaloutdir = os.path.join(config['output_folder'],'e5land_reanalysis')
    #if (config['delete_existing']==False) & ((os.path.isfile(os.path.join(scratchoutdir,'ta.nc'))==True) | (os.path.isfile(os.path.join(finaloutdir,'ta.nc'))==True)):
     #   print('Already processed, skipping this scenario')
      #  continue
        
    # Check if all input variables are present
    #varnames = ['tas','tasmin','tasmax','hurs','sfcwind','ps','rsds','rlds','pr']
    varnames = ['ws', 'ta','td','rn','rgd']
    nfiles = {}
    for vr in varnames:
        files = glob.glob(os.path.join(config['e5land_folder'], namefiles + "_ta_*.nc"))  
        nfiles[vr] = len(files)
    nfiles_arr = np.array(list(nfiles.values()))
    if (np.max(nfiles_arr)==0) | (any(nfiles_arr<np.max(nfiles_arr))):
        print(nfiles)
        raise Exception('Not all input variables present')
        sys.exit()

    # Initialize output files
    if os.path.isdir(scratchoutdir)==False:
        os.makedirs(scratchoutdir)
    ncfile_ta = initialize_netcdf(os.path.join(scratchoutdir,'ta.nc'),template_lat,template_lon,'ta','degree_Celsius',compression,1)
    ncfile_pr = initialize_netcdf(os.path.join(scratchoutdir,'tp.nc'),template_lat,template_lon,'tp','mm d-1',compression,1)
    
    if petc=="1":
        ncfile_et = initialize_netcdf(os.path.join(scratchoutdir,'et.nc'),template_lat,template_lon,'et','mm d-1',compression,1)
        ncfile_ew = initialize_netcdf(os.path.join(scratchoutdir,'ew.nc'),template_lat,template_lon,'ew','mm d-1',compression,1)
        ncfile_es = initialize_netcdf(os.path.join(scratchoutdir,'es.nc'),template_lat,template_lon,'es','mm d-1',compression,1)
    else: 
        ncfile_rn = initialize_netcdf(os.path.join(scratchoutdir,'rn.nc'),template_lat,template_lon,'rn','J m-2 d',compression,1)
        ncfile_ws = initialize_netcdf(os.path.join(scratchoutdir,'ws.nc'),template_lat,template_lon,'ws','m s-1',compression,1)
        ncfile_rgd = initialize_netcdf(os.path.join(scratchoutdir,'rgd.nc'),template_lat,template_lon,'rgd','J m-2 d',compression,1)
        ncfile_td = initialize_netcdf(os.path.join(scratchoutdir,'td.nc'),template_lat,template_lon,'td','degree_Celsius',compression,1)

    # Loop over input files (MFDataset doesn't work properly)
    print(os.path.join(config['e5land_folder'],namefiles + "_ta_*.nc"))
    files = glob.glob(os.path.join(config['e5land_folder'], namefiles + "_ta_*.nc"))
    

    for file in files:

        splitname=os.path.basename(file).split('_')
        yrloc=len(splitname)-1 #the year is always at the end of the filename
        file_yearnc = os.path.basename(file).split('_')[yrloc]
        file_year = int(file_yearnc.split('.')[0])
        #file_year_end = int(os.path.basename(file).split('_')[6][:-3])
        file_dates_dly = pd.date_range(start=datetime(file_year,1,1), end=datetime(file_year+1,1,1)-pd.Timedelta(days=1), freq='D')
        hits = np.sum((out_dates_dly.year==file_year))
        if hits==0:
          continue

        # Open input files
        print('Processing '+os.path.basename(file))
        t0 = time.time()
        dset_tmean = xr.open_dataset(file,diskless=True) # degrees C

        #file2 = glob.glob(os.path.join('Z:\ClimateRun3\ERA5-land',"ta","0.1_deg", "e5ldxc_ta_1981.nc"))[0]
        #compare_dset= xr.open_dataset(file2,diskless=True)

        #ta_vals=dset_tmean['ta']
        #ta_vals=ta_vals.rio.write_crs(4326)
        #ta_vals=ta_vals.rio.write_nodata('nan')
        #ta_vals=ta_vals.rio.set_spatial_dims('lon','lat', inplace=True)
        dset_tdew = xr.open_dataset(file.replace('ta','td'),diskless=True) 
        dset_wind = xr.open_dataset(file.replace('ta','ws'),diskless=True) # m/s
        dset_tdew = xr.open_dataset(file.replace('ta','td'),diskless=True) # Pa
        dset_swd = xr.open_dataset(file.replace('ta','rgd'),diskless=True) # W/m2
        dset_lwd = xr.open_dataset(file.replace('ta','rn'),diskless=True) # W/m2
        dset_pr = xr.open_dataset(file.replace('ta','tp'),diskless=True) # mm/d
        
        # Loop over days of input file
        for ii in np.arange(len(file_dates_dly)):

            if file_dates_dly[ii] not in out_dates_dly:
                continue
   
            print('Processing year '+ str(file_year) +', date: '+str(file_dates_dly[ii])+', time stamp: '+datetime.now().strftime("%d/%m/%Y, %H:%M:%S")+')')                    

            # Read data from input files
            data = {}
            #index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]                            
            data['ta'] = dset_tmean['ta'][ii,:,:] # degrees C
            if petc=="1":
               data['ws'] = dset_wind.variables['ws'][ii,:,:]*0.75 # m/s (factor 0.75 to translate from 10-m to 2-m height)
            else:
                data['ws'] = dset_wind['ws'][ii,:,:] # m/s (the factor will be applied in LISVAP)
                data['rgd'] = dset_swd['rgd'][ii,:,:] # J/m2/d
                data['rn'] = dset_lwd['rn'][ii,:,:] # J/m2/d
                data['tp'] = dset_pr['tp'][ii,:,:] # mm/d
                data['td'] = dset_tdew['td'][ii,:,:]
                
     
            
            # Switch eastern and western hemispheres
            #for key in data.keys():
            #    data[key] = np.roll(data[key],int(data[key].shape[1]/2),axis=1)
        
            # Simple lapse rate downscaling of temperature and air pressure, nearest-neighbor resampling of other vars

            for key in data.keys():
                data[key]=data[key].rio.write_crs(4326)
                data[key]=data[key].rio.write_nodata('nan')
                data[key]=data[key].rio.set_spatial_dims('lon','lat')
                data[key]=data[key].rio.interpolate_na(method='nearest')
                if (key=='ta')| (key=='td'):
                    data[key] = resize(data[key],mapsize_europe,order=1,mode='constant',anti_aliasing=False)
                    data[key] = np.around(data[key]+temp_delta,1)
                else:
                    data[key] = resize(data[key],mapsize_europe,order=0,mode='edge',anti_aliasing=False)    
                    
            # Subset data to template region
            for key in data.keys():
                if cover=="global":
                    data[key] = data[key][row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                if cover=="europe":
                    data[key] = data[key][row_ue:row_ue+len(template_lat),col_lue:col_lue+len(template_lon)]

        # fill the values with NA where condition is false, where there should be NaN
                data[key]= np.where(condition,np.nan,data[key])
                if compression=="1":
                    scale_factor=meteo_vars_config[key][KEY_SCALE_FACTOR]   
                    add_offset=meteo_vars_config[key][KEY_OFFSET] 
                    data[key][np.isnan(data[key])] = (-9999 - add_offset) * scale_factor

         
            #Potential evapotranspiration is not implemented yet for the set of input from ERA5-land
            if petc=="1":
            # Compute potential evaporation
                albedo = {'et':0.23,'ew':0.05,'es':0.15}
                factor = {'et':1,'ew':0.5,'es':0.75}
                doy = int(datetime.strftime(file_dates_dly[ii],'%j'))
                pet = potential_evaporation(data,albedo,factor,doy,lat_template,elev_template)
            
            # Write data to output netCDFs
            time_value = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400                 
            index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]
            
            #daily variables need to be the accumulation of the previous day
            #no need to move date as the variable is already an accumulation of the previous day
            
            ncfile_pr.variables['time'][index] = time_value
            ncfile_pr.variables['tp'][index,:,:] = data['tp']
            
            #+1 as LISVAP and LISFLOOD use the value of the previous day
            ncfile_ta.variables['time'][index] = time_value+1
            ncfile_ta.variables['ta'][index,:,:] = data['ta']
      
            
            if petc=="1": 
                ncfile_et.variables['time'][index] = time_value
                ncfile_et.variables['et'][index,:,:] = pet['et']
                ncfile_ew.variables['time'][index] = time_value
                ncfile_ew.variables['ew'][index,:,:] = pet['ew']
                ncfile_es.variables['time'][index] = time_value
                ncfile_es.variables['es'][index,:,:] = pet['es']
            else:
                ncfile_ws.variables['time'][index] = time_value+1 
                ncfile_ws.variables['ws'][index,:,:] = data['ws']
                ncfile_rn.variables['time'][index] = time_value
                ncfile_rn.variables['rn'][index,:,:] = data['rn']
                ncfile_rgd.variables['time'][index] = time_value
                ncfile_rgd.variables['rgd'][index,:,:] = data['rgd']
                ncfile_td.variables['time'][index] = time_value+1
                ncfile_td.variables['td'][index,:,:] = data['td']
                    
            # Generate figures to verify output
            #plt.imshow(data['ta'][1600:1700,1400:1500])
            #plt.show()

            if ii==0:
                makefig('figures','ta',data['ta'],0,12)
                #makefig('figures','pr',data['pr'],0,12)
                for key in data.keys():
                    makefig('figures',key,data[key],np.min(data[key]),np.max(data[key]))
                makefig('figures','elev_template',elev_template,0,6000)
        
        # Close input files
        dset_tmean.close()
        dset_wind.close()
        dset_swd.close()
        dset_lwd.close()
        dset_pr.close()
        dset_tdew.close()
        

        print("Time elapsed is "+str(time.time()-t0)+" sec") 
    # Close output files
    ncfile_pr.close()
    ncfile_ta.close()
    
    if "petc"==1: 
        ncfile_et.close()
        ncfile_ew.close()
        ncfile_es.close()
    else:  
        ncfile_ws.close()
        ncfile_rn.close()
        ncfile_rgd.close()
        ncfile_td.close()
    # Move output from scratch folder to output folder
    print('-------------------------------------------------------------------------------')
    if os.path.isdir(finaloutdir)==False:
        os.makedirs(finaloutdir)
    for file in glob.glob(os.path.join(scratchoutdir,'*')):
        t0 = time.time()
        print('Moving '+os.path.basename(file)+' ('+str(round(os.path.getsize(file)/10**9))+' GB) to '+finaloutdir)
        shutil.copy(file, finaloutdir)
        print("Time elapsed is "+str(time.time()-t0)+" sec")
    shutil.rmtree(scratchoutdir)

if __name__ == '__main__':
    main()