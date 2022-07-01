
"""
Created on Mon Apr  4 17:09:58 2022

@author: Alois Tilloy and Hylke Beck
"""

#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Alois Tilloy"
__date__ = "June 2022"


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
from skimage.transform import rescale
from skimage.io import imread
import matplotlib.pyplot as plt
#import rasterio
import shutil
from datetime import timedelta

# Load configuration file 
config = load_config(sys.argv[1])
namefiles=config['namefiles']
petc = config['petc']
compression = config['compression']

eu_area = [72.25, -25.25, 22.25, 50.25,]
input_res= config['input_res']
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
    dset.set_auto_maskandscale(False)

    template_lat = np.asarray(dset.variables['lat'])
    template_lon = np.asarray(dset.variables['lon'])
    template_res = round(template_lon[1]-template_lon[0],13)
    varname = list(dset.variables.keys())[-1]
    template_np = np.array(dset.variables[varname][:])
    condition = template_np==0

    # Determine map sizes
    mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
   
    # this is the map size for pan-european hydrological anaysis
    scaleR=(0.1/template_res)
    
    # The grid need to be shifted if border has 2 significan digits
    shift=round(eu_area[0]-round(eu_area[0],1),2)
    mapsize_ereu =((np.round((eu_area[0]-eu_area[2])/input_res)+1).astype(int),(np.round((eu_area[3]-eu_area[1])/input_res)+1).astype(int))
    size_eulon = round(mapsize_ereu[0]*scaleR)
    size_eulat = round(mapsize_ereu[1]*scaleR)
    mapsize_europe = (size_eulon,size_eulat)
    mapsize_europe2 = (np.round((eu_area[0]-eu_area[2])/template_res).astype(int),np.round((eu_area[3]-eu_area[1]+1)/template_res).astype(int))
    row_ue,col_lue = latlon2rowcol(template_lat[0],template_lon[0],template_res,eu_area[0]+shift,eu_area[1]-shift)

    mapsize_template = template_np.shape

    # locate the european and template domain in the global domain
    row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)
    row_uppeu,col_lefeu = latlon2rowcol(eu_area[0]+shift,eu_area[1]-shift,template_res,90,-180)

    # Load elevation data, append zeros to top and bottom to make global, and resample to template resolution
    elev = np.zeros((21600,43200),dtype=np.single)

    srcz=imread(os.path.join(config['dem_folder'],'elevation_1KMmn_GMTEDmn.tif'), plugin="tifffile")
    
    # Removed rastrio because it is conflicting with xarray
    #src = rasterio.open(os.path.join(config['dem_folder'],'elevation_1KMmn_GMTEDmn.tif'))

    elev[720:17520,:] = srcz

    elev_global = imresize_mean(elev,mapsize_global)
    elev_europe = elev_global[row_uppeu:row_uppeu+mapsize_europe[0],col_lefeu:col_lefeu+mapsize_europe[1]]
    elev_template = elev_global[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]

    # Prepare temperature downscaling
    tmp = imresize_mean(elev_global,(1800,3600)) # Resample to dimensions of input data - global
    
    #Rescale to output resolution
    tmp = rescale(tmp,scaleR,order=1,mode='edge',anti_aliasing=False)

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
    lat_europe = lat_global[row_uppeu:row_uppeu+mapsize_europe[0],col_lefeu:col_lefeu++mapsize_europe[1]]
   
    # Check if output already exists in scratch folder or output folder
    scratchoutdir = os.path.join(config['scratch_folder'],'e5land_reanalysis')
    finaloutdir = os.path.join(config['output_folder'],'e5land_reanalysis')
        
    # Check if all input variables are present
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
                    #3 possible methods
                    data[key] = resize(data[key],mapsize_europe,order=0,mode='edge',anti_aliasing=False) 
                    
                    #If this method is used no need to subset data to template region
                    #data[key] = data[key].interp_like(obj, method='nearest')

            # Subset data to template region
            for key in data.keys():
                if cover=="global":
                    data[key] = data[key][row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
                if cover=="europe":
                    data[key] = data[key][row_ue:row_ue+len(template_lat),col_lue:col_lue+len(template_lon)]

        # fill the values with NA where condition is false, where there should be NaN

                data[key]= np.where(condition,np.nan,data[key])
                
            #compression part
                if compression=="1":
                    scale_factor=meteo_vars_config[key][KEY_SCALE_FACTOR]   
                    add_offset=meteo_vars_config[key][KEY_OFFSET] 
                    data[key][np.isnan(data[key])] = (-9999 - add_offset) * scale_factor

         
            # Write data to output netCDFs
            time_value = (file_dates_dly[ii]-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400                 
            index = np.where(out_dates_dly==file_dates_dly[ii])[0][0]

            #Potential evapotranspiration is not implemented yet for the set of input from ERA5-land
            #daily variables need to be the accumulation of the previous day
            #no need to move date as the variable is already an accumulation of the previous day
            
            ncfile_pr.variables['time'][index] = time_value
            ncfile_pr.variables['tp'][index,:,:] = data['tp']
            
            #+1 as LISVAP and LISFLOOD use the value of the previous day
            ncfile_ta.variables['time'][index] = time_value+1
            ncfile_ta.variables['ta'][index,:,:] = data['ta']
      
            ncfile_ws.variables['time'][index] = time_value+1 
            ncfile_ws.variables['ws'][index,:,:] = data['ws']
            ncfile_rn.variables['time'][index] = time_value
            ncfile_rn.variables['rn'][index,:,:] = data['rn']
            ncfile_rgd.variables['time'][index] = time_value
            ncfile_rgd.variables['rgd'][index,:,:] = data['rgd']
            ncfile_td.variables['time'][index] = time_value+1
            ncfile_td.variables['td'][index,:,:] = data['td']
                    
            # Generate figures to verify output
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