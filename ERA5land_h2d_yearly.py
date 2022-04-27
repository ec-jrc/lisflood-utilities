# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:29:27 2021

@author: tilloal
"""
# Import the os module
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
# Change working directory
os.chdir(dir_path)

from nco import Nco
nco = Nco()
import re
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import time
import pandas as pd
from xarray import concat

from tools import *
#%%

## Script 2/4

# script that creates yearly files of daily ERA5-land variables and from aggregated hourly values
# and creates yealry files from monthly files in the specified years.
 
# Requires:
# 1) monthly files of hourly data for the specified variables in the specified years
# (Outputs from "ERA5land_downloader.py")
 
# Output:
# 1) separate netCDF file for chosen daily variable for each year
#%%

#config = load_config(sys.argv[1])
config = load_config("cds_config.cfg")
# Print the current working directory
#print("Current working directory: {0}".format(os.getcwd()))

# select your variable(s); name must be a valid ERA5 CDS API name 
longnames = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature','surface_net_thermal_radiation', 'surface_net_solar_radiation','2m_dewpoint_temperature']

#varnames = ['u10', 'v10','2t','str','ssrd','2d','t2m', 'd2m']
input_string = input('Enter varialbes that need to be processed from the list [u10,v10,2t,str,ssrd,2d,t2m,d2m], enter "all" for all variables')
print("\n")
vin = input_string.split()
if vin == ["all"]:
    varnames = ['u10', 'v10','2t','str','ssrd','2d','t2m', 'd2m','tp']
else:
    varnames=vin
    
# print list
print('list: ', varnames)
# define names of new variables
tvar = ['ws','ta','rn','td','tp','rgd','u10','v10']
var = tvar

file_dir= config['download_folder']
namefile="e5l"

scratchoutdir = config['scratch_folder']
files = glob.glob(os.path.join(config['download_folder'],"hourly",namefile+"_*_*.nc"))    
nfiles= len(files)
files.sort()
yr=0 
yrlist=[]   
data ={}  
#read all the monthly files downloaded from the CDS
for file in files:
#format of namefile must be : xx_xxx_yr_month.nc 
    yrprev=yr
    lf=len(file)
    splitname=os.path.basename(file).split('_')
    yrloc=len(splitname)-2 #the year is always at the end of the filename
    moloc=len(splitname)-1 #the month is always at the end of the filename
    file_yearnc = os.path.basename(file).split('_')[yrloc]
    file_monthnc = os.path.basename(file).split('_')[moloc]
            
    yr = int(file_yearnc.split('.')[0])
    mo = int(file_monthnc.split('.')[0])
    hourly_v = xr.open_dataset(os.path.join(file))
    yrlist.append(yr)
    wc = list(hourly_v.keys())
    
    print("the variables in this file are " + ' | '.join(wc)) 
    common=set(wc).intersection(varnames)
    print("the variables which will be treated are " + ' | '.join(common)) 
    vc=list(common)
    print('month ' + str(mo) + " year " + str(yr))
    tan=len(hourly_v['time'])
    koudur=range(0,tan,24)
    start = time.time() 
          
    if 'u10' in vc and 'v10' in vc:
        v10= hourly_v['v10']
        u10 = hourly_v['u10']
        #wind= wind_uv_to_spd(u10,v10)
        u10=u10.rename('u10')
        daily_u = u10.resample(time='D').mean('time')
        v10=v10.rename('v10')
        daily_v = v10.resample(time='D').mean('time')
        wind= wind_uv_to_spd(u10,v10)
        wind=wind.rename({'ws'})
        daily_w = wind.resample(time='D').mean('time')  
        
        if mo==1:
            data['u10']=daily_u
            data['v10']=daily_v
            data['ws']=daily_w
        else:
            data['u10'] = concat([data['u10'],daily_u],dim='time')
            data['v10'] = concat([data['v10'],daily_v],dim='time')
            data['ws'] = concat([data['ws'],daily_w],dim='time')

    # daily accumulation of precipitation   
    if 'tp' in vc:
        tp= hourly_v['tp'] 
        #convertion to mm
        tp=tp*1000
        daily_pr=tp.resample(time='D').max('time')  
        dptest= daily_pr
        #dptest.plot.surface(yincrease=True)
        daily_pr=daily_pr.rename('tp')
        if mo==1:
            data['tp']=daily_pr
        else:
            data['tp'] = concat([data['tp'],daily_pr],dim='time')
        
      # precipitation: calculate sum with frequency of 24h and multiply by 1000
      # precipitation value is for the day before
      
    # daily mean tempearature
    if '2t' in vc or 't2m'in vc:
        t2m=hourly_v['t2m']
        #convert Kelvin to degrees C
        daily_t2m = t2m.resample(time='D').mean('time')
        daily_t2m=daily_t2m-273.15

        daily_t2m=daily_t2m.rename('ta')
        if mo==1:
            data['ta']=daily_t2m
        else:
            data['ta'] = concat([data['ta'],daily_t2m],dim='time')
            
    # daily means of surface net thermal radiation and surface net solar radiation    
    if 'str'in vc:
        sstr=hourly_v['str']
        
        daily_rn = sstr[koudur,:,:] 
        daily_rn=daily_rn.rename('rn')
        if mo==1:
            data['rn']=daily_rn
        else:
            data['rn'] = concat([data['rn'],daily_rn],dim='time')
            
            
    if 'ssrd'in vc:
        ssrd=hourly_v['ssrd']
        daily_rgd = ssrd[koudur,:,:] 
        daily_rgd= daily_rgd.rename('rgd')
        if mo==1:
            data['rgd']=daily_rgd
        else:
            data['rgd'] = concat([data['rgd'],daily_rgd],dim='time')
            
 
    # daily mean of dew point temperature
    if '2d' in vc or 'd2m' in vc:
        print(vc)
        d2m=hourly_v['d2m']
        #convert Kelvin to degrees C
        d2m=d2m-273.15
        daily_d2m = d2m.resample(time='D').mean('time')
        daily_d2m=daily_d2m.rename({'td'})
        if mo==1:
            data['td']=daily_d2m
        else:
            data['td']= concat([data['td'],daily_d2m],dim='time')

    end=time.time()
    print(end - start)
 #%% 
    if mo==12:
        compression='0'
        file_dates_dly = pd.date_range(start=datetime(yr,1,1), end=datetime(yr+1,1,1)-pd.Timedelta(days=1), freq='D')
        template_lat = np.array(hourly_v['latitude'][:])
        template_lon = np.array(hourly_v['longitude'][:])
        
        if 'v10' in vc:
            vr=var[7]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_v10 = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_v10_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_v10.variables['time'][:] = timef
            ncfile_v10.variables['v10'][:] = data['v10']
            ncfile_v10.close()

        if 'u10' in vc:
            vr=var[6]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_u10 = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_u10_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_u10.variables['time'][:] = timef
            ncfile_u10.variables['u10'][:] = data['u10']
            ncfile_u10.close()

            
        if '2t' in vc or 't2m' in vc:
            vr=var[1]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_ta = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_ta_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_ta.variables['time'][:] = timef
            ncfile_ta.variables['ta'][:] = data['ta']
            ncfile_ta.close()
            
        if '2d' in vc or 'd2m' in vc:
            vr=var[3]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_td = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_td_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_td.variables['time'][:] = timef
            ncfile_td.variables['td'][:] = data['td']
            ncfile_td.close()
            
        if 'ssrd' in vc:
            vr=var[5]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_rgd = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_rgd_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_rgd.variables['time'][:] = timef
            ncfile_rgd.variables['rgd'][:] = data['rgd']
            ncfile_rgd.close()
           
            
        if 'str' in vc:
            vr=var[2]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_rn = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_rn_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_rn.variables['time'][:] = timef
            ncfile_rn.variables['rn'][:] = data['rn']
            ncfile_rn.close()
            
        if 'tp' in vc:
            vr=var[4]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_tp = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_tp_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_tp.variables['time'][:] = timef
            ncfile_tp.variables['tp'][:] = data['tp']
            ncfile_tp.close()
            
        if 'u10' in vc and 'v10' in vc:
            vr=var[0]
            print ('Start generating netcdf file for variable: '+ vr)
            ncfile_ws = initialize_netcdf(os.path.join(scratchoutdir,'e5ld_01deg_ws_' + str(yr) + '.nc'),template_lat,template_lon,vr,'m s-1',compression,1)
            time_value =(file_dates_dly-pd.to_datetime(datetime(1979, 1, 1)))
            ti0=time_value.astype('timedelta64[D]')
            tlist=ti0.tolist()
            timef=np.asarray(list(map(int,tlist)))
            ncfile_ws.variables['time'][:] = timef
            ncfile_ws.variables['ws'][:] = data['ws']
            ncfile_ws.close()
            
            