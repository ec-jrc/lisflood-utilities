#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "July 2021"

import os, sys, glob, time, h5py, scipy.io, pdb, csv
import pandas as pd
import numpy as np
import pcraster as pcr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import dateutil.parser
from datetime import datetime, timedelta

def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()

def readmatfile(filepath,var):
    try:
        f = h5py.File(filepath,'r') 
        data = f.get(var)[()]
        data = data.transpose() 
        f.close()
    except:
        try:
            mat = scipy.io.loadmat(filepath)
            data = eval("mat['"+var.replace("/","'][0,0]['")+"'][:]")
        except:
            pass    
    return data

def save_netcdf(file, varname, data, least_sig_dig, lat, lon):

    if os.path.isfile(file)==True: 
        os.remove(file)

    ncfile = Dataset(file, 'w', format='NETCDF4')

    ncfile.createDimension('lon', len(lon))
    ncfile.createDimension('lat', len(lat))

    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'

    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'][:] = lat
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'].long_name = 'latitude'
    
    ncfile.createVariable(varname, data.dtype, ('lat', 'lon'), zlib=True, chunksizes=(32,32,), fill_value=-9999, least_significant_digit=least_sig_dig)

    ncfile.variables[varname][:,:] = data
    
    ncfile.close()
    
t1 = time.time()
    
# Load configuration file
config = pd.read_csv('config.cfg',header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")
    string = string.replace("=","=r")
    exec(string)

date_start = dateutil.parser.parse(date_start)
date_end = dateutil.parser.parse(date_end)
min_record_length = float(min_record_length)

# Date stuff
dates_ref = pd.date_range(start=datetime(1900,1,1), end=datetime(2049,12,31), freq='D').round('D')
ind = np.where((dates_ref>=date_start) & (dates_ref<=date_end))

# Load LDD map
dset = Dataset(ldd_path)
ldd_np = np.array(dset.variables['ldd'][:])
ldd_np[np.isnan(ldd_np)] = -9999 
ldd_np = ldd_np.astype(np.int16)
lat = np.array(dset.variables['lat'][:])
lon = np.array(dset.variables['lon'][:])
res = np.diff(lon)[0]
lat_upper = lat[0]+res/2
lon_left = lon[0]-res/2

# Set pcraster clone map
pcr.setclone(ldd_np.shape[0],ldd_np.shape[1],res,lon[0]-res/2,lat[0]-res/2)

# Convert LDD to pcraster format
ldd_pcr = pcr.numpy2pcr(pcr.Ldd,ldd_np,mv=-9999)

# Load upstream map
dset = Dataset(ups_path)
upstreamarea_np = np.array(dset.variables['ups'][:])

# Initialize empty stations and Qtss csvs
# If you add/remove/shift fields, be sure to change the "stations.loc[count] = [ ...]" line as well
stations = pd.DataFrame(columns = ["ObsID","StationName","Provider ID","Country code","StationLat","StationLon","DrainingArea.km2.Provider","River","DrainingArea.km2.LDD"])
Qtss = pd.DataFrame(index=dates_ref[ind])

# Create csv output folder
if os.path.isdir(csv_dir)==False:
    os.mkdir(csv_dir)
if os.path.isdir(os.path.join(csv_dir,'catch_masks'))==False:
    os.mkdir(os.path.join(csv_dir,'catch_masks'))
    

############################################################################
#   Loop over catchments to automatically snap station locations to the 
#   'correct' grid cell. If snapping doesn't work, a window is opened to 
#   manually select the correct location
############################################################################

catchment_dirs = glob.glob(os.path.join(database_dir,'*',""))
catchment_dirs = sorted(catchment_dirs)
count = 0
for ii in np.arange(len(catchment_dirs)):

    catchment_dir = catchment_dirs[ii]  
    ID = os.path.split(os.path.dirname(catchment_dir))[-1]
    print("===============================================================================")
    print(str(ii)+" "+ID)

    if 'Discharge' in globals(): del Discharge
    if 'StatLat' in globals(): del StatLat
    if 'StatLon' in globals(): del StatLon
    if 'Station' in globals(): del Station
    if 'Area' in globals(): del Area
    if 'corr_loc' in globals(): del corr_loc
    
    
    ############################################################################
    #   Load data
    ############################################################################
    
    # Load corrected station location
    if os.path.isfile(os.path.join(corrected_locations_dir,ID+'.txt'))==False:
        print('Corrected station location not found, skipping')
        continue
    corr_loc = pd.read_csv(os.path.join(corrected_locations_dir,ID+'.txt'),header=None,index_col=False)
    StatRowCorr,StatColCorr = int(corr_loc.iloc[1]),int(corr_loc.iloc[0])
    StatLatCorr,StatLonCorr = rowcol2latlon(StatRowCorr,StatColCorr,res,lat_upper,lon_left)
    Area_LDD = upstreamarea_np[StatRowCorr,StatColCorr]
    
    # Load discharge data
    try:
        Discharge = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Discharge')
        StatLat = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lat')[0][0]
        StatLon = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lon')[0][0]
        Station = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Station')[0]
        Station = ''.join(map(chr,Station))
    except:
        print("Unable to load DISCHARGE.mat, skipping")
        continue
    Q_ts = Discharge[ind]
    if np.sum(np.isnan(Q_ts)==False)/365.25 < min_record_length:
        print("Discharge record too short, skipping")
        continue

    # Load catchment boundaries    
    try:
        Area = float(readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/Area'))        
    except:
        Area = np.NaN
        print("Unable to load BOUNDARIES.mat")
        
    
    ############################################################################
    #   Insert data into dataframes
    ############################################################################

    stations.loc[count] = [count,Station,ID,"",StatLatCorr,StatLonCorr,Area,"",Area_LDD]
    Qtss[count] = Q_ts
    
        
    ############################################################################
    #   Save catchment mask
    ############################################################################
    
    point_np = np.zeros(ldd_np.shape,dtype=bool)
    point_np[StatRowCorr,StatColCorr] = True
    point_pcr = pcr.numpy2pcr(pcr.Boolean,point_np,mv=-9999)
    catch_pcr = pcr.catchment(ldd_pcr, point_pcr)        
    catch_np = pcr.pcr2numpy(catch_pcr,mv=-9999)
    catch_np = (catch_np==1).astype(np.uint8)
    filename = str(count).zfill(5)+'_'+ID+'.nc'
    save_netcdf(os.path.join(csv_dir,'catch_masks',filename), 'mask', catch_np, 3, lat, lon)
    
    count = count+1
    

############################################################################
#   Save dataframes to csv
############################################################################

print('Saving dataframes to csv')
stations.to_csv(os.path.join(csv_dir,'stations.csv'),index=False,quoting=csv.QUOTE_NONNUMERIC)
Qtss.to_csv(os.path.join(csv_dir,'Qtss.csv'))

print('Total time elapsed is ' + str(time.time() - t1) + ' sec')

pdb.set_trace()