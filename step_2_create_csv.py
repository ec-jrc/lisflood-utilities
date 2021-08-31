#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "July 2021"

import os, sys, glob, time, h5py, scipy.io, pdb
import pandas as pd
import numpy as np
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

# Load upstream map
dset = Dataset(ups_path)
upstreamarea_np = np.array(dset.variables['ups'][:])

# Initialize empty stations and Qtss csvs
stations = pd.DataFrame(columns = ["ObsID","StationName","Provider ID","Country code","StationLat","StationLon","DrainingArea.km2.Provider","River","DrainingArea.km2.LDD"])
#Qtss = pd.DataFrame(columns = column_names)

# Create csv output folder
if os.path.isdir(csv_dir)==False:
    os.mkdir(csv_dir)


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
    t1 = time.time()

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
    StatLatCorr,StatLonCorr = rowcol2latlon(corr_loc.iloc[0],corr_loc.iloc[1],res,lat_upper,lon_left)
    Area_LDD = upstreamarea_np[corr_loc.iloc[0],corr_loc.iloc[1]]
    
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
    dates_ref = pd.date_range(start=datetime(1900,1,1), end=datetime(2049,12,31), freq='D').round('D')
    ind = np.where((dates_ref>=date_start) & (dates_ref<=date_end))
    Q_ts = Discharge[ind]
    if np.sum(np.isnan(Q_ts)==False)/365.25 < min_record_length:
        print("Discharge record too short, skipping")
        continue

    # Load catchment boundaries    
    try:
        Area = float(readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/Area'))        
    except:
        print("Unable to load BOUNDARIES.mat")
        
    
    ############################################################################
    #   Insert data into dataframes
    ############################################################################
    
    #stations.append(ignore_index=True)
    station.loc[count] = [count,Station,ID,"",StatLatCorr,StatLonCorr,Area,"",Area_LDD]
    pdb.set_trace()
    
    
    Q_ts
    
    count = count+1
    
    print('Time elapsed is ' + str(time.time() - t1) + ' sec')