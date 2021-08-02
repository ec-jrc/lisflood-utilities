#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "July 2021"

import os, sys, glob, time, h5py
import scipy.io
import pdb
import pandas as pd
import numpy as np
import pcraster as pcr
from netCDF4 import Dataset
import matplotlib.pyplot as plt

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

# Load LDD file
dset = Dataset(lddpath)
ldd_np = np.array(dset.variables['ldd'][:])
ldd_np[np.isnan(ldd_np)] = -9999 
ldd_np = ldd_np.astype(np.int16)
lat = np.array(dset.variables['lat'][:])
lon = np.array(dset.variables['lon'][:])
res = np.diff(lon)[0]
lat_upper = lat[0]+res/2
lon_left = lon[0]-res/2

# ##### LOAD UPS MAP
dset = Dataset(r"//zeus/hylkeb/Temp/ups.nc")
ups_np = np.array(dset.variables['ups'][:])

# Set pcraster clone map
pcr.setclone(ldd_np.shape[0],ldd_np.shape[1],res,lon[0]-res/2,lat[0]-res/2)

# Gridcell area map
xi, yi = np.meshgrid(lon, lat)
area_np = (40075*res/360)**2*np.cos(np.deg2rad(yi))
area_pcr = pcr.numpy2pcr(pcr.Scalar,area_np,mv=-9999)

# Upstream area map
ldd_pcr = pcr.numpy2pcr(pcr.Ldd,ldd_np,mv=-9999)
upstreamarea_pcr = pcr.accuflux(ldd_pcr,area_pcr)
upstreamarea_np = pcr.pcr2numpy(upstreamarea_pcr,mv=-9999)

# Initialize qgis csv
df = pd.DataFrame(columns = ["ObsID","StationName","Provider ID","Country code","StataaaaLat","StataaaaLon","Height","Height Units","DrainingArea.km2.Provider","Catchment Area Units","Added Date","River","Catchment","EC_Catchments","Calibration ID","DrainingArea.km2.LDD","LisfloodX","LisfloodY","RealTime","PostProcess","FixedRepPoint","HasWL","HasD","Alive D","StartDate_hist_24","EndDate_hist_24","StartDate_nrt_24","EndDate_nrt_24","StartDate_hist_6","EndDate_hist_6","StartDate_nrt_6","EndDate_nrt_6","DiffDays_hist_24","DiffDays_nrt_24","DiffDays_hist_6","DiffDays_nrt_6","cal_hist_24h","cal_nrt_24h","cal_hist_6h","cal_nrt_6h","CAL_TYPE","Notes","EC_calib","Dam/Lake","EnoughQdata","Val_Start","Val_End","Cal_Start","Cal_End","CatchmentArea","SamplingFrequency"])

# Loop over catchments and snap station locations to the 'correct' grid cell
catchment_dirs = glob.glob(os.path.join(databasedir,'*',""))
catchment_dirs = sorted(catchment_dirs)
for ii in np.arange(len(catchment_dirs)):
    if ii!=421: continue
    catchment_dir = catchment_dirs[ii]  
    print("===============================================================================")
    print(str(ii)+" "+os.path.split(os.path.dirname(catchment_dir))[-1])    
    t1 = time.time()
    done = False
    toosmall = False
    snapped = False
    
    if 'Discharge' in globals(): del Discharge
    if 'StatLat' in globals(): del StatLat
    if 'StatLon' in globals(): del StatLon
    if 'Area' in globals(): del Area
    if 'CatchCentroidLat' in globals(): del CatchCentroidLat
    if 'CatchCentroidLon' in globals(): del CatchCentroidLon

    
    ############################################################################
    #   Load discharge and catchment boundary data
    ############################################################################
    
    # Load discharge data
    try:
        Discharge = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Discharge')
        StatLat = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lat')[0][0]
        StatLon = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lon')[0][0]
        StatRow,StatCol = latlon2rowcol(StatLat,StatLon,res,lat_upper,lon_left)
    except:
        print("Unable to load DISCHARGE.mat, skipping")
        continue
        
    # Load catchment boundaries    
    try:
        Area = float(readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/Area'))        
        CatchBoundsLat = readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/CatchBounds/Lat')[0]
        CatchBoundsLon = readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/CatchBounds/Lon')[0]
        CatchCentroidLat = (np.nanmin(CatchBoundsLat)+np.nanmax(CatchBoundsLat))/2
        CatchCentroidLon = (np.nanmin(CatchBoundsLon)+np.nanmax(CatchBoundsLon))/2
    except:
        print("Unable to load BOUNDARIES.mat")
    
        
    ############################################################################
    #   If we have area and catchment boundary estimates, we will automatically 
    #   determine the most suitable grid cell
    ############################################################################
    
    if 'Area' in globals():
    
        if Area<10000:             
            print("Catchment too small")
            continue

        # Snap settings
        permissible_area_error = 20 # Maximum allowed catchment area error (percent)
        permissible_centroid_error = 0.5 # Maximum allowed centroid location error (degrees)
        margin = 2 # Search box size (grid cells around station)
        
        # Compute area error
        area_error_map = 100*np.abs(upstreamarea_np-Area)/Area
        mask = np.zeros(ldd_np.shape,dtype=bool)
        mask[StatRow-margin:StatRow+margin+1,StatCol-margin:StatCol+margin+1] = True
        area_error_map[~mask] = 9999
        
        # Compute centroid error
        centroid_error_map = np.zeros(ldd_np.shape,dtype=np.single)+9999
        indices = zip(*np.where(mask))
        for i, j in indices:
            point_np = np.zeros(ldd_np.shape,dtype=bool)
            point_np[i,j] = True
            point_pcr = pcr.numpy2pcr(pcr.Boolean,point_np,mv=-9999)
            catch_pcr = pcr.catchment(ldd_pcr, point_pcr)        
            catch_np = pcr.pcr2numpy(catch_pcr,mv=-9999)
            catch_np = catch_np==1
            if np.sum(catch_np)==0: continue
            ind = np.where(catch_np)
            centroid_row = (ind[0][0]+ind[0][-1])/2
            centroid_col = (ind[1][0]+ind[1][-1])/2
            centroid_lat,centroid_lon = rowcol2latlon(centroid_row,centroid_col,res,lat_upper,lon_left)
            centroid_error_map[i,j] = np.sqrt((centroid_lat-CatchCentroidLat)**2+(centroid_lon-CatchCentroidLon)**2)                       
        
        # Determine the most appropriate station grid cell
        area_error = area_error_map[StatRow,StatCol]
        centroid_error = centroid_error_map[StatRow,StatCol]
        if (area_error<permissible_area_error) & (centroid_error<permissible_centroid_error):
            StatRowCorr,StatColCorr = StatRow,StatCol
            print('Station snapping not necessary')
        else:
            rank_area_error_map = area_error_map.ravel().argsort().argsort().reshape(area_error_map.shape)
            rank_centroid_error_map = centroid_error_map.ravel().argsort().argsort().reshape(centroid_error_map.shape)
            rank_cum_map = rank_area_error_map+rank_centroid_error_map            
            StatRowCorr,StatColCorr = np.where(rank_cum_map==np.min(rank_cum_map))
            StatRowCorr,StatColCorr = StatRowCorr[0],StatColCorr[0]
            StatLatCorr,StatLonCorr = rowcol2latlon(StatRowCorr,StatColCorr,res,lat_upper,lon_left)
            print('Station snapped')

        area_error = area_error_map[StatRowCorr,StatColCorr]
        centroid_error = centroid_error_map[StatRowCorr,StatColCorr]    
        print("Area error: "+str(area_error)+" %, centroid error: "+str(centroid_error)+" degrees")
        
        # Plot stuff
        if (area_error>permissible_area_error) | (centroid_error>permissible_centroid_error):
            point_np = np.zeros(ldd_np.shape,dtype=bool)
            point_np[StatRowCorr,StatColCorr ] = True
            point_pcr = pcr.numpy2pcr(pcr.Boolean,point_np,mv=-9999)
            catch_pcr = pcr.catchment(ldd_pcr, point_pcr)        
            catch_np = pcr.pcr2numpy(catch_pcr,mv=-9999)
            catch_np = catch_np==1
            plt.clf()
            extent = [lon_left,lon_left+len(lon)*res,lat_upper-len(lat)*res,lat_upper] # (left, right, bottom, top)
            plt.imshow(np.single(catch_np)+np.log10(upstreamarea_np)/4, extent=extent)
            plt.plot(CatchBoundsLon,CatchBoundsLat)
            plt.gca().set_xlim(np.nanmin(CatchBoundsLon)-6,np.nanmax(CatchBoundsLon)+6)
            plt.gca().set_ylim(np.nanmin(CatchBoundsLat)-6,np.nanmax(CatchBoundsLat)+6)
            plt.title(str(ii))
            plt.show()
        
        
    pdb.set_trace()
    
    
    ############################################################################
    #   If we don't have area and catchment boundary estimates, we have to 
    #   pick the most suitable grid cell manually
    ############################################################################
            
    #if 'Area' not in globals():     
    #    continue
    
    
    ############################################################################
    #   If we don't have area and catchment boundary estimates, we have to 
    #   pick the most suitable grid cell manually
    ############################################################################        



'''
# pysheds is broken...
affine = Affine(res,0,lon[0],0,-res,lat[0])
grid = Grid()
grid.add_gridded_data(data=ldd, data_name='dir',
                          affine=affine,
                          crs=grid.crs,
                          nodata=np.NaN)

#         N  NE E  SE S SW  W  NW
dirmap = (8, 9, 6, 3, 2, 1, 4, 7)

grid.accumulation(data='dir', dirmap=dirmap, out_name='acc')

plt.figure(0)
plt.imshow(grid.acc[208:213,579:584])

plt.figure(1)
plt.imshow(ldd[208:213,579:584])

plt.show()


# Specify pour point
x, y = 34.91584983, 11.2406754833

# Delineate the catchment
grid.catchment(data='dir', x=x, y=y, dirmap=dirmap, out_name='catch',
               recursionlimit=15000, xytype='label', nodata_out=0)

grid.clip_to('catch')
grid.view('catch')

plt.figure(1)
plt.imshow(ldd)

plt.figure(2)
plt.imshow(grid.catch)

plt.show()

pdb.set_trace()
'''



'''
# Loop over ADHI stations
ADHI_stations = pd.read_csv(os.path.join(rawdir,'ADHI_restricted_raw_data','ADHI_stations.csv'))
for ii in np.arange(len(ADHI_stations)):

    # Basic station info
    ObsID = str(ii)
    StationName = ADHI_stations['Name'][0]
    StatLat = ADHI_stations['Latitude'][0]
    StatLon = ADHI_stations['Longitude'][0]
    df = df.append({'ObsID':ObsID, 'StationName':StationName, 'StatLat':StatLat,'StatLon':StatLon},ignore_index=True)
    pdb.set_trace()

    # Catchmemt area and boundaries
	clear BOUNDARIES	
	shp_filepath = [rawdir filesep 'ADHI_restricted_raw_data' filesep 'CatchmentBoundaries' filesep 'ADHIcatch_' ADHI_co '.shp'];
	if ~exist(shp_filepath), continue; end
	shp = shaperead(shp_filepath);
	BOUNDARIES.CatchBounds.Lat = single(shp.Y);
	BOUNDARIES.CatchBounds.Lon = single(shp.X);
	BOUNDARIES.Area = shp.Area;
'''