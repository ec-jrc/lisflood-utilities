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
    return row,col

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat,lon
    
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
df = pd.DataFrame(columns = ["ObsID","StationName","Provider ID","Country code","StationLat","StationLon","Height","Height Units","DrainingArea.km2.Provider","Catchment Area Units","Added Date","River","Catchment","EC_Catchments","Calibration ID","DrainingArea.km2.LDD","LisfloodX","LisfloodY","RealTime","PostProcess","FixedRepPoint","HasWL","HasD","Alive D","StartDate_hist_24","EndDate_hist_24","StartDate_nrt_24","EndDate_nrt_24","StartDate_hist_6","EndDate_hist_6","StartDate_nrt_6","EndDate_nrt_6","DiffDays_hist_24","DiffDays_nrt_24","DiffDays_hist_6","DiffDays_nrt_6","cal_hist_24h","cal_nrt_24h","cal_hist_6h","cal_nrt_6h","CAL_TYPE","Notes","EC_calib","Dam/Lake","EnoughQdata","Val_Start","Val_End","Cal_Start","Cal_End","CatchmentArea","SamplingFrequency"])

# Loop over catchments and snap station locations to the 'correct' grid cell
catchment_dirs = glob.glob(os.path.join(databasedir,'*',""))
catchment_dirs = sorted(catchment_dirs)
for ii in np.arange(len(catchment_dirs)):
    catchment_dir = catchment_dirs[ii]    
    print(str(ii)+" "+os.path.split(os.path.dirname(catchment_dir))[-1])    
    t1 = time.time()
    done = False
   
    # Load discharge data
    try:
        Discharge = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Discharge')
        StationCoordsLat = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lat')[0]
        StationCoordsLon = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lon')[0]
    except:
        print("Unable to load DISCHARGE.mat, skipping")
        continue
        
    # Load catchment boundaries    
    if 'Area' in globals(): del Area
    try:
        Area = np.float(readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/Area'))
        
        '''
        CatchBoundsLat = readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/CatchBounds/Lat')
        CatchBoundsLon = readmatfile(os.path.join(catchment_dir,"BOUNDARIES.mat"),'BOUNDARIES/CatchBounds/Lon')
        
        # Create catchment boundaries array
        if len(CatchBoundsLat[0,:])!=1: CatchBoundsLat = CatchBoundsLat.transpose()
        if len(CatchBoundsLon[0,:])!=1: CatchBoundsLon = CatchBoundsLon.transpose()    
        CatchBounds = np.concatenate((CatchBoundsLon,CatchBoundsLat),axis=1)

        # Plot for testing purposes
        plt.plot(CatchBoundsLon,CatchBoundsLat,color='k',linewidth=0.5)
        
        # Create grid mesh
        lattop = np.ceil(np.nanmax(CatchBoundsLat)/resolution)*resolution-resolution/2
        lattop_ind = np.int(np.round((90-lattop)/resolution-0.5))
        latbottom = np.floor(np.nanmin(CatchBoundsLat)/resolution)*resolution+resolution/2
        latbottom_ind = np.int(np.round((90-latbottom)/resolution-0.5))+1
        lonright = np.ceil(np.nanmax(CatchBoundsLon)/resolution)*resolution-resolution/2
        lonright_ind = np.int(np.round((lonright+180)/resolution-0.5))+1
        lonleft = np.floor(np.nanmin(CatchBoundsLon)/resolution)*resolution+resolution/2
        lonleft_ind = np.int(np.round((lonleft+180)/resolution-0.5))
        latsub = lat[lattop_ind:latbottom_ind]
        lonsub = lon[lonleft_ind:lonright_ind]
        x, y = np.meshgrid(lonsub, latsub)
        x, y = x.flatten(), y.flatten()
        point = np.vstack((x,y)).T 

        # Make catchment mask
        p = Path(CatchBounds) # Make polygon
        grid = p.contains_point(point)
        mask = grid.reshape(len(latsub),len(lonsub)) 
        '''

    except:
        print("Unable to load BOUNDARIES.mat")
        continue

    # Settings
    search_margin = 1 # Grid cells
    area_permissible_error = 20 # Percent        
    
    # Station row col (initial estimate)
    row,col = latlon2rowcol(StationCoordsLat,StationCoordsLon,res,lat[0]+res/2,lon[0]-res/2)
    
    # If we have area estimate
    if 'Area' in globals(): 
        error_map = 100*np.abs(upstreamarea_np/Area-1)
        mask = np.zeros(ldd_np.shape,dtype=np.bool)
        mask[row-search_margin,row+search_margin+1,col-search_margin,col+search_margin+1] = True
        error_map[~mask] = 9999
        if error_map[row,col]<area_permissible_error: 
            done = True
        else:
            error_map
            np.where(error_map==np.min(error_map))
            pdb.set_trace()
            
            
        
        
    if done==False:
    
    
    
    upstreamarea_pcr
    
    
    
    for rr in np.arange(row-search_margin,row+search_margin+1):
        for cc in np.arange(col-search_margin,col+search_margin+1):
            
    point_np = np.zeros(ldd_np.shape,dtype=np.bool)
    point_np[row,col] = True
    point_pcr = pcr.numpy2pcr(pcr.Boolean,point_np,mv=-9999)
    catch_pcr = pcr.catchment(ldd_pcr, point_pcr)
    
    pcr.aguila(catch_pcr)
    

    pdb.set_trace()


pdb.set_trace()
pcr.setclone(ldd)



pcr.catchment(ldd,point)
pdb.set_trace()


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
    StationLat = ADHI_stations['Latitude'][0]
    StationLon = ADHI_stations['Longitude'][0]
    df = df.append({'ObsID':ObsID, 'StationName':StationName, 'StationLat':StationLat,'StationLon':StationLon},ignore_index=True)
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