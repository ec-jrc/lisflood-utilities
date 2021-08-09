#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "July 2021"

import os, sys, glob, time, h5py, scipy.io, pdb
import pandas as pd
import numpy as np
import pcraster as pcr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import geopandas as gpd

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
    
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    global coords
    coords = np.round([ix, iy]).astype(int)
    print("Selected location "+str(coords))
    plt.close()
    
# Load configuration file
config = pd.read_csv('config.cfg',header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")
    string = string.replace("=","=r")
    exec(string)

# Large rivers shapefile
# http://ihp-wins.unesco.org/layers/geonode:world_rivers
rivers_shp = gpd.read_file(r'\\zeus\hylkeb\RESEARCH\Data\world_rivers\world_rivers.shp')

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

# Set pcraster clone map
pcr.setclone(ldd_np.shape[0],ldd_np.shape[1],res,lon[0]-res/2,lat[0]-res/2)

# Create gridcell area map
xi, yi = np.meshgrid(lon, lat)
area_np = (40075*res/360)**2*np.cos(np.deg2rad(yi))
area_pcr = pcr.numpy2pcr(pcr.Scalar,area_np,mv=-9999)

# Create upstream area map
ldd_pcr = pcr.numpy2pcr(pcr.Ldd,ldd_np,mv=-9999)
upstreamarea_pcr = pcr.accuflux(ldd_pcr,area_pcr)
upstreamarea_np = pcr.pcr2numpy(upstreamarea_pcr,mv=999999999)

# Create folder with corrected station locations
if os.path.isdir(corrected_locations_dir)==False:
    os.mkdir(corrected_locations_dir)


############################################################################
#   Loop over catchments to automatically snap station locations to the 
#   'correct' grid cell. If snapping doesn't work, a window is opened to 
#   manually select the correct location
############################################################################

catchment_dirs = glob.glob(os.path.join(databasedir,'*',""))
catchment_dirs = sorted(catchment_dirs)
for ii in np.arange(len(catchment_dirs)):

    catchment_dir = catchment_dirs[ii]  
    ID = os.path.split(os.path.dirname(catchment_dir))[-1]
    print("===============================================================================")
    print(str(ii)+" "+ID)
    t1 = time.time()
    
    if os.path.isfile(os.path.join(corrected_locations_dir,ID+'.txt')):
        print('Corrected station location already found, skipping')
        continue
    
    if 'Discharge' in globals(): del Discharge
    if 'StatLat' in globals(): del StatLat
    if 'StatLon' in globals(): del StatLon
    if 'Area' in globals(): del Area
    if 'CatchCentroidLat' in globals(): del CatchCentroidLat
    if 'CatchCentroidLon' in globals(): del CatchCentroidLon
    if 'coords' in globals(): del coords

    
    ############################################################################
    #   Load discharge and catchment boundary data
    ############################################################################
    
    # Load discharge data
    try:
        Discharge = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Discharge')
        StatLat = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lat')[0][0]
        StatLon = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/StationCoords/Lon')[0][0]
        StatRow,StatCol = latlon2rowcol(StatLat,StatLon,res,lat_upper,lon_left)
        Station = readmatfile(os.path.join(catchment_dir,"DISCHARGE.mat"),'DISCHARGE/Station')[0]
        Station = ''.join(map(chr,Station))
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
        
    snap_success = False
    

    ############################################################################
    #   If we have area and catchment boundary estimates, we will automatically 
    #   determine the most suitable grid cell
    ############################################################################
    
    if ('Area' in globals()) & ('CatchCentroidLat' in globals()) & ('CatchCentroidLat' in globals()):
    
        if Area<100000:             
            print("Catchment too small")
            continue

        # Snap settings
        permissible_area_error = 20 # Maximum allowed catchment area error (percent)
        permissible_centroid_error = 0.5 # Maximum allowed centroid location error (degrees)        
        margin = (np.ceil(np.sqrt(Area*0.1)/2)/(40000*0.08333/360)).astype(int) # Search box size (grid cells around station)
        if margin<2: margin = 2
        if margin>4: margin = 4
        print('Station snap margin set to '+str(margin))
        
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
            StatLatCorr,StatLonCorr = rowcol2latlon(StatRowCorr,StatColCorr,res,lat_upper,lon_left)
            print('Automatic snapping not necessary')
        else:
            rank_area_error_map = area_error_map.ravel().argsort().argsort().reshape(area_error_map.shape)
            rank_centroid_error_map = centroid_error_map.ravel().argsort().argsort().reshape(centroid_error_map.shape)
            rank_cum_map = rank_area_error_map+rank_centroid_error_map            
            StatRowCorr,StatColCorr = np.where(rank_cum_map==np.min(rank_cum_map))
            StatRowCorr,StatColCorr = StatRowCorr[0],StatColCorr[0]
            StatLatCorr,StatLonCorr = rowcol2latlon(StatRowCorr,StatColCorr,res,lat_upper,lon_left)
            print('Station automatically snapped')

        area_error = area_error_map[StatRowCorr,StatColCorr]
        centroid_error = centroid_error_map[StatRowCorr,StatColCorr]    
        print("Area error: "+str(area_error)+" %, centroid error: "+str(centroid_error)+" degrees")
        
        if (area_error<permissible_area_error):
            snap_success = True
            print('Automatic snap successful')
        else:
            print('Automatic snap failed')
        
    
    ############################################################################
    #   If we couldn't automatically snap, or if the snapping was unsuccessful,
    #   we have to manually select the correct grid cell
    ############################################################################
    
    plt.close('all')
    
    if snap_success==False:
        
        extent = [lon_left,lon_left+len(lon)*res,lat_upper-len(lat)*res,lat_upper] # (left, right, bottom, top)
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        
        fig.suptitle(str(ii)+' '+ID+' '+Station+', mean discharge '+str(np.nanmean(Discharge))+' m3/s')
        
        ax1.imshow(np.log10(upstreamarea_np), extent=extent)
        ax1.scatter(StatLon,StatLat,c='r')
        ax1.set_xlim(StatLon-6,StatLon+6)
        ax1.set_ylim(StatLat-6,StatLat+6)
        rivers_shp.plot(ax=ax1,column=None,cmap=None,facecolor="none",edgecolor='r',linewidth=0.35) 
                 
        ax2.imshow(np.log10(upstreamarea_np), extent=extent)
        ax2.scatter(StatLon,StatLat,c='r')
        ax2.set_xlim(StatLon-2,StatLon+2)
        ax2.set_ylim(StatLat-2,StatLat+2)
        rivers_shp.plot(ax=ax2,column=None,cmap=None,facecolor="none",edgecolor='r',linewidth=0.35) 

        ax3.imshow(np.log10(upstreamarea_np), extent=extent)
        ax3.scatter(StatLon,StatLat,c='r')
        ax3.set_xlim(StatLon-0.5,StatLon+0.5)
        ax3.set_ylim(StatLat-0.5,StatLat+0.5)
        rivers_shp.plot(ax=ax3,column=None,cmap=None,facecolor="none",edgecolor='r',linewidth=0.35) 
        
        margin = 6        
        ax4.imshow(np.log10(upstreamarea_np))
        ax4.scatter(StatCol,StatRow,c='r')
        ax4.set_xlim(StatCol-margin,StatCol+margin)
        ax4.set_ylim(StatRow-margin,StatRow+margin)
        plt.gca().invert_yaxis()
        plt.title('Click, then close')        
        
        # Select correct grid cell manually
        print('Click on correct grid cell, then close plot')
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        if 'coords' not in globals(): 
            print('No location selected, skipping')
            continue
        StatColCorr,StatRowCorr = coords[0],coords[1]
        StatLatCorr,StatLonCorr = rowcol2latlon(StatRowCorr,StatColCorr,res,lat_upper,lon_left)
        
    
    ############################################################################
    #   Show catchment of corrected location
    ############################################################################
          
    point_np = np.zeros(ldd_np.shape,dtype=bool)
    point_np[StatRowCorr,StatColCorr] = True
    point_pcr = pcr.numpy2pcr(pcr.Boolean,point_np,mv=-9999)
    catch_pcr = pcr.catchment(ldd_pcr, point_pcr)        
    catch_np = pcr.pcr2numpy(catch_pcr,mv=-9999)
    catch_np = catch_np==1
    extent = [lon_left,lon_left+len(lon)*res,lat_upper-len(lat)*res,lat_upper] # (left, right, bottom, top)
    plt.imshow(np.single(catch_np)+np.log10(upstreamarea_np)/5, extent=extent,vmin=0,vmax=3)
    #plt.plot(CatchBoundsLon,CatchBoundsLat)
    plt.scatter(StatLon,StatLat,c='r')
    plt.gca().set_xlim(StatLon-10,StatLon+10)
    plt.gca().set_ylim(StatLat-6,StatLat+6)        
    plt.title(str(ii)+' '+ID)
    plt.show(block=False)
    plt.pause(0.01)        


    ############################################################################
    #   Save corrected location
    ############################################################################
    
    np.savetxt(os.path.join(corrected_locations_dir,ID+'.txt'), np.array([StatColCorr,StatRowCorr]), delimiter=',',fmt='%1.3f')
    
    print('Time elapsed is ' + str(time.time() - t1) + ' sec')