#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "January 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from skimage.transform import resize
from skimage.transform import downscale_local_mean
from datetime import datetime, timedelta
from scipy import ndimage as nd
import rasterio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
        
def load_country_code_map(filepath,mapsize):
    
    src = rasterio.open(filepath)
    country_code_map = src.read(1).astype(np.single)
    country_code_map_refmatrix = src.get_transform()
    src.close()

    # Check if country border raster covers entire globe
    assert (country_code_map_refmatrix[0]==-180) & (country_code_map_refmatrix[3]==90)

    country_code_map = resize(country_code_map,mapsize,order=0,mode='edge',anti_aliasing=False)
    country_code_map[country_code_map==158] = 156 # Add Taiwan to China
    country_code_map[country_code_map==736] = 729 # South Sudan missing from map
    country_code_map[country_code_map==0] = np.NaN
    
    return country_code_map
    
def load_us_state_code_map(filepath,mapsize):
    
    src = rasterio.open(filepath)
    state_code_map = src.read(1).astype(np.single)
    state_code_map = resize(state_code_map,mapsize,order=0,mode='edge',anti_aliasing=False)
    state_code_map_refmatrix = src.get_transform()
    src.close()

    # Check if state border raster covers entire globe
    assert (state_code_map_refmatrix[0]==-180) & (state_code_map_refmatrix[3]==90)

    return state_code_map
    
def latlon2rowcol(lat,lon,res,lat_upper,lon_left):
    row = np.round((lat_upper-lat)/res-0.5).astype(int)
    col = np.round((lon-lon_left)/res-0.5).astype(int)    
    return row.squeeze(),col.squeeze()

def rowcol2latlon(row,col,res,lat_upper,lon_left):
    lat = lat_upper-row*res-res/2
    lon = lon_left+col*res+res/2
    return lat.squeeze(),lon.squeeze()
            
def imresize_mean(oldarray,newshape):
    '''Resample an array using simple averaging'''
    
    oldshape = oldarray.shape
    
    factor = oldshape[0]/newshape[0]
    
    if factor==int(factor):
        factor = int(factor)
        newarray = downscale_local_mean(oldarray,(factor,factor))
        
    else:
        factor = 1
        while newshape[0]*factor<oldshape[0]:
            factor = factor+1        
        intarray = resize(oldarray,(int(newshape[0]*factor),int(newshape[1]*factor)),order=0,mode='constant',anti_aliasing=False)
        newarray = downscale_local_mean(intarray,(factor,factor))

    return newarray
   
def fill(data, invalid=None):
    '''Nearest neighbor interpolation gap fill by Juh_'''
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]    

def load_config(filepath):
    '''Load configuration file into dict'''
    
    df = pd.read_csv(filepath,header=None,index_col=False)
    config = {}
    for ii in np.arange(len(df)): 
        string = df.iloc[ii,0].replace(" ","")
        varname = string.rpartition('=')[0]
        varcontents = string.rpartition('=')[2]
        try:
            varcontents = float(varcontents)
        except:
            pass
        config[varname] = varcontents
    return config
    
def initialize_netcdf(outfile,lat,lon,varname,units,least_significant_digit):
    
    ncfile = Dataset(outfile, 'w', format='NETCDF4')
    ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')
    ncfile.createDimension('lon', len(lon))
    ncfile.createDimension('lat', len(lat))
    ncfile.createDimension('time', None)
    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'
    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'][:] = lat
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'].long_name = 'latitude'
    ncfile.createVariable('time', 'f8', 'time')
    ncfile.variables['time'].units = 'days since 1979-01-02 00:00:00'
    ncfile.variables['time'].long_name = 'time'
    ncfile.variables['time'].calendar = 'proleptic_gregorian'
    ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'),zlib=True,
        chunksizes=(1,450,450,), fill_value=-9999,
        least_significant_digit=least_significant_digit)
    ncfile.variables[varname].units = units
    
    return ncfile

def potential_evaporation(data,albedo,factor,doy,lat,elev):
    """
    Calculate potential evaporation (mm/d) using an approach based on Penman-
    Monteith. More details provided in the LISVAP documentation (Van der 
    Knijff, 2006).
    
    Van der Knijff, J., 2006. LISVAP – Evaporation Pre-Processor for the 
        LISFLOOD Water Balance and Flood Simulation Model, User Manual. EUR 
        22639 EN, Office for Official Publications of the European 
        Communities, Luxembourg, 31 pp.
    
    INPUTS
    data:   Dict with grids of tmean, tmin, tmax (all in degrees Celsius), 
            relhum (%), wind (m/s), pres(mbar), swd (W/m2), and lwd (W/m2)
    albedo: Albedo (0 to 1)
    factor: Empirical factor related to land cover (>0)
    doy:    Day of year (1 to 366)
    lat:    Latitude grid (degrees)
    elev:   Elevation grid (m asl)
    
    OUTPUTS
    pet:    Potential evaporation (mm/d)
    """

    # Difference between daily maximum and minimum temperature (degrees C)
    DeltaT = data['tmax']-data['tmin']
    DeltaT[DeltaT<0] = 0
    
    # Empirical constant in windspeed formula (if DeltaT is less than 12 
    # degrees C, BU=0.54)
    BU = 0.54+0.35*((DeltaT-12)/4)
    BU[BU<0.54] = 0.54
    
    # Goudriaan equation (1977) to calculate saturated vapour pressure (mbar)
    ESat = 6.10588*np.exp((17.32491*data['tmean'])/(data['tmean']+238.102))
    
    # Actual vapor pressure calculated from relative humidity (mbar)
    EAct = data['relhum']*ESat/100
    
    # Vapour pressure deficit (mbar)
    VapPressDef = ESat-EAct
    VapPressDef[VapPressDef<0] = 0

    # Evaporative demand (mm/d)
    EA = {}
    for key in factor.keys():
        EA[key] = 0.26*VapPressDef*(factor[key]+BU*data['wind'])
        
    # Latent heat of vaporization (MJ/kg)
    LatHeatVap = 2.501-0.002361*data['tmean']
        
    # Allen et al. (1998) equation 8 (mbar/degrees C)
    Psychro = 10*(1.013*10**-3*data['pres']/10)/(0.622*LatHeatVap)
    
    # Slope of saturated vapour pressure curve (mbar/degrees C)
    Delta = (238.102*17.32491*ESat)/((data['tmean']+238.102)**2)
    
    
    #--------------------------------------------------------------------------
    #   Extra-terrestrial radiation
    #--------------------------------------------------------------------------

    # Solar declination (rad)
    Declin_rad = np.arcsin(0.39795*np.cos(0.2163108+2*np.arctan(0.9671396*np.tan(0.00860*(doy-186)))))
    
    # Convert latitude from degrees to radians
    lat_rad = lat*np.pi/180
    
    # Solar constant at top of the atmosphere (J/m2/s)
    SolarConstant = 1370*(1+0.033*np.cos(2*np.pi*doy/365))
    
    # Day length (h) equation from Forsythe et al. (1995; https://doi.org/10.1016/0304-3800(94)00034-F)
    sinLD = np.sin(Declin_rad)*np.sin(lat_rad)
    cosLD = np.cos(Declin_rad)*np.cos(lat_rad)
    DayLength = 24-(24/np.pi)*np.arccos((np.sin(0.8333*np.pi/180)+sinLD)/cosLD)
    DayLength = np.tile(fill(DayLength[:,:1]),(1,DayLength.shape[1])) # Nearest-neighbor gap filling
    
    # Integral of solar height over the day (s)
    IntSolarHeight = 3600*(DayLength*sinLD+(24/np.pi)*cosLD*np.sqrt(1-(sinLD/cosLD)**2))
    IntSolarHeight[IntSolarHeight<0] = 0
    IntSolarHeight = np.tile(fill(IntSolarHeight[:,:1]),(1,IntSolarHeight.shape[1])) # Nearest-neighbor gap filling
    
    # Daily extra-terrestrial radiation (J/m2/d)
    Ra = IntSolarHeight*SolarConstant
    
    # plt.imshow(Ra)
    # plt.colorbar()
    # plt.savefig('Ra.png',dpi=300)
    # plt.close()
    
    #--------------------------------------------------------------------------
    #   Net absorbed radiation
    #--------------------------------------------------------------------------

    # Clear-sky radiation (J/m2/d) from Allen et al. (1998; equation 37)
    Rso = Ra*(0.75+(2*10**-5*elev))
    
    # plt.imshow(Rso)
    # plt.colorbar()
    # plt.savefig('Rso.png',dpi=300)
    # plt.close()
    
    # Adjustment factor for cloud cover
    TransAtm_Allen = (data['swd']*86400+1)/(Rso+1)
    AdjCC = 1.8*TransAtm_Allen-0.35
    AdjCC[AdjCC<0.05] = 0.05
    AdjCC[AdjCC>1] = 1
    
    # plt.imshow(AdjCC)
    # plt.colorbar()
    # plt.savefig('AdjCC.png',dpi=300)    
    # plt.close()
    
    # plt.imshow(data['swd'])
    # plt.colorbar()
    # plt.savefig('swd.png',dpi=300)    
    # plt.close()    
    
    # Net emissivity
    EmNet = 0.56-0.079*np.sqrt(EAct)
    
    # Net longwave radiation (J/m2/d)
    StefBoltzConstant = 4.903*10**-3 # J/K4/m2/d
    RN = StefBoltzConstant*((data['tmean']+273.15)**4)*EmNet*AdjCC    

    # Net absorbed radiation of reference vegetation canopy (mm/d)
    RNA = {}
    for key in albedo.keys():
        RNA[key] = ((1-albedo[key])*data['swd']*86400-RN)/(10**6*LatHeatVap)
        RNA[key] = RNA[key].clip(0,None)
        
    # Potential reference evapotranspiration rate (mm/d)
    pet = {}
    for key in albedo.keys():
        pet[key] = ((Delta*RNA[key])+(Psychro*EA[key]))/(Delta+Psychro)
    
    # plt.imshow(pet)
    # plt.colorbar()
    # plt.savefig('pet.png',dpi=300)
    # plt.close()
    
    return pet
    
def makefig(folder,title,data,vmin,vmax):
    if os.path.isdir(folder)==False:
        os.makedirs(folder)
    plt.figure()
    ax = plt.gca()
    im = ax.imshow(data,vmin=vmin,vmax=vmax)
    ax.set_axis_off()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)   
    plt.colorbar(im, cax=cax)
    plt.title(title)
    plt.savefig(os.path.join(folder,title+'.png'),dpi=300,bbox_inches='tight')
    plt.close()