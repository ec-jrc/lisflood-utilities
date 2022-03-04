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
    
def potential_evaporation(data,albedo,factor,doy,lat):

    # difference between daily maximum and minimum temperature [deg C]
    DeltaT = data['tmax']-data['tmin']
    DeltaT[DeltaT<0] = 0
    
    # empirical constant in windspeed formula
    # if DeltaT is less than 12 degrees, BU=0.54
    BU = 0.54+0.35*((DeltaT-12)/4)
    BU[BU<0.54] = 0.54
    
    # Goudriaan equation (1977)
    # saturated vapour pressure [mbar]
    # TAvg [deg Celsius]
    # exp is correct (e-power) (Van Der Goot, pers. comm 1999)    
    ESat = 6.10588*np.exp((17.32491*data['tmean'])/(data['tmean']+238.102))
    EAct = data['relhum']*ESat/100
    
    # Vapour pressure deficit [mbar]
    VapPressDef = ESat-EAct
    VapPressDef[VapPressDef<0] = 0

    # evaporative demand of reference vegetation canopy [mm/d]
    EA = 0.26*VapPressDef*(factor+BU*data['wind'])
    
    # latent heat of vaporization [MJ/kg]
    LatHeatVap = 2.501-0.002361*data['tmean']
    
    '''
    # psychrometric constant at sea level [mbar/deg C]
    # Corrected constant, was wrong originally
    # Psychro0 should be around 0.67 mbar/ deg C
    Psychro0 = 0.00163*(self.Press0/LatHeatVap)

    # Correction for altitude (FAO, http://www.fao.org/docrep/X0490E/x0490e00.htm )
    # Note that previously some equation from Supit et al was used,
    # but this produced complete rubbish!
    Psychro = Psychro0*((293-0.0065*self.Dem)/293) ** 5.26
    '''
    
    # Allen et al. (1994) equation 8 [mbar/deg C]
    Psychro = 10*(1.013*10**-3*data['pres']/10)/(0.622*LatHeatVap)
    
    # slope of saturated vapour pressure curve [mbar/deg C]
    Delta = (238.102*17.32491*ESat)/((data['tmean']+238.102)**2)
    

    # ************************************************************
    # ***** ANGOT RADIATION **************************************
    # ************************************************************

    # Latitude degrees to radians
    lat = lat*np.pi/180

    # solar declination [radians]
    Declin = -23.45*np.cos(np.deg2rad((360*(doy+10))/365))*np.pi/180
    
    # solar constant at top of the atmosphere [J/m2/s]
    SolarConstant = 1370*(1+(0.033*np.cos(np.deg2rad(2*np.pi*doy/365))))
    
    # daylength [hour]
    tmp1 = ((-np.sin(np.deg2rad(-2.65/np.pi)))+np.sin(Declin)*np.sin(lat))/(np.cos(Declin)*np.cos(lat))
    tmp2 = (tmp1<0)*np.rad2deg(np.arcsin(tmp1)-360)+(tmp1>=0)*np.rad2deg(np.arcsin(tmp1))
    DayLength = 12+(24/180)*tmp2
    pdb.set_trace()

    # Daylength equation can produce NaN at high latitudes,
    # replace NaNs with zeros
    DayLength = cover(DayLength, 0.0)

    # integral of solar height [s] over the day
    IntSolarHeight = 3600*(DayLength*np.sin(Declin)*np.sin(lat)+(24./np.pi)*np.cos(Declin)*np.cos(lat)*np.sqrt(1-sqr(np.tan(Declin)*np.tan(lat))))

    # Integral of solar height cannot be negative,
    # so truncate at 0
    IntSolarHeight = maximum(IntSolarHeight, 0.0)

    # Replace NaNs with zeros
    IntSolarHeight = cover(IntSolarHeight, 0.0)

    # daily extra-terrestrial radiation (Angot radiation) [J/m2/d]
    RadiationAngot = IntSolarHeight*SolarConstant
    

    # ************************************************************
    # ***** NET ABSORBED RADIATION *******************************
    # ************************************************************

    # equation Allen et al. 1994
    # using the digital elevation model
    # from:  An Update for the Definition of Reference Evapotranspiration  Allen et al. 1994
    Rso = RadiationAngot*(0.75+(2*10**-5*elev))
    TransAtm_Allen = data['swd']/Rso
    TransAtm_Allen = cover(TransAtm_Allen, 0)
    AdjCC = 1.8*TransAtm_Allen-0.35
    AdjCC[AdjCC<0] = 0.05
    AdjCC[AdjCC>1] = 1
    
    # Net emissivity
    EmNet = (0.56-0.079*np.sqrt(EAct))
    
    # net  longwave radiation [J/m2/day]
    RN = 4.903*10**-3*((data['tmean']+273)**4)*EmNet*AdjCC    

    # net absorbed radiation of reference vegetation canopy [mm/d]
    RNA = maximum(((1-albedo)*data['swd']-RN)/(10**6*LatHeatVap), 0.0)    

    # potential reference evapotranspiration rate [mm/day]
    return ((Delta*RNA)+(Psychro*EA))/(Delta+Psychro)
    