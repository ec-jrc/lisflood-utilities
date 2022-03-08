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
import math
import matplotlib.pyplot as plt
        
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
    
def potential_evaporation(data,albedo,factor,doy,lat,elev):
    """
    # Inputs
    data = dict with grids of tmean, tmin, tmax, relhum, wind, pres, swd, lwd
    albedo = albedo (0 to 1)
    factor = empirical factor related to land cover (>0)
    doy = day of year (1 to 366)
    lat = latitude (degrees)
    elev = elevation (m asl)
    
    # Output
    Potential evaporation (mm/d)
    """

    #from pcraster.operations import sin,cos,tan,asin,sqrt,sqr,max,ifthenelse,exp
    #import pcraster as pcr
    #pcr.setclone(3600,7200,0.05,-180,90)
    
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
        
    # Allen et al. (1994) equation 8 (mbar/deg C)
    Psychro = 10*(1.013*10**-3*data['pres']/10)/(0.622*LatHeatVap)
    
    # Slope of saturated vapour pressure curve (mbar/deg C)
    Delta = (238.102*17.32491*ESat)/((data['tmean']+238.102)**2)
    

    # ************************************************************
    # ***** ANGOT RADIATION **************************************
    # ************************************************************

    # Solar declination (rad)
    Declin_rad = math.asin(0.39795 * np.cos(0.2163108 + 2 * math.atan(0.9671396 * math.tan(.00860 * (doy - 186)))))
    Declin_deg = Declin_rad*180/np.pi
    #Declin_deg_pcr = pcr.numpy2pcr(pcr.Scalar,Declin_deg,mv=-9999)
    
    lat_deg = lat
    lat_rad = lat*np.pi/180    
    #lat_deg_pcr = pcr.numpy2pcr(pcr.Scalar,lat_deg,mv=-9999)
    
    # Solar constant at top of the atmosphere (J/m2/s)
    SolarConstant = 1370*(1+0.033*np.cos(2*np.pi*doy/365))
    
    # Day length (h)
    # Ecological Modeling, volume 80 (1995) pp. 87-95
    # "A Model Comparison for Daylength as a Function of Latitude and Day of the Year"    
    DayLength = 24 - (24 / np.pi) * np.arccos((np.sin(0.8333 * np.pi / 180) + np.sin(lat_rad) * np.sin(Declin_rad)) / (np.cos(lat_rad) * np.cos(Declin_rad)))
    DayLength = np.tile(fill(DayLength[:,:1]),(1,DayLength.shape[1]))
    
    plt.imshow(DayLength)
    plt.savefig('DayLength.png',dpi=300)
    plt.close()          
    
    #DayLength_pcr = pcr.numpy2pcr(pcr.Scalar,DayLength,mv=-9999)
    
    # Integral of solar height over the day (s)
    #int_solar_height = 3600. * (DayLength_pcr * sin(Declin_deg_pcr) * sin(lat_deg_pcr) + (24./3.14) * cos(Declin_deg_pcr) * cos(lat_deg_pcr) * sqrt(1 - sqr(tan(Declin_deg_pcr) * tan(lat_deg_pcr))))
    sinLD = np.sin(Declin_rad)*np.sin(lat_rad)
    cosLD = np.cos(Declin_rad)*np.cos(lat_rad)
    IntSolarHeight = 3600*(DayLength*sinLD+(24/np.pi)*cosLD*np.sqrt(1-(sinLD/cosLD)**2))
    IntSolarHeight[IntSolarHeight<0] = 0
    IntSolarHeight = np.tile(fill(IntSolarHeight[:,:1]),(1,IntSolarHeight.shape[1]))    
    
    # Daily extra-terrestrial radiation (Angot radiation) (J/m2/d)
    RadiationAngot = IntSolarHeight*SolarConstant
    
    plt.imshow(RadiationAngot)
    plt.savefig('RadiationAngot.png',dpi=300)
    plt.close()          
    

    # ************************************************************
    # ***** NET ABSORBED RADIATION *******************************
    # ************************************************************

    # equation Allen et al. 1994
    # using the digital elevation model
    # from:  An Update for the Definition of Reference Evapotranspiration  Allen et al. 1994
    Rso = RadiationAngot*(0.75+(2*10**-5*elev))/86400
    TransAtm_Allen = (data['swd']+0.001)/(Rso+0.001)
    AdjCC = 1.8*TransAtm_Allen-0.35
    AdjCC[AdjCC<0] = 0.05
    AdjCC[AdjCC>1] = 1
    
    
    plt.imshow(TransAtm_Allen,vmin=-1,vmax=2)
    plt.savefig('TransAtm_Allen.png',dpi=300)
    plt.close()          
    
    plt.imshow(AdjCC,vmin=0,vmax=1)
    plt.savefig('AdjCC.png',dpi=300)
    plt.close()          
    
    # Net emissivity
    EmNet = (0.56-0.079*np.sqrt(EAct))
    
    # net  longwave radiation [J/m2/day]
    RN = 4.903*10**-3*((data['tmean']+273.15)**4)*EmNet*AdjCC    

    # net absorbed radiation of reference vegetation canopy [mm/d]
    RNA = ((1-albedo)*data['swd']-RN)/(10**6*LatHeatVap)
    RNA[RNA<0] = 0

    # potential reference evapotranspiration rate [mm/day]
    pet = ((Delta*RNA)+(Psychro*EA))/(Delta+Psychro)
    
    plt.imshow(pet,vmin=0,vmax=10)
    plt.savefig('pet.png',dpi=300)
    plt.close()          
    
    
    pdb.set_trace()
    
    
    
    return pet
    