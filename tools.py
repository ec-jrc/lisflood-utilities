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
    
def potential_evaporation(data,albedo,cropf):


    DeltaT = maximum(data['tmax'] - data['tmin'], 0.0)
    # difference between daily maximum and minimum temperature [deg C]

    BU = maximum(0.54 + 0.35 * ((DeltaT - 12) / 4), 0.54)
    # empirical constant in windspeed formula
    # if DeltaT is less than 12 degrees, BU=0.54

    ESat = 6.10588 * exp((17.32491 * data['tmean']) / (data['tmean'] + 238.102))
    # Goudriaan equation (1977)
    # saturated vapour pressure [mbar]
    # TAvg [deg Celsius]
    # exp is correct (e-power) (Van Der Goot, pers. comm 1999)

    VapPressDef = maximum(ESat - self.EAct, 0.0)
    # Vapour pressure deficit [mbar]

    # Evaporative demand is calculated for three reference surfaces:
    #
    # 1. Reference vegetation canopy
    # 2. Bare soil surface
    # 3. Open water surface
    EA = 0.26 * VapPressDef * (self.FactorCanopy + BU * data['wind'])
    # evaporative demand of reference vegetation canopy [mm/d]
    EASoil = 0.26 * VapPressDef * (self.FactorSoil + BU * data['wind'])
    # evaporative demand of bare soil surface [mm/d]
    EAWater = 0.26 * VapPressDef * (self.FactorWater + BU * data['wind'])
    # evaporative demand of water surface [mm/d]

    LatHeatVap = 2.501 - 0.002361 * data['tmean']
    # latent heat of vaporization [MJ/kg]

    Psychro0 = 0.00163 * (self.Press0 / LatHeatVap)
    # psychrometric constant at sea level [mbar/deg C]
    # Corrected constant, was wrong originally
    # Psychro0 should be around 0.67 mbar/ deg C

    Psychro = Psychro0 * ((293 - 0.0065 * self.Dem) / 293) ** 5.26
    # Correction for altitude (FAO, http://www.fao.org/docrep/X0490E/x0490e00.htm )
    # Note that previously some equation from Supit et al was used,
    # but this produced complete rubbish!

    Delta = (238.102 * 17.32491 * ESat) / ((data['tmean'] + 238.102) ** 2)
    # slope of saturated vapour pressure curve [mbar/deg C]

    # ************************************************************
    # ***** ANGOT RADIATION **************************************
    # ************************************************************

    Declin = -23.45 * cos((360. * (self.calendar_day + 10)) / 365.)
    # solar declination [degrees]

    SolarConstant = self.AvSolarConst * (1 + (0.033 * np.cos(2 * self.Pi * self.calendar_day / 365.)))
    # solar constant at top of the atmosphere [J/m2/s]

    tmp1 = ((-sin(self.PD / self.Pi)) + sin(Declin) * sin(self.Lat)) / (cos(Declin) * cos(self.Lat))
    tmp2 = ifthenelse(tmp1 < 0, scalar(asin(tmp1)) - 360., scalar(asin(tmp1)))
    DayLength = 12. + (24. / 180.) * tmp2
    # daylength [hour]

    DayLength = cover(DayLength, 0.0)
    # Daylength equation can produce MV at high latitudes,
    # this statements sets day length to 0 in that case

    IntSolarHeight = 3600. * (DayLength * sin(Declin) * sin(self.Lat) + (24./self.Pi) * cos(Declin) * cos(self.Lat) * sqrt(1 - sqr(tan(Declin) * tan(self.Lat))))
    # integral of solar height [s] over the day

    IntSolarHeight = maximum(IntSolarHeight, 0.0)
    # Integral of solar height cannot be negative,
    # so truncate at 0
    IntSolarHeight = cover(IntSolarHeight, 0.0)

    RadiationAngot = IntSolarHeight * SolarConstant
    # daily extra-terrestrial radiation (Angot radiation) [J/m2/d]

    # ************************************************************
    # ***** NET ABSORBED RADIATION *******************************
    # ************************************************************

    # equation Allen et al. 1994
    # using the digital elevation model
    # from:  An Update for the Definition of Reference Evapotranspiration  Allen et al. 1994

    Rds = self.Rds
    Rso = RadiationAngot * (0.75 + (2 * 10 ** -5 * self.Dem))
    TransAtm_Allen = Rds/Rso
    TransAtm_Allen = cover(TransAtm_Allen, 0)
    AdjCC = 1.8 * TransAtm_Allen - 0.35
    AdjCC = ifthenelse(AdjCC < 0, 0.05, AdjCC)
    AdjCC = ifthenelse(AdjCC > 1, 1, AdjCC)

    EmNet = (0.56 - 0.079 * sqrt(self.EAct))
    # Net emissivity
    RN = self.StefBolt * ((data['tmean'] + 273) ** 4) * EmNet * AdjCC
    # net  longwave radiation [J/m2/day]

    RNA = maximum(((1 - self.AlbedoCanopy) * Rds - RN) / (1E6 * LatHeatVap), 0.0)
    # net absorbed radiation of reference vegetation canopy [mm/d]
    RNASoil = maximum(((1 - self.AlbedoSoil) * Rds - RN) / (1E6 * LatHeatVap), 0.0)
    # net absorbed radiation of bare soil surface
    RNAWater = maximum(((1 - self.AlbedoWater) * Rds - RN) / (1E6 * LatHeatVap), 0.0)
    # net absorbed radiation of water surface

    # ************************************************************
    # ***** EA: EVAPORATIVE DEMAND *******************************
    # ************************************************************
    # Evaporative demand is calculated for three reference surfaces:
    # 1. Reference vegetation canopy
    # 2. Bare soil surface
    # 3. Open water surface
    self.ETRef = ((Delta * RNA) + (Psychro * EA)) / (Delta + Psychro)
    # potential reference evapotranspiration rate [mm/day]
    self.ESRef = ((Delta * RNASoil) + (Psychro * EASoil)) / (Delta + Psychro)
    # potential evaporation rate from a bare soil surface [mm/day]
    self.EWRef = ((Delta * RNAWater) + (Psychro * EAWater)) / (Delta + Psychro)
    # potential evaporation rate from water surface [mm/day]
