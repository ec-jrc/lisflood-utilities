# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:31:52 2022

@author: tilloal
"""

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
  

__OUTPUT_FILE_EXT = '.nc'

__NETCDF_DATASET_FORMAT = 'NETCDF4_CLASSIC'
__NETCDF_CONVENTIONS = 'CF-1.6'
__NETCDF_SOURCE_SOFTWARE = 'Python netCDF4'

__NETCDF_VAR_TIME_DIMENSION = None
__NETCDF_VAR_TIME_CALENDAR_TYPE = 'proleptic_gregorian'

__NETCDF_VAR_DATA_TYPE = 'f8'
__NETCDF_VALUE_DATA_TYPE = 'f4'
__NETCDF_COORDINATES_DATA_TYPE = 'i4'

__KEY_STANDARD_NAME = 'value_standard_name'
__KEY_LONG_NAME = 'value_long_name'
__KEY_UNIT = 'value_unit'
__KEY_OFFSET = 0
__KEY_SCALE_FACTOR = 1
__KEY_VMIN = -400
__KEY_VMAX = 400

__meteo_vars_config = {
   
    'tp' : {__KEY_UNIT : 'mm', __KEY_STANDARD_NAME : 'tp', __KEY_LONG_NAME : 'total_precipitation',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -1, __KEY_VMAX : 7000},
    'ta' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'ta', __KEY_LONG_NAME : 'mean_temperature',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
    'td' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'td', __KEY_LONG_NAME : 'mean_dewpoint_temperature',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
    'ws' : {__KEY_UNIT : 'm/s', __KEY_STANDARD_NAME : 'ws', __KEY_LONG_NAME : 'avg_wind_speed',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : 0, __KEY_VMAX : 45},
    'u10' : {__KEY_UNIT : 'm/s', __KEY_STANDARD_NAME : 'u10', __KEY_LONG_NAME : 'avg_u_component_wind',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : 0, __KEY_VMAX : 45},
    'v10' : {__KEY_UNIT : 'm/s', __KEY_STANDARD_NAME : 'ws', __KEY_LONG_NAME : 'avg_v_component_wind',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : 0, __KEY_VMAX : 45},
    'rgd' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'ssrd', __KEY_LONG_NAME : 'surface_downward_solar_radiation',
        __KEY_OFFSET :0.0, __KEY_SCALE_FACTOR : 10000.0},
    'rn' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'str', __KEY_LONG_NAME : 'surface_net_thermal_radiation',
            __KEY_OFFSET : 0.0, __KEY_SCALE_FACTOR : 10000.0},
}

      
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
        least_significant_digit=least_significant_digit,complevel=4)
    ncfile.variables[varname].units = units
    
    ncfile.scale_factor=__meteo_vars_config[varname][__KEY_SCALE_FACTOR]       
    ncfile.add_offset=__meteo_vars_config[varname][__KEY_OFFSET]
    #add projection system 
    proj = ncfile.createVariable('wsg_1984', 'i4')
    proj.grid_mapping_name = 'latitude_longitude'
    #proj.false_easting= ''
    #proj.false_northing= ''
    #proj.longitude_of_projection_origin= ''
    #proj.latitude_of_projection_origin= ''
    proj.semi_major_axis= '6378137.0'
    proj.inverse_flattening='298.257223563'
    proj.proj4_params='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    proj.EPSG_code='EPSG:4326'
    
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
    
    
    #--------------------------------------------------------------------------
    #   Net absorbed radiation
    #--------------------------------------------------------------------------

    # Clear-sky radiation (J/m2/d) from Allen et al. (1998; equation 37)
    Rso = Ra*(0.75+(2*10**-5*elev))
    
    # Adjustment factor for cloud cover
    TransAtm_Allen = (data['swd']*86400+1)/(Rso+1)
    AdjCC = 1.8*TransAtm_Allen-0.35
    AdjCC[AdjCC<0.05] = 0.05
    AdjCC[AdjCC>1] = 1
    
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