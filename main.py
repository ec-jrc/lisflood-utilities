#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tools import *
import rasterio
import gc

print('===============================================================================')
        
# Load configuration file
config = pd.read_csv(sys.argv[1],header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    print(string)
    string = string.replace(" ","")    
    try:
        exec(string.replace("=","=r"))     
    except:
        exec(string)         

# Create output folder
if os.path.isdir(output_folder)==False:
    os.mkdir(output_folder)
   
# Load template map
dset = Dataset(templatemap_path)
template_lat = dset.variables['lat'][:]
template_lon = dset.variables['lon'][:]
template_res = template_lon[1]-template_lon[0]
varname = list(dset.variables.keys())[-1]
template_np = np.array(dset.variables[varname][:])

# Determine map sizes
mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
mapsize_template = template_np.shape
row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)

# List of years with HYDE data
hyde_files = glob.glob(os.path.join(hyde_folder,'baseline','zip','cropland*'))
hyde_years = np.array([int(os.path.basename(hyde_file)[8:12]) for hyde_file in hyde_files])

# List of years with VCF data
vcf_files = glob.glob(os.path.join(vcf_folder,'*'))
vcf_years = np.array([int(os.path.basename(vcf_file)[8:12]) for vcf_file in vcf_files])

# List of years with GAIA data
# Do not use data prior to 1990
gaia_files = glob.glob(os.path.join(gaia_folder,'*.nc'))
gaia_years = np.array([int(os.path.basename(gaia_file)[0:4]) for gaia_file in gaia_files])
gaia_years = gaia_years[gaia_years>=1990]

# List of years with GSWE data
gswe_files = glob.glob(os.path.join(gswe_folder,'*'))
gswe_years = np.array([int(os.path.basename(gswe_file)[:4]) for gswe_file in gswe_files])
gswe_years = np.unique(gswe_years)

vars = ['fracwater','fracforest','fracsealed','fracrice','fracirrigation','fracother']
        

############################################################################
#   Loop over years and months to calculate fractions
############################################################################

for year in np.arange(year_start,year_end+1):
    #if year!=1979: continue
    print('-------------------------------------------------------------------------------')
    print('Year: '+str(year))
    
    t0 = time.time()
    idx = (np.abs(np.array(gaia_years)-year)).argmin()
    gaia_year = gaia_years[idx]
    print('Loading GAIA data ('+os.path.join(gaia_folder,str(gaia_year)+'.nc')+')')
    dset = Dataset(os.path.join(gaia_folder,str(gaia_year)+'.nc'))
    gaia_raw = np.array(dset.variables['impervious_fraction'][:]).clip(0,1)
    dset.close()
    print("Time elapsed is "+str(time.time()-t0)+" sec")

    t0 = time.time()
    print('Loading HILDA+ data ('+str(year)+')')
    ind = year-1899
    dset = Dataset(os.path.join(hildaplus_folder,'hildaplus_vGLOB-1.0-f_states.nc'))
    hilda_raw = np.array(dset.variables['LULC_states'][ind,:,:])
    dset.close()
    fracwater_hilda = imresize_mean(np.single((hilda_raw==0) | (hilda_raw==77)),mapsize_global).astype(np.double) 
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    for month in np.arange(1,13):
        #if month!=1: continue
        print('-------------------------------------------------------------------------------')
        print('Year: '+str(year)+' Month: '+str(month))
        t0 = time.time()
        
        t0 = time.time()
        print('Resampling GAIA data')
        fracsealed = 0.75*imresize_mean(gaia_raw,mapsize_global).astype(np.double)
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        t0 = time.time()
        idx = (np.abs(np.array(gswe_years)-year)).argmin()
        gswe_year = gswe_years[idx]
        gswe_month = month.copy()
        if (gswe_year==2018) & (gswe_month==10): gswe_month = 9 # October 2018 GSWE data are erroneous
        print('Loading and resampling GSWE fracwater data ('+os.path.join(gswe_folder,str(gswe_year)+'_'+str(gswe_month).zfill(2)+'_B.nc')+')')
        dset = Dataset(os.path.join(gswe_folder,str(gswe_year)+'_'+str(gswe_month).zfill(2)+'_B.nc'))
        gswe_lats = np.array(dset.variables['lat'][:])
        gswe_lons = np.array(dset.variables['lon'][:])
        gswe_res = gswe_lats[0]-gswe_lats[1]
        add_top = int(np.round((90-gswe_lats[0])/gswe_res))
        add_bottom = int(np.round((90+gswe_lats[-1])/gswe_res))
        gswe_shape = (int(len(gswe_lons)/2),len(gswe_lons))
        gswe_raw = np.zeros(gswe_shape,dtype=np.single)*np.NaN
        varname = list(dset.variables.keys())[-1] # Variable name differs for some years...
        gswe_raw[add_top+1:gswe_shape[0]-add_bottom,:] = dset.variables[varname][:].clip(0,1)
        fracwater = imresize_mean(gswe_raw,mapsize_global).astype(np.double)
        del gswe_raw
        dset.close()
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        print('Filling gaps in GSWE fracwater (at high latitudes) with HILDA+ fracwater')
        fracwater[np.isnan(fracwater)] = fracwater_hilda[np.isnan(fracwater)]        
        
        print('Making sure fracsealed+fracwater is <1')
        totals = fracsealed+fracwater
        mask = totals>1
        fracsealed[mask] = fracsealed[mask]*1/totals[mask]
        fracsealed = fracsealed.clip(0,1)
        fracwater[mask] = fracwater[mask]*1/totals[mask]
        fracwater = fracwater.clip(0,1)
        
        print('Computing fracother as residual of fracsealed and fracwater')
        fracother_init = 1-fracsealed-fracwater
        fracother_init = fracother_init.clip(0,1)
        
        t0 = time.time()
        idx = (np.abs(np.array(vcf_years)-year)).argmin()
        vcf_year = vcf_years[idx]
        vcf_path = glob.glob(os.path.join(vcf_folder,'*'+str(vcf_year)+'001*.tif'))[0]
        print('Loading MEaSUREs VCF data ('+vcf_path+')')
        src = rasterio.open(vcf_path)
        fracforest = src.read(1).astype(np.double)
        fracforest = (fracforest/100).clip(0,1)
        src.close()
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        t0 = time.time()
        idx = (np.abs(np.array(hyde_years)-year)).argmin()
        hyde_year = hyde_years[idx]
        print('Loading and resampling HYDE data ('+str(hyde_year)+')')
        hyde_garea_cr = np.array(pd.read_csv(os.path.join(hyde_folder,'general_files','garea_cr.asc'), 
            header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double) # total gridcell area in km2
        hyde_garea_cr = hyde_garea_cr[:,:-1] # Rogue last column due to spaces after last value in asc file...
        hyde_ir_norice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_norice'+str(hyde_year)+'AD.asc'), 
            header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # irrigated other crops area (no rice) area
        hyde_ir_norice[np.isnan(hyde_ir_norice)] = 0
        hyde_ir_norice = imresize_mean(hyde_ir_norice,mapsize_global)    
        hyde_ir_rice = np.array(pd.read_csv(os.path.join(hyde_folder,'baseline','zip','ir_rice'+str(hyde_year)+'AD.asc'), 
            header=None,sep=' ',skiprows=6,na_values=-9999).values,dtype=np.double)/hyde_garea_cr # irrigated rice area (no rice)
        hyde_ir_rice[np.isnan(hyde_ir_rice)] = 0
        hyde_ir_rice = imresize_mean(hyde_ir_rice,mapsize_global)
        fracrice = hyde_ir_rice.clip(0,1)
        fracirrigation = hyde_ir_norice.clip(0,1)
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        print('Reducing fracforest, fracirrigation, and fracrice if sum exceeds fracother')
        totals = fracforest+fracrice+fracirrigation
        mask = totals>fracother_init
        fracforest[mask] = fracforest[mask]*1/totals[mask]
        fracrice[mask] = fracrice[mask]*1/totals[mask]
        fracirrigation[mask] = fracirrigation[mask]*1/totals[mask]
        
        print('Making sure sum of fractions (without other) is <=1')
        total = fracwater+fracsealed+fracforest+fracrice+fracirrigation
        mask = total>1
        fracwater[mask] = fracwater[mask]/total[mask]
        fracsealed[mask] = fracsealed[mask]/total[mask]
        fracforest[mask] = fracforest[mask]/total[mask]
        fracrice[mask] = fracrice[mask]/total[mask]
        fracirrigation[mask] = fracirrigation[mask]/total[mask]
        
        print('Recalculating fracother as residual')
        fracother = 1-fracwater-fracsealed-fracforest-fracrice-fracirrigation
        
        print('Subsetting data to template map area')
        fracwater = fracwater[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        fracforest = fracforest[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        fracsealed = fracsealed[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        fracrice = fracrice[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        fracirrigation = fracirrigation[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        fracother = fracother[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)]
        
        print('Saving data to in Numpy format (folder '+output_folder+')')
        t0 = time.time()
        for vv in np.arange(len(vars)):            
            data = eval(vars[vv]) #.astype(np.single)
            #data = np.round(data*100)/100
            np.savez_compressed(os.path.join(output_folder,str(year)+str(month).zfill(2)+'_'+vars[vv]),data=data)
        print("Time elapsed is "+str(time.time()-t0)+" sec")

        gc.collect()


############################################################################
#   Convert to netCDF
############################################################################

for vv in np.arange(len(vars)):

    file = os.path.join(output_folder,vars[vv]+'.nc')
    varname = vars[vv]
    
    if os.path.isfile(file):
        os.remove(file)
        
    ncfile = Dataset(file, 'w', format='NETCDF4')
    ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

    ncfile.createDimension('lon', len(template_lon))
    ncfile.createDimension('lat', len(template_lat))
    ncfile.createDimension('time', None)

    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = template_lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'

    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'][:] = template_lat
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'].long_name = 'latitude'

    ncfile.createVariable('time', 'f8', 'time')
    ncfile.variables['time'].units = 'days since 1979-01-02 00:00:00'
    ncfile.variables['time'].long_name = 'time'
    ncfile.variables['time'].calendar = 'proleptic_gregorian'

    ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999) #, least_significant_digit=9
    ncfile.variables[varname].units = 'fraction'

    for year in np.arange(year_start,year_end+1):
        for month in np.arange(1,13):
            print('-------------------------------------------------------------------------------')
            print('Saving to netCDF var: '+varname+' year: '+str(year)+' month: '+str(month))
            t0 = time.time()
            
            data = np.load(os.path.join(output_folder,str(year)+str(month).zfill(2)+'_'+varname+'.npz'))['data']
            #data = np.round(data,decimals=3)
            
            index = (year-year_start)*12+month-1
             
            ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
            ncfile.variables[varname][index,:,:] = data
            
            print("Time elapsed is "+str(time.time()-t0)+" sec")
            
    ncfile.close()


############################################################################
#   Verify that sum of fractions are 1
############################################################################

print('-------------------------------------------------------------------------------')
print('Sum verification')

sum = np.zeros(mapsize_template)

for vv in np.arange(len(vars)):

    varname = vars[vv]        
    
    file = os.path.join(output_folder,varname+'.nc')    
    ncfile = Dataset(file)
    data = np.array(ncfile.variables[varname][0,:,:])
    
    #data = np.load(os.path.join(output_folder,'197901_'+varname+'.npz'))['data']
    
    sum = sum+data

print('Max sum of fractions: '+str(np.nanmax(sum)))
print('Min sum of fractions: '+str(np.nanmin(sum)))