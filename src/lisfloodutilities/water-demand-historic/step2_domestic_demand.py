#!/usr/bin/env python  
# -*- coding: utf-8 -*-
"""
Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""
__author__ = "Hylke E. Beck"
__date__ = "January 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
from skimage.transform import resize
from tools import *
import rasterio
from calendar import monthrange

# Load configuration file
config = load_config(sys.argv[1])

def main():
    print('===============================================================================')
           
    # Create output folders
    if os.path.isdir(os.path.join(config['output_folder'],'step2_domestic_demand','tables'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step2_domestic_demand','tables'))
    if os.path.isdir(os.path.join(config['output_folder'],'step2_domestic_demand','figures'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step2_domestic_demand','figures'))
       
    # Load template map
    dset = Dataset(config['templatemap_path'])
    template_lat = dset.variables['lat'][:]
    template_lon = dset.variables['lon'][:]
    template_res = template_lon[1]-template_lon[0]
    varname = list(dset.variables.keys())[-1]
    template_np = np.array(dset.variables[varname][:])

    # Determine map sizes
    mapsize_global = (np.round(180/template_res).astype(int),np.round(360/template_res).astype(int))
    mapsize_template = template_np.shape
    row_upper,col_left = latlon2rowcol(template_lat[0],template_lon[0],template_res,90,-180)

    # Compute area for each grid-cell (includes correction because lat/lon 
    # values in templates are often not rounded...
    xi, yi = np.meshgrid(np.arange(-180+template_res/2,180+template_res/2,template_res), np.arange(90-template_res/2,-90-template_res/2,-template_res))
    if yi.shape[0]>np.round(180/template_res):
        yi, xi = yi[:-1,:], xi[:-1,:]
    if yi.shape[1]>np.round(360/template_res):
        yi, xi = yi[:,:-1], xi[:,:-1]
    area_map = (40075*template_res/360)**2*np.cos(np.deg2rad(yi))

    # List of years
    years = np.arange(config['year_start'],config['year_end']+1).astype(int)

    # List of countries and US states
    country_codes = pd.read_csv(os.path.join(config['ancillary_data_folder'],'un_country_codes.csv'))
    state_codes = pd.read_csv(os.path.join(config['ancillary_data_folder'],'us-state-ansi-fips.csv'))
    

    ############################################################################
    #   Load and reproject Huang et al. (2018) water demand data (only used for
    #   comparison)
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and reprojecting Huang et al. (2018) water demand data')
    t0 = time.time()

    dset = Dataset(os.path.join(config['huang_folder'],'withd_dom.nc'))
    raw_data = np.array(dset.variables['withd_dom'][:])
    raw_lat = np.array(dset.variables['lat'][:])
    raw_lon = np.array(dset.variables['lon'][:])
    dset.close()

    Huang_withdrawal = np.zeros((360,720,raw_data.shape[0]),dtype=np.single)*np.NaN
    rows,cols = latlon2rowcol(raw_lat,raw_lon,0.5,90,-180)
    for ii in np.arange(raw_data.shape[0]):    
        reprojected = np.zeros((360,720),dtype=np.single)*np.NaN
        reprojected[rows,cols] = raw_data[ii,:]
        Huang_withdrawal[:,:,ii] = reprojected
        
    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load country border raster
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and resampling country border raster')
    t0 = time.time()
    country_code_map = load_country_code_map(os.path.join(config['world_borders_folder'],'CNTR_RG_01M_2024_4326_rasterized.tif'),mapsize_global)
    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load US states raster
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and resampling US state border raster')
    t0 = time.time()
    state_code_map = load_us_state_code_map(os.path.join(config['us_states_folder'],'cb_2023_us_state_500k_rasterized.tif'),mapsize_global)
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    

    ############################################################################
    #   Map of R parameter: measures the relative difference of domestic water
    #   withdrawal between the warmest and coldest months in a given year.
    #   Values from Huang et al. (2018).
    ############################################################################

    R_map = np.zeros(mapsize_global,dtype=np.single)+0.45
    R_map[country_code_map==124] = 0.36 # Canada
    R_map[country_code_map==840] = 0.52 # US
    R_map[country_code_map==36] = 0.8 # Australia
    R_map[country_code_map==356] = 0.29 # India
    R_map[country_code_map==156] = 0.2 # China
    R_map[country_code_map==392] = 0.1 # Japan
    R_map[country_code_map==724] = 0.1 # Spain


    ############################################################################
    #   Compute time series of population for each country
    ############################################################################

    filepath = os.path.join(config['output_folder'],'step2_domestic_demand','tables','population.csv')

    if os.path.isfile(filepath):
        table_population = pd.read_csv(filepath,index_col=0).values

    else:    
        print('-------------------------------------------------------------------------------')
        print('Computating population for each country and year')
        t0 = time.time()

        table_population = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
        for ii in np.arange(len(years)):
            print(str(years[ii]))
            npz = np.load(os.path.join(config['output_folder'],'step1_population_density',str(years[ii])+'.npz'))    
            pop_map = npz['data']
            for jj in np.arange(country_codes.shape[0]):
                country_code = country_codes.iloc[jj]['country-code']
                mask = country_code_map==country_code
                table_population[jj,ii] = np.sum(pop_map[mask]*area_map[mask])+1 # +1 to avoid divide by zero

        pd.DataFrame(table_population,index=country_codes['name'],columns=years).to_csv(filepath)

        print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Temporally inter- and extrapolate AQUASTAT data using population 
    #   estimates
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Temporally inter- and extrapolating AQUASTAT data')
    t0 = time.time()

    aquastat = pd.read_csv(os.path.join(config['aquastat_folder'],'aquastat_clean.csv'),index_col=False)

    table_aquastat_domestic_withdrawal_interp = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']

        # Find location of data in AQUASTAT table
        sel = np.where((aquastat['m49']==country_code) & (aquastat['Variable Name']=='Municipal water withdrawal'))[0]
        if len(sel)==0: 
            continue
            
        fao_years = aquastat['Year'][sel].values
        fao_withdrawals = aquastat['Value'][sel].values # km3
        keep = ~np.isnan(fao_withdrawals)
        fao_years = fao_years[keep]
        fao_withdrawals = fao_withdrawals[keep]
        if len(fao_withdrawals)==0:
            continue
        
        # Interpolate between years
        table_aquastat_domestic_withdrawal_interp[jj,:] = np.interp(years,fao_years,fao_withdrawals) # km3/year
        
        # If necessary, extrapolate backwards based on population
        if fao_years[0]>years[0]:
            ind = np.argwhere(years==fao_years[0])[0][0]
            table_aquastat_domestic_withdrawal_interp[jj,:ind] = table_population[jj,:ind]*table_aquastat_domestic_withdrawal_interp[jj,ind]/table_population[jj,ind]
        
        # If necessary, extrapolate forwards based on population
        if fao_years[-1]<years[-1]:
            ind = np.argwhere(years==fao_years[-1])[0][0]
            table_aquastat_domestic_withdrawal_interp[jj,ind:] = table_population[jj,ind:]*table_aquastat_domestic_withdrawal_interp[jj,ind]/table_population[jj,ind]

    # Save to csv
    pd.DataFrame(table_aquastat_domestic_withdrawal_interp,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','aquastat_domestic_withdrawal_interp.csv'))
    
    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load USGS water use data
    ############################################################################

    # I FOUND DATA TABLES FOR 1985 AND 1990 ON THE USGS WEBSITE, DELETE THIS?
    # The USGS National Water Information System (NWIS) does not provide 
    # thermoelectric withdrawal estimates prior to 1995. We therefore take 
    # estimates from Table 7 of USGS Circular 1001 (available at 
    # https://pubs.usgs.gov/circ/1983/1001/report.pdf). Manufacturing withdrawal
    # estimates are not available before 1985 from the NWIS, which we have to 
    # accept, because manufacturing estimates from USGS Circular 1001 appear to 
    # be inconsistent. Ideally, we would also use thermoelectric withdrawal 
    # estimates for 1985 and 1990 but this should do for now. NOT USED BECAUSE
    #USGS_1980_table_7_part_1 = pd.read_csv(os.path.join(config['ancillary_data_folder'],'1980_table_7_industrial_water_use_part_1.csv'),skiprows=2,header=None)
    #USGS_1980_table_7_part_2 = pd.read_csv(os.path.join(config['ancillary_data_folder'],'1980_table_7_industrial_water_use_part_2.csv'),skiprows=2,header=None)
    
    # Load USGS data tables for 1985 and 1990. These values correspond to the 
    # values published in the USGS circulars. However, the total livestock 
    # withdrawal, for example, is much higher than the year 2000 total shown in a 
    # graph on the USGS website 
    # (https://www.usgs.gov/mission-areas/water-resources/science/livestock-water-use),
    # for some reason...
    USGS_1985_data_table = pd.read_csv(os.path.join(config['usgs_water_use_folder'],'us85st.txt'),sep='\t',index_col=None)
    USGS_1990_data_table = pd.read_csv(os.path.join(config['usgs_water_use_folder'],'us90st.txt'),sep='\t',index_col=None)
    
    # List of USGS NWIS water use files
    files = glob.glob(os.path.join(config['usgs_water_use_folder'],'water_use*'))
    for ii in np.arange(len(files)):
        file = files[ii]
        df_state = pd.read_csv(file,comment='#',skiprows=2, sep='\t')
        df_state.drop(df_state.index[0],inplace=True)
        df_state.replace("-", "", regex=True, inplace=True)
        df_state.replace(r'^\s*$', np.nan,regex=True, inplace=True) # Turn empty strings into NaN
        df_state = df_state.astype(object)
        if ii==0: 
            df_country = df_state
        else:
            df_country = pd.merge(df_country,df_state,how="outer")
    
    # Loop over states and load water use data
    table_usgs_domestic_withdrawal = np.zeros((state_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_usgs_manufacturing_withdrawal = np.zeros((state_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_usgs_thermoelectric_withdrawal = np.zeros((state_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_usgs_thermoelectric_generated = np.zeros((state_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_usgs_livestock_withdrawal = np.zeros((state_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(state_codes.shape[0]):
        state_id = state_codes['st'][jj]
        state_name = state_codes['stname'][jj]
        state_acr = state_codes['stusps'][jj]
        sel = df_country['state_cd'].astype(int)==state_id
        
        # Load NWIS water use data
        values_dom = df_country['Domestic total self-supplied withdrawals, fresh, in Mgal/d'][sel].values.astype(float)+df_country['Public Supply total self-supplied withdrawals, fresh, in Mgal/d'][sel].values.astype(float)
        values_manuf = df_country['Industrial total self-supplied withdrawals, fresh, in Mgal/d'][sel].values.astype(float)
        values_thermo = df_country['Total Thermoelectric Power total self-supplied withdrawals, fresh, in Mgal/d'][sel].values.astype(float)
        values_gen = df_country['Total Thermoelectric Power power generated, in gigawatt-hours'][sel].values.astype(float)
        values_liv = df_country['Livestock total self-supplied withdrawals, fresh, in Mgal/d'][sel].values.astype(float)
        yrs = df_country['year'][sel].values.astype(int)    
        for yr, value_dom, value_manuf, value_thermo, value_gen, value_liv in zip(yrs,values_dom,values_manuf,values_thermo,values_gen,values_liv):
            try:
                ii = years.tolist().index(yr)
            except:
                continue                
            table_usgs_domestic_withdrawal[jj,ii] = 365.25*value_dom*4.54609*10**-6 # Mgal/d to km3/yr
            table_usgs_manufacturing_withdrawal[jj,ii] = 365.25*value_manuf*4.54609*10**-6 # Mgal/d to km3/yr
            table_usgs_thermoelectric_withdrawal[jj,ii] = 365.25*value_thermo*4.54609*10**-6 # Mgal/d to km3/yr
            table_usgs_thermoelectric_generated[jj,ii] = value_gen
            table_usgs_livestock_withdrawal[jj,ii] = 365.25*value_liv*4.54609*10**-6 # Mgal/d to km3/yr
        
        '''
        # Load thermoelectric withdrawal estimates for 1980 from USGS circular
        ii = years.tolist().index(1980)
        row = USGS_1980_table_7_part_2.iloc[:,0].tolist().index(state_name)
        table_usgs_thermoelectric_withdrawal[jj,ii] = 365.25*USGS_1980_table_7_part_2.iloc[row,4]*4.54609*10**-6 # Mgal/d to km3/yr
        '''
        
        # Load withdrawal estimates for 1985 from USGS table
        ii = years.tolist().index(1985)
        row = USGS_1985_data_table.iloc[:,0].tolist().index(state_acr)
        col1 = USGS_1985_data_table.columns.tolist().index('do-sstot')
        col2 = USGS_1985_data_table.columns.tolist().index('do-total')        
        table_usgs_domestic_withdrawal[jj,ii] = 365.25*(USGS_1985_data_table.iloc[row,col1]+USGS_1985_data_table.iloc[row,col2])*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1985_data_table.columns.tolist().index('in-wtofr')
        table_usgs_manufacturing_withdrawal[jj,ii] = 365.25*USGS_1985_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1985_data_table.columns.tolist().index('pt-frtot')
        table_usgs_thermoelectric_withdrawal[jj,ii] = 365.25*USGS_1985_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1985_data_table.columns.tolist().index('lv-total')
        table_usgs_livestock_withdrawal[jj,ii] = 365.25*USGS_1985_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr
        
        # Load withdrawal estimates for 1990 from USGS table
        ii = years.tolist().index(1990)
        row = USGS_1990_data_table.iloc[:,0].tolist().index(state_acr)
        col1 = USGS_1990_data_table.columns.tolist().index('do-sstot')
        col2 = USGS_1990_data_table.columns.tolist().index('do-total')        
        table_usgs_domestic_withdrawal[jj,ii] = 365.25*(USGS_1990_data_table.iloc[row,col1]+USGS_1990_data_table.iloc[row,col2])*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1990_data_table.columns.tolist().index('in-wtofr')
        table_usgs_manufacturing_withdrawal[jj,ii] = 365.25*USGS_1990_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1990_data_table.columns.tolist().index('pt-frtot')
        table_usgs_thermoelectric_withdrawal[jj,ii] = 365.25*USGS_1990_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr
        col = USGS_1990_data_table.columns.tolist().index('lv-total')
        table_usgs_livestock_withdrawal[jj,ii] = 365.25*USGS_1990_data_table.iloc[row,col]*4.54609*10**-6 # Mgal/d to km3/yr

    # Save to csv
    pd.DataFrame(table_usgs_domestic_withdrawal,index=state_codes['stname'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_domestic_withdrawal.csv'))
    pd.DataFrame(table_usgs_manufacturing_withdrawal,index=state_codes['stname'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_manufacturing_withdrawal.csv'))
    pd.DataFrame(table_usgs_thermoelectric_withdrawal,index=state_codes['stname'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_thermoelectric_withdrawal.csv'))
    pd.DataFrame(table_usgs_thermoelectric_generated,index=state_codes['stname'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_thermoelectric_generated.csv'))
    pd.DataFrame(table_usgs_livestock_withdrawal,index=state_codes['stname'],columns=years).to_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_livestock_withdrawal.csv'))

    # Linearly interpolate and nearest-neighbor extrapolate USGS domestic withdrawals
    table_usgs_domestic_withdrawal_interp = pd.DataFrame(table_usgs_domestic_withdrawal).interpolate(axis=1,method="linear",fill_value="extrapolate", limit_direction="both").values

        
    ############################################################################
    #   Disaggregate country scale withdrawal estimates using population data
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Disaggregating country withdrawal estimates using population data')

    for ii in np.arange(len(years)):
        t0 = time.time()
        year = years[ii]
        
        #if os.path.isfile(os.path.join(config['output_folder'],'step2_domestic_demand',str(year)+'12.npz')): continue
            
        print('Processing '+str(year))

        # Load population data
        npz = np.load(os.path.join(config['output_folder'],'step1_population_density',str(year)+'.npz'))
        pop_map = npz['data']

        # Load MSWX air temperature data
        mswx_monthly = np.zeros((mapsize_global[0],mapsize_global[1],12),dtype=np.single)*np.NaN
        for month in np.arange(1,13):
            dset = Dataset(os.path.join(config['mswx_folder'],'Past','Temp','Monthly',str(year)+str(month).zfill(2)+'.nc'))
            data = np.squeeze(np.array(dset.variables['air_temperature']))
            mswx_monthly[:,:,month-1] = imresize_mean(data,mapsize_global)
        mswx_avg = np.mean(mswx_monthly,axis=2)
        mswx_max = np.max(mswx_monthly,axis=2)
        mswx_min = np.min(mswx_monthly,axis=2)
        
        
        #--------------------------------------------------------------------------
        #   Produce map of water demand per person
        #--------------------------------------------------------------------------    
    
        withdrawal_per_capita_map = np.zeros(mapsize_global,dtype=np.single)*np.NaN
        for jj in np.arange(country_codes.shape[0]):
            country_code = country_codes.iloc[jj]['country-code']
            mask = country_code_map==country_code
            withdrawal_per_capita_map[mask] = 10**12*table_aquastat_domestic_withdrawal_interp[jj,ii]/(np.sum(pop_map[mask]*area_map[mask])+10) # liters/person/year (+10 to prevent divide by zero)

        for jj in np.arange(state_codes.shape[0]):
            state_code = state_codes.iloc[jj]['st']
            mask = state_code_map==state_code
            withdrawal_per_capita_map[mask] = 10**12*table_usgs_domestic_withdrawal_interp[jj,ii]/(np.sum(pop_map[mask]*area_map[mask])+10) # liters/person/year (+10 to prevent divide by zero)
            
        # Cap unrealistic values
        withdrawal_per_capita_map = withdrawal_per_capita_map.clip(0,5000*365.25)
                    
        # Fill gaps (Sudan for example)
        withdrawal_per_capita_map = fill(withdrawal_per_capita_map)


        #--------------------------------------------------------------------------
        #   Spatial disaggregation
        #--------------------------------------------------------------------------    
    
        data_annual_map = np.zeros(mapsize_global,dtype=np.single)
        for jj in np.arange(country_codes.shape[0]):
            country_code = country_codes.iloc[jj]['country-code']
            mask = country_code_map==country_code
            data_annual_map[mask] = np.mean(withdrawal_per_capita_map[mask])*pop_map[mask]/10**6 # mm/year
        
        for jj in np.arange(state_codes.shape[0]):
            state_code = state_codes.iloc[jj]['st']
            mask = state_code_map==state_code
            data_annual_map[mask] = np.mean(withdrawal_per_capita_map[mask])*pop_map[mask]/10**6 # mm/year

        
        #--------------------------------------------------------------------------
        #   Temporal disaggregation (Huang et al., 2018, equation 2)
        #--------------------------------------------------------------------------
        
        for month in np.arange(1,13):
            data_monthly_map = (data_annual_map/12)*(((mswx_monthly[:,:,month-1]-mswx_avg)/(mswx_max-mswx_min))*R_map+1) # mm/month
            month_ndays = monthrange(year,month)[1]
            data = data_monthly_map/month_ndays # Convert from mm/month to mm/day        
            np.savez_compressed(os.path.join(config['output_folder'],'step2_domestic_demand',str(years[ii])+str(month).zfill(2)),data=data)
        
        
        #--------------------------------------------------------------------------
        #   Plot figure
        #--------------------------------------------------------------------------
    
        # Initialize figure
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        
        # Subpanel 1
        ax1.imshow(np.sqrt(imresize_mean(data_annual_map,(360,720))),vmin=0,vmax=12,cmap=plt.get_cmap('YlGnBu'))
        ax1.set_title('Beck et al. (2022) withdrawal (mm/month)')

        # Subpanel 2
        try:
            timeindex = year-1971
            ax2.imshow(np.sqrt(np.sum(Huang_withdrawal[:,:,np.arange(timeindex*12,timeindex*12+12)],axis=2)),vmin=0,vmax=12,cmap=plt.get_cmap('YlGnBu'))
            ax2.set_title('Huang et al. (2018) withdrawal (mm/month)')
        except:
            pass
    
        # Save figure
        f.set_size_inches(10, 10)
        plt.savefig(os.path.join(config['output_folder'],'step2_domestic_demand','figures','dom_'+str(year)+'.png'),dpi=150)
        plt.close()
        
        print("Time elapsed is "+str(time.time()-t0)+" sec")

            
    ############################################################################
    #   Convert to netCDF
    ############################################################################

    varname = 'dom'
    file = os.path.join(config['output_folder'],'step2_domestic_demand',varname+'.nc')

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

    ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,200,200,), fill_value=-9999, least_significant_digit=2)
    ncfile.variables[varname].units = 'mm/d'

    for year in years:
        for month in np.arange(1,13):
            print('-------------------------------------------------------------------------------')
            print('Saving as netCDF year: '+str(year)+' month: '+str(month))
            t0 = time.time()
            
            data = np.load(os.path.join(config['output_folder'],'step2_domestic_demand',str(year)+str(month).zfill(2)+'.npz'))['data']
            data = np.round(data,decimals=2)
            data[np.isnan(data)] = 0
            data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)] # Subset to template map area        
            
            index = (year-config['year_start'])*12+month-1
             
            ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
            ncfile.variables[varname][index,:,:] = data
            
            print("Time elapsed is "+str(time.time()-t0)+" sec")
            
    ncfile.close()


if __name__ == '__main__':
    main()