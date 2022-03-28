#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "January 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from skimage.transform import resize
from tools import *
import rasterio
from calendar import monthrange

# Load configuration file
config = load_config(sys.argv[1])

def main():
    print('===============================================================================')

    # Create output folder
    if os.path.isdir(os.path.join(config['output_folder'],'step3b_industrial_demand','figures'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step3b_industrial_demand','figures'))
       
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
    #   Load country border raster
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and resampling country border raster')
    t0 = time.time()
    country_code_map = load_country_code_map(os.path.join(config['world_borders_folder'],'TM_WORLD_BORDERS_UN_rasterized.tif'),mapsize_global)
    country_code_map = fill(country_code_map)
    country_code_map_1800x3600 = resize(country_code_map,(1800,3600),order=0,mode='edge',anti_aliasing=False)
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    

    ############################################################################
    #   Load US states raster
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and resampling US state border raster')
    t0 = time.time()
    state_code_map = load_us_state_code_map(os.path.join(config['us_states_folder'],'cb_2018_us_state_500k_rasterized.tif'),mapsize_global)
    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load and reproject Huang et al. (2018) water demand data (only used for
    #   comparison)
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and reprojecting Huang et al. (2018) water demand data')
    t0 = time.time()

    # Loop over vars
    vars = ['elec','mfg']
    Huang_withdrawal = np.zeros((len(vars),360,720,480),dtype=np.single)*np.NaN
    for vv in np.arange(len(vars)):
        var = vars[vv]
        
        # Load raw data
        dset = Dataset(os.path.join(config['huang_folder'],'withd_'+var+'.nc'))
        raw_data = np.array(dset.variables['withd_'+var][:])
        raw_lat = np.array(dset.variables['lat'][:])
        raw_lon = np.array(dset.variables['lon'][:])
        dset.close()

        # Reproject data
        rows,cols = latlon2rowcol(raw_lat,raw_lon,0.5,90,-180)
        for ii in np.arange(raw_data.shape[0]):    
            reprojected = np.zeros((360,720),dtype=np.single)
            reprojected[rows,cols] = raw_data[ii,:]
            Huang_withdrawal[vv,:,:,ii] = reprojected

    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load all tables from preceding script
    ############################################################################

    files1 = glob.glob(os.path.join(config['output_folder'],'step2_domestic_demand','tables','*.csv'))
    files2 = glob.glob(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','*.csv'))
    files = files1+files2
    tables = {}
    for ii in np.arange(len(files)):
        fname = os.path.basename(files[ii])[:-4]
        tables[fname] = pd.read_csv(files[ii],index_col=0).values

    files = glob.glob(os.path.join(config['output_folder'],'step2_domestic_demand','tables','*.csv'))
    for ii in np.arange(len(files)):
        fname = os.path.basename(files[ii])[:-4]
        tables[fname] = pd.read_csv(files[ii],index_col=0).values


    ############################################################################
    #   Downscale manufacturing and thermoelectric withdrawals
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Downscaling manufacturing and thermoelectric withdrawals')
    t0 = time.time()

    varnames = ['ind','ene']

    for varname in varnames:

        # Load country withdrawals from preceding script
        if varname=='ene':
            country_table = tables['withdrawal_thermoelectric']
            state_table = tables['usgs_thermoelectric_withdrawal']
        if varname=='ind':
            country_table = tables['withdrawal_manufacturing']
            state_table = tables['usgs_manufacturing_withdrawal']
        
        
        #--------------------------------------------------------------------------
        #   Initialize netCDF 
        #--------------------------------------------------------------------------

        file = os.path.join(config['output_folder'],'step3b_industrial_demand',varname+'.nc')

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

        ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,200,200,), fill_value=-9999, least_significant_digit=1)
        ncfile.variables[varname].units = 'mm/d'

        for ii in np.arange(len(years)):    
            year = years[ii]
            print(varname+' year '+str(year))
            t1 = time.time()
                    
            # Load population data
            npz = np.load(os.path.join(config['output_folder'],'step1_population_density',str(year)+'.npz'))
            pop_map = npz['data']+10**-6

            
            #--------------------------------------------------------------------------
            #   Spatial downscaling
            #--------------------------------------------------------------------------
        
            # Spatially downscale annual country withdrawals using population data
            data_annual_map = np.zeros(mapsize_global,dtype=np.single)*np.NaN
            for jj in np.arange(country_codes.shape[0]):
                country_code = country_codes.iloc[jj]['country-code']
                mask = country_code_map==country_code
                country_val = country_table[jj,ii] # Country withdrawal estimate from preceding script
                data_annual_map[mask] = 10**6*pop_map[mask]*country_val/(np.sum(pop_map[mask]*area_map[mask]))
                #np.sum(area_map[mask]*data_annual_map[mask]/1000/1000)

            # Spatially downscale annual US state withdrawals using population data
            for jj in np.arange(state_codes.shape[0]):
                state_code = state_codes.iloc[jj]['st']
                mask = state_code_map==state_code
                state_val = pd.Series(state_table[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values[ii] # State withdrawal estimate from preceding script
                data_annual_map[mask] = 10**6*pop_map[mask]*state_val/(np.sum(pop_map[mask]*area_map[mask]))
                #np.sum(area_map[mask]*data_annual_map[mask]/1000/1000)
                    
            # Check if there are too many missing values
            nan_percentage = 100*np.sum(np.isnan(data_annual_map))/(data_annual_map.shape[0]*data_annual_map.shape[1])
            assert nan_percentage<2
                
            # Replace NaNs with zeros
            data_annual_map[np.isnan(data_annual_map)] = 0
            
            
            #--------------------------------------------------------------------------
            #   Temporal downscaling
            #--------------------------------------------------------------------------
        
            # Compute p coefficients for temporally downscaling thermoelectric (see Huang et al., 2018)        
            p_b, p_it = np.zeros(mapsize_global,dtype=np.single)*np.NaN,np.zeros(mapsize_global,dtype=np.single)*np.NaN
            p_h, p_c, p_u = np.zeros(mapsize_global,dtype=np.single)*np.NaN,np.zeros(mapsize_global,dtype=np.single)*np.NaN,np.zeros(mapsize_global,dtype=np.single)*np.NaN
            for jj in np.arange(country_codes.shape[0]):
                country_code = country_codes.iloc[jj]['country-code']
                country_name = country_codes.iloc[jj]['name']
                country_acronym = country_codes.iloc[jj]['alpha-3']
                kw = dict(method="linear",fill_value="extrapolate", limit_direction="both")
                gcam_elec_building = pd.Series(tables['gcam_elec_building'][jj,:]).interpolate(**kw).values[ii]+10**-6 # +10**-6 to avoid divide by zero
                gcam_elec_trans_ind = pd.Series(tables['gcam_elec_trans_ind'][jj,:]).interpolate(**kw).values[ii]+10**-6
                gcam_elec_heating = pd.Series(tables['gcam_elec_heating'][jj,:]).interpolate(**kw).values[ii]+10**-6
                gcam_elec_cooling = pd.Series(tables['gcam_elec_cooling'][jj,:]).interpolate(**kw).values[ii]+10**-6
                gcam_elec_other = pd.Series(tables['gcam_elec_other'][jj,:]).interpolate(**kw).values[ii]+10**-6
                mask = country_code_map==country_code    
                p_b[mask] = gcam_elec_building/(gcam_elec_building+gcam_elec_trans_ind)
                p_it[mask] = gcam_elec_trans_ind/(gcam_elec_building+gcam_elec_trans_ind)
                p_h[mask] = gcam_elec_heating/(gcam_elec_heating+gcam_elec_cooling+gcam_elec_other)
                p_c[mask] = gcam_elec_cooling/(gcam_elec_heating+gcam_elec_cooling+gcam_elec_other)
                p_u[mask] = gcam_elec_other/(gcam_elec_heating+gcam_elec_cooling+gcam_elec_other)
            p_b, p_it = fill(p_b), fill(p_it)
            p_h, p_c, p_u = fill(p_h), fill(p_c), fill(p_u)
            
            # Load HDD and CDD data
            hdd_monthly_maps = np.zeros((mapsize_global[0],mapsize_global[1],12),dtype=np.single)*np.NaN
            cdd_monthly_maps = np.zeros((mapsize_global[0],mapsize_global[1],12),dtype=np.single)*np.NaN
            for month in np.arange(1,13):
                npz = np.load(os.path.join(config['output_folder'],'step3a_industrial_demand','hdd',str(year)+str(month).zfill(2)+'.npz'))
                hdd_monthly_maps[:,:,month-1] = imresize_mean(npz['data'],mapsize_global)+10**-6
                npz = np.load(os.path.join(config['output_folder'],'step3a_industrial_demand','cdd',str(year)+str(month).zfill(2)+'.npz'))
                cdd_monthly_maps[:,:,month-1] = imresize_mean(npz['data'],mapsize_global)+10**-6

            # Temporally downscale annual withdrawals
            data_monthly_maps = np.zeros((mapsize_global[0],mapsize_global[1],12),dtype=np.single)*np.NaN
            for month in np.arange(1,13):
                month_ndays = monthrange(year,month)[1]
                if varname=='ene':
                    
                    # Compute annual HDD and CDD sums
                    hdd_annual_sum = np.sum(hdd_monthly_maps,axis=2)
                    cdd_annual_sum = np.sum(cdd_monthly_maps,axis=2)
                    
                    # Temporally downscale withdrawals following Huang et al. (2018) equations 7 to 10
                    data_monthly_map = np.zeros(mapsize_global,dtype=np.single)*np.NaN
                    mask = (hdd_annual_sum<650) & (cdd_annual_sum<450)
                    data_monthly_map[mask] = month_ndays*data_annual_map[mask]/365.25                
                    mask = (hdd_annual_sum>=650) & (cdd_annual_sum<450)
                    data_monthly_map[mask] = data_annual_map[mask]*(p_b[mask]*((p_h[mask]+p_c[mask])*hdd_monthly_maps[:,:,month-1][mask]/hdd_annual_sum[mask]+p_u[mask]*month_ndays/365.25)+p_it[mask]*month_ndays/365.25)
                    mask = (hdd_annual_sum<650) & (cdd_annual_sum>=450)
                    data_monthly_map[mask] = data_annual_map[mask]*(p_b[mask]*((p_h[mask]+p_c[mask])*cdd_monthly_maps[:,:,month-1][mask]/cdd_annual_sum[mask]+p_u[mask]*month_ndays/365.25)+p_it[mask]*month_ndays/365.25)
                    mask = (hdd_annual_sum>=650) & (cdd_annual_sum>=450)
                    data_monthly_map[mask] = data_annual_map[mask]*(p_b[mask]*(p_h[mask]*hdd_monthly_maps[:,:,month-1][mask]/hdd_annual_sum[mask]+p_c[mask]*cdd_monthly_maps[:,:,month-1][mask]/cdd_annual_sum[mask]+p_u[mask]*month_ndays/365.25)+p_it[mask]*month_ndays/365.25)
                    data_monthly_map[np.isnan(data_monthly_map)] = 0
                    data_monthly_maps[:,:,month-1] = data_monthly_map
                    
                if varname=='ind':
                    data_monthly_maps[:,:,month-1] = month_ndays*data_annual_map/365.25
           
            
            #--------------------------------------------------------------------------
            #   Plot figure
            #--------------------------------------------------------------------------
        
            # Initialize figure
            f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
            
            # Subpanel 1
            ax1.imshow(np.sqrt(imresize_mean(data_annual_map,(360,720))),vmin=0,vmax=15,cmap=plt.get_cmap('YlGnBu'))
            ax1.set_title('Beck et al. (2022) withdrawal (mm/month)')

            # Subpanel 2
            try:
                if varname=='ene': varindex = 0
                if varname=='ind': varindex = 1
                timeindex = year-1971
                ax2.imshow(np.sqrt(np.sum(Huang_withdrawal[varindex,:,:,np.arange(timeindex*12,timeindex*12+12)],axis=0)),vmin=0,vmax=15,cmap=plt.get_cmap('YlGnBu'))
                ax2.set_title('Huang et al. (2018) withdrawal (mm/month)')
            except:
                pass
        
            # Save figure
            f.set_size_inches(10, 10)
            plt.savefig(os.path.join(config['output_folder'],'step3b_industrial_demand','figures',varname+'_'+str(year)+'.png'),dpi=150)
            plt.close()
                
            
            #--------------------------------------------------------------------------
            #   Save to netCDF
            #--------------------------------------------------------------------------
        
            for month in np.arange(1,13):
                month_ndays = monthrange(year,month)[1]
                data = data_monthly_maps[:,:,month-1]/month_ndays
                data = np.round(data,decimals=2)
                data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)] # Subset to template map area
                index = (year-config['year_start'])*12+month-1
                ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
                ncfile.variables[varname][index,:,:] = data
                   
            print("Time elapsed is "+str(time.time()-t1)+" sec")
                
        ncfile.close()

    print("Total time elapsed is "+str(time.time()-t0)+" sec")


    '''
    ############################################################################
    #   Load HTAP data
    ############################################################################

    dset = Dataset(os.path.join(config['htap_folder'],'htapv2.2.emisso2.surface.x3600y1800t12.2016.integrate.nc4'))
    raw = np.array(dset.variables['sanl1'][:])
    mean = np.flipud(np.nanmean(raw,axis=0)*10**9)
    lat = np.array(dset.variables['lat'][:])[::-1]
    lon = np.array(dset.variables['lon'][:])
    dset.close()

    file = os.path.join(config['htap_folder'],'htapv2.2.emisso2.surface.x3600y1800t12.2016.integrate_mean_hylke.nc')

    ncfile = Dataset(file, 'w', format='NETCDF4')
    ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

    ncfile.createDimension('lon', len(lon))
    ncfile.createDimension('lat', len(lat))

    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'

    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'][:] = lat
    ncfile.variables['lat'].long_name = 'latitude'

    ncfile.createVariable(varname, np.single, ('lat', 'lon'), zlib=True, chunksizes=(32,32,), fill_value=-9999) #, least_significant_digit=9

    ncfile.variables[varname][:,:] = mean
           
    ncfile.close()

    '''


if __name__ == '__main__':
    main()