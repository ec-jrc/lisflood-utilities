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
from skimage.transform import resize
from tools import *
import rasterio
from calendar import monthrange
from scipy import stats

# Load configuration file
config = load_config(sys.argv[1])

def main():
    print('===============================================================================')

    # Create output folders
    if os.path.isdir(os.path.join(config['output_folder'],'step4_livestock_demand','tables'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step4_livestock_demand','tables'))
    if os.path.isdir(os.path.join(config['output_folder'],'step4_livestock_demand','figures'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step4_livestock_demand','figures'))
       
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

    # Compute area for each grid-cell
    _, yi = np.meshgrid(np.arange(-180+template_res/2,180+template_res/2,template_res), np.arange(90-template_res/2,-90-template_res/2,-template_res))
    area_map = (40075*template_res/360)**2*np.cos(np.deg2rad(yi))

    # List of years
    years = np.arange(config['year_start'],config['year_end']+1).astype(int)

    # List of countries and US states
    country_codes = pd.read_csv(os.path.join(config['ancillary_data_folder'],'un_country_codes.csv'))
    state_codes = pd.read_csv(os.path.join(config['ancillary_data_folder'],'us-state-ansi-fips.csv'))

    # Load USGS livestock withdrawal estimates
    table_usgs_livestock_withdrawal = pd.read_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','usgs_livestock_withdrawal.csv'),index_col=0).values
    

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
    #   Load AQUASTAT livestock withdrawal data (calculated from difference
    #   between agriculture and irrigation)
    ############################################################################

    aquastat = pd.read_csv(os.path.join(config['aquastat_folder'],'aquastat_clean.csv'),index_col=False)

    table_aquastat_irrigation_withdrawal = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_aquastat_agriculture_withdrawal = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']    
        for ii in np.arange(len(years)):
            sel_ir = (aquastat['Area Id']==country_code) & (aquastat['Variable Name']=='Irrigation water withdrawal') & (aquastat['Year']==years[ii])
            if np.sum(sel_ir)>0:
                table_aquastat_irrigation_withdrawal[jj,ii] = aquastat['Value'][sel_ir].values # km3/year
            sel_ag = (aquastat['Area Id']==country_code) & (aquastat['Variable Name']=='Agricultural water withdrawal') & (aquastat['Year']==years[ii])
            if np.sum(sel_ag)>0:
                table_aquastat_agriculture_withdrawal[jj,ii] = aquastat['Value'][sel_ag].values # km3/year
    
    # Compute livestock withdrawal as difference between agriculture and irrigation
    table_aquastat_withdrawal_livestock = table_aquastat_agriculture_withdrawal-table_aquastat_irrigation_withdrawal
    table_aquastat_withdrawal_livestock[table_aquastat_withdrawal_livestock<0] = np.NaN
    
    # Save to csv
    pd.DataFrame(table_aquastat_irrigation_withdrawal,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','aquastat_irrigation_withdrawal.csv'))
    pd.DataFrame(table_aquastat_agriculture_withdrawal,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','aquastat_agriculture_withdrawal.csv'))
    pd.DataFrame(table_aquastat_withdrawal_livestock,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','aquastat_livestock_withdrawal.csv'))


    ############################################################################
    #   Load Gridded Livestock of the World (GLW) 3 data and estimate the total
    #   mass
    ############################################################################

    # Table with GLW animal types and average weight in kg (rough guesses based 
    # on quick Google searches)
    glw_table = np.array([\
        ['Bf','Buffaloes',700],\
        ['Ch','Chicken',3],\
        ['Ct','Cattle',630],\
        ['Dk','Ducks',1.5],\
        ['Gt','Goats',30],\
        ['Ho','Horses',500],\
        ['Pg','Pigs',100],\
        ['Sh','Sheep',75],\
        ],dtype=object)

    # Load data for each species
    glw_mass_maps = np.zeros((mapsize_global[0],mapsize_global[1],len(glw_table)),dtype=np.single)*np.NaN
    for ii in np.arange(len(glw_table)):
        src = rasterio.open(os.path.join(config['glw_folder'],'5_'+glw_table[ii][0]+'_2010_Da.tif'))
        data = np.array(src.read(1).astype(np.single))
        src.close()
        data[data==-np.Inf] = np.NaN
        data[np.isnan(data)] = 0
        mapsize_native = data.shape
        data = resize(data,mapsize_global)*mapsize_native[0]/mapsize_global[0] # Number of animals per grid-cell
        glw_mass_maps[:,:,ii] = data*glw_table[ii][2] # kg per grid-cell
        #print(glw_table[ii][1]+' total mass '+str(np.round(np.mean(glw_mass_maps[:,:,ii])))+' kg')

    # Compute total mass
    glw_mass_map = np.sum(glw_mass_maps,axis=2) # kg per grid-cell
    

    ############################################################################
    #   Compute total livestock mass for each country
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Computating GLW total livestock mass for each country')
    t0 = time.time()
    table_glw_livestock_mass = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']
        mask = country_code_map==country_code
        table_glw_livestock_mass[jj] = np.sum(glw_mass_map[mask])+1 # +1 to avoid divide by zero
    pd.DataFrame(table_glw_livestock_mass,index=country_codes['name'],columns=['Livestock mass [kg]']).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','glw_livestock_mass.csv'))
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    
    ############################################################################
    #   Load GCAM withdrawals and rescale for each country using GLW total 
    #   livestock mass estimates
    ############################################################################

    gcam_output = pd.read_csv(os.path.join(config['gcam_folder'],'water_demand.csv'))

    gcam_regions = np.array([\
        [['Africa_Eastern'],['Burundi', 'Comoros', 'Djibouti', 'Eritrea', 'Ethiopia', 'Kenya', 'Madagascar', 'Mauritius', 'Réunion', 'Rwanda', 'Sudan', 'Somalia', 'Uganda','South Sudan']],\
        [['Africa_Northern'],['Algeria', 'Egypt', 'Western Sahara', 'Libya', 'Morocco', 'Tunisia']],\
        [['Africa_Southern'],['Angola', 'Botswana', 'Lesotho', 'Mozambique', 'Malawi', 'Namibia', 'Eswatini', 'Tanzania, United Republic of', 'Zambia', 'Zimbabwe']],\
        [['Africa_Western'],['Benin', 'Burkina Faso', 'Central African Republic', "Côte d'Ivoire", 'Cameroon', 'Congo, Democratic Republic of the', 'Congo', 'Cabo Verde', 'Gabon', 'Ghana', 'Guinea', 'Gambia', 'Guinea-Bissau', 'Equatorial Guinea', 'Liberia', 'Mali', 'Mauritania', 'Niger', 'Nigeria', 'Senegal', 'Sierra Leone', 'Sao Tome and Principe', 'Chad', 'Togo']],\
        [['Argentina'],['Argentina']],\
        [['Australia_NZ'],['Australia', 'New Zealand']],\
        [['Brazil'],['Brazil']],\
        [['Canada'],['Canada']],\
        [['Central America and the Caribbean'],['Aruba', 'Anguilla', 'Aruba', 'Sint Maarten (Dutch part)', 'Curaçao', 'Bonaire, Sint Eustatius and Saba', 'Antigua and Barbuda', 'Bahamas', 'Belize', 'Bermuda', 'Barbados', 'Costa Rica', 'Cuba', 'Cayman Islands', 'Dominica', 'Dominican Republic', 'Guadeloupe', 'Grenada', 'Guatemala', 'Honduras', 'Haiti', 'Jamaica', 'Saint Kitts and Nevis', 'Saint Lucia', 'Montserrat', 'Martinique', 'Nicaragua', 'Panama', 'El Salvador', 'Trinidad and Tobago', 'Saint Vincent and the Grenadines', 'Puerto Rico', 'Virgin Islands (U.S.)', 'Virgin Islands (British)']],\
        [['Central Asia'],['Armenia', 'Azerbaijan', 'Georgia', 'Kazakhstan', 'Kyrgyzstan', 'Mongolia', 'Tajikistan', 'Turkmenistan', 'Uzbekistan']],\
        [['China'],['China', 'Hong Kong', 'Macao', 'Taiwan, Province of China']],\
        [['Colombia'],['Colombia']],\
        [['EU-12'],['Bulgaria', 'Cyprus', 'Czechia', 'Estonia', 'Hungary', 'Lithuania', 'Latvia', 'Malta', 'Poland', 'Romania', 'Slovakia', 'Slovenia']],\
        [['EU-15'],['Andorra', 'Austria', 'Belgium', 'Denmark', 'Finland', 'France', 'Germany', 'Greece', 'Greenland', 'Ireland', 'Italy', 'Luxembourg', 'Monaco', 'Netherlands', 'Portugal', 'Sweden', 'Spain', 'United Kingdom of Great Britain and Northern Ireland','San Marino']],\
        [['Europe_Eastern'],['Belarus', 'Moldova, Republic of', 'Ukraine']],\
        [['European Free Trade Association'],['Iceland', 'Norway', 'Switzerland']],\
        [['Europe_Non_EU'],['Albania', 'Bosnia and Herzegovina', 'Croatia', 'North Macedonia', 'Montenegro', 'Serbia', 'Turkey']],\
        [['India'],['India']],\
        [['Indonesia'],['Indonesia']],\
        [['Japan'],['Japan']],\
        [['Mexico'],['Mexico']],\
        [['Middle East'],['United Arab Emirates', 'Bahrain', 'Iran (Islamic Republic of)', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Lebanon', 'Oman', 'Palestine, State of', 'Qatar', 'Saudi Arabia', 'Syrian Arab Republic', 'Yemen']],\
        [['Pakistan'],['Pakistan']],\
        [['Russia'],['Russian Federation']],\
        [['South Africa'],['South Africa']],\
        [['South America_Northern'],['French Guiana', 'Guyana', 'Suriname', 'Venezuela (Bolivarian Republic of)']],\
        [['South America_Southern'],['Bolivia (Plurinational State of)', 'Chile', 'Ecuador', 'Peru', 'Paraguay', 'Uruguay']],\
        [['South Asia'],['Afghanistan', 'Bangladesh', 'Bhutan', 'Sri Lanka', 'Maldives', 'Nepal']],\
        [['Southeast Asia'],['American Samoa', 'Brunei Darussalam', 'Cocos (Keeling) Islands', 'Cook Islands', 'Christmas Island', 'Fiji', 'Micronesia (Federated States of)', 'Guam', 'Cambodia', 'Kiribati', "Lao People's Democratic Republic", 'Marshall Islands', 'Myanmar', 'Northern Mariana Islands', 'Malaysia', 'Mayotte', 'New Caledonia', 'Norfolk Island', 'Niue', 'Nauru', 'Northern Mariana Islands','Palau','Guam','Northern Mariana Islands', 'Pitcairn', 'Philippines', 'Palau', 'Papua New Guinea', "Korea (Democratic People's Republic of)", 'French Polynesia', 'Singapore', 'Solomon Islands', 'Seychelles', 'Thailand', 'Tokelau', 'Timor-Leste', 'Tonga', 'Tuvalu', 'Viet Nam', 'Vanuatu', 'Samoa']],\
        [['South Korea'],['Korea, Republic of']],\
        [['Taiwan'],['']],\
        [['USA'],['United States of America']],\
        ],dtype=object)

    # Values for each country and year
    table_gcam_glw_withdrawal_livestock = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        
        # GCAM withdrawals represent regional totals. To obtain country estimates,
        # we need to know the livestock mass of GCAM region and country
        
        # Which region are we in?
        region_ind = -1
        for kk in np.arange(len(gcam_regions)):
            if country_name in gcam_regions[kk][1]:
                region_ind = kk
                break
        if region_ind==-1: continue
        region_name = gcam_regions[region_ind][0][0]

        # What is the total livestock mass of the region?
        region_mass = 0
        region_countries = gcam_regions[region_ind][1]
        for region_country in region_countries:
            hh = country_codes['name'].tolist().index(region_country)
            region_mass += table_glw_livestock_mass[hh]
        
        # Aggregate different sectors
        sel = ["Reference" in s for s in gcam_output['scenario']] &\
            (gcam_output['region']==region_name) &\
            (gcam_output['sector'].isin(['SheepGoat','Beef','Dairy','Pork','Poultry']))
        gcam_livestock_withdrawal = np.sum(gcam_output.iloc[np.where(sel)[0],3:-1].values,axis=0) # km3/year
        gcam_years = np.array([int(x) for x in gcam_output.columns[3:-1].values])    
        
        # Rescale regional withdrawals based on country population
        mass_frac = table_glw_livestock_mass[jj]/region_mass
        for kk in np.arange(len(gcam_years)):
            try:
                ii = years.tolist().index(gcam_years[kk])
            except:
                continue
            table_gcam_glw_withdrawal_livestock[jj,ii] = gcam_livestock_withdrawal[kk]*mass_frac
        
    # Save to csv
    pd.DataFrame(table_gcam_glw_withdrawal_livestock,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','gcam_glw_withdrawal_livestock.csv'))


    ############################################################################
    #   Verification of estimates. Huang et al. (2018) is supposed to represent 
    #   difference between AQUASTAT agriculture and irrigation, which is the 
    #   same approach as we use, but there are large differences.
    ############################################################################

    # Load Huang et al. (2018) livestock estimates
    table_Huang_withdrawal_livestock = pd.read_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','Huang_withdrawal_livestock.csv'),index_col=0).values

    # What are the estimates around 2010 for each data source?
    ii_2010 = years.tolist().index(2010)
    aquastat = table_aquastat_withdrawal_livestock[:,ii_2010+2]
    gcam = table_gcam_glw_withdrawal_livestock[:,ii_2010]
    huang = table_Huang_withdrawal_livestock[:,ii_2010]
    sel = np.isnan(aquastat+gcam+huang)==False
    
    # Print estimates for a few major countries
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        if aquastat[jj]>0.2:
            print('-------------------------------------------------------------------------------')        
            print(country_name)
            print('aquastat '+str(np.round(aquastat[jj],2)))
            print('gcam '+str(np.round(gcam[jj],2)))
            print('huang '+str(np.round(huang[jj],2)))
    
    # How are the estimates correlated?
    stats.spearmanr(aquastat[sel], gcam[sel])
    stats.spearmanr(aquastat[sel], huang[sel])
    stats.spearmanr(gcam[sel], huang[sel])
    
    # What do the different distributions look like?
    plt.figure(figsize=(8,6))
    plt.hist(aquastat[sel], bins=300, alpha=0.5, range=(0,15), label="aquastat")
    plt.hist(gcam[sel], bins=300, alpha=0.5, range=(0,15), label="gcam")
    plt.hist(huang[sel], bins=300, alpha=0.5, range=(0,15), label="huang")
    plt.legend()
    plt.show(block=False)    
        
    # Can we derive a simple correction factor to improve the GCAM+GLW estimates? Not really because mean and median give opposite results...
    np.mean(aquastat[sel])/np.mean(gcam[sel])
    np.median(aquastat[sel])/np.median(gcam[sel])
    

    ############################################################################
    #   Load and reproject Huang et al. (2018) water demand data (only used for
    #   comparison)
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and reprojecting Huang et al. (2018) water demand data')
    t0 = time.time()

    # Loop over vars
    vars = ['liv']
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
    #   Loop through countries and compute final livestock withdrawal time 
    #   series. We will use AQUASTAT just for rescaling, not to derive time 
    #   series, due to the uncertainty in the AQUASTAT estimates (which 
    #   represent the difference between agriculture and irrigation estimates).
    ############################################################################

    table_livestock_industry = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']
        country_name = country_codes.iloc[jj]['name']
        country_acronym = country_codes.iloc[jj]['alpha-3']    
        print('jj='+str(jj)+' '+country_name)
        
        # Linearly interpolate and nearest-neighbor extrapolate GCAM+GLW time series
        ts_gcam = pd.Series(table_gcam_glw_withdrawal_livestock[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values
        
        # If AQUASTAT estimate is available, rescale GCAM+GLW time series (using median 
        # rescaling factor to reduce uncertainty)
        factors = (table_aquastat_withdrawal_livestock[jj,:]+10**-6)/(ts_gcam+10**-6)
        factor = np.nanmedian(factors)
        if np.isnan(factor)==False:
            ts_gcam = ts_gcam*factor
        
        # If AQUASTAT estimate is not available, use GCAM+GLW estimate
        table_livestock_industry[jj,:] = ts_gcam

    pd.DataFrame(table_livestock_industry,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step4_livestock_demand','tables','livestock_industry.csv'))        
        
        
    ############################################################################
    #   Downscale livestock withdrawals
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Downscaling livestock withdrawals')
    t0 = time.time()
    
    varname = 'liv'
    
    #--------------------------------------------------------------------------
    #   Initialize netCDF 
    #--------------------------------------------------------------------------

    file = os.path.join(config['output_folder'],'step4_livestock_demand',varname+'.nc')

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

    ncfile.createVariable(varname, np.single, ('time', 'lat', 'lon'), zlib=True, chunksizes=(1,32,32,), fill_value=-9999, least_significant_digit=1)
    ncfile.variables[varname].units = 'mm/d'

    for ii in np.arange(len(years)):    
        year = years[ii]
        print(varname+' year '+str(year))
        t1 = time.time()
                
        
        #--------------------------------------------------------------------------
        #   Spatial downscaling
        #--------------------------------------------------------------------------
    
        # Spatially downscale annual country withdrawals using GLW total livestock mass grid
        data_annual_map = np.zeros(mapsize_global,dtype=np.single)*np.NaN
        for jj in np.arange(country_codes.shape[0]):
            country_code = country_codes.iloc[jj]['country-code']
            mask = country_code_map==country_code
            country_val = table_livestock_industry[jj,ii] # Country withdrawal estimate
            data_annual_map[mask] = 10**6*glw_mass_map[mask]*country_val/np.sum(glw_mass_map[mask]*area_map[mask])
            #np.sum(area_map[mask]*data_annual_map[mask]/1000/1000)            
            
        # Spatially downscale annual US state withdrawals using population data
        for jj in np.arange(state_codes.shape[0]):
            state_code = state_codes.iloc[jj]['st']
            mask = state_code_map==state_code
            state_val = pd.Series(table_usgs_livestock_withdrawal[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values[ii] # State withdrawal estimate from preceding script
            data_annual_map[mask] = 10**6*glw_mass_map[mask]*state_val/np.sum(glw_mass_map[mask]*area_map[mask])
            #np.sum(area_map[mask]*data_annual_map[mask]/1000/1000)
            
        # Check if there are too many missing values
        nan_percentage = 100*np.sum(np.isnan(data_annual_map))/(data_annual_map.shape[0]*data_annual_map.shape[1])
        assert nan_percentage<37
            
        # Replace NaNs with zeros
        data_annual_map[np.isnan(data_annual_map)] = 0
        
        
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
            varindex = 0
            timeindex = year-1971
            ax2.imshow(np.sqrt(np.sum(Huang_withdrawal[varindex,:,:,np.arange(timeindex*12,timeindex*12+12)],axis=0)),vmin=0,vmax=15,cmap=plt.get_cmap('YlGnBu'))
            ax2.set_title('Huang et al. (2018) withdrawal (mm/month)')
        except:
            pass
    
        # Save figure
        f.set_size_inches(10, 10)
        plt.savefig(os.path.join(config['output_folder'],'step4_livestock_demand','figures',varname+'_'+str(year)+'.png'),dpi=150)
        plt.close()
            
        
        #--------------------------------------------------------------------------
        #   Save to netCDF
        #--------------------------------------------------------------------------
    
        for month in np.arange(1,13):
            month_ndays = monthrange(year,month)[1]
            data = data_annual_map/365.25
            data = np.round(data,decimals=2)
            data = data[row_upper:row_upper+len(template_lat),col_left:col_left+len(template_lon)] # Subset to template map area
            index = (year-config['year_start'])*12+month-1
            ncfile.variables['time'][index] = (pd.to_datetime(datetime(year,month,1))-pd.to_datetime(datetime(1979, 1, 1))).total_seconds()/86400    
            ncfile.variables[varname][index,:,:] = data
               
        print("Time elapsed is "+str(time.time()-t1)+" sec")
            
    ncfile.close()

    print("Total time elapsed is "+str(time.time()-t0)+" sec")

    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()