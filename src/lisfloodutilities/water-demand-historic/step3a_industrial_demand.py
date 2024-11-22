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
    if os.path.isdir(os.path.join(config['output_folder'],'step3a_industrial_demand','figures'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step3a_industrial_demand','figures'))
    if os.path.isdir(os.path.join(config['output_folder'],'step3a_industrial_demand','tables'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step3a_industrial_demand','tables'))
    if os.path.isdir(os.path.join(config['output_folder'],'step3a_industrial_demand','cdd'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step3a_industrial_demand','cdd'))
    if os.path.isdir(os.path.join(config['output_folder'],'step3a_industrial_demand','hdd'))==False:
        os.makedirs(os.path.join(config['output_folder'],'step3a_industrial_demand','hdd'))

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

    # List of years
    years = np.arange(config['year_start'],config['year_end']+1).astype(int)
    ii_2000 = years.tolist().index(2000)

    # List of countries
    country_codes = pd.read_csv(os.path.join(config['ancillary_data_folder'],'un_country_codes.csv'))

    # Load country population estimates from preceding script
    table_population = pd.read_csv(os.path.join(config['output_folder'],'step2_domestic_demand','tables','population.csv'),index_col=0).values


    ############################################################################
    #   Load country border raster
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and resampling country border raster')
    t0 = time.time()
    country_code_map = load_country_code_map(os.path.join(config['world_borders_folder'],'CNTR_RG_01M_2024_4326_rasterized.tif'),mapsize_global)
    country_code_map_360x720 = resize(country_code_map,(360,720),order=0,mode='edge',anti_aliasing=False)
    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Load and reproject Huang et al. (2018) water demand data (only used for
    #   comparison)
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Loading and reprojecting Huang et al. (2018) water demand data')
    t0 = time.time()

    # Spatiotemporal information
    Huang_dates = pd.date_range(start='1/1/1971', end='12/1/2010', freq='MS')  
    Huang_res = 0.5
    _, yi = np.meshgrid(np.arange(-180+Huang_res/2,180+Huang_res/2,Huang_res), np.arange(90-Huang_res/2,-90-Huang_res/2,-Huang_res))
    area_map_small = (40075*Huang_res/360)**2*np.cos(np.deg2rad(yi))

    # Loop over vars
    vars = ['elec','mfg','liv']
    table_Huang_withdrawal = np.zeros((len(vars),country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for vv in np.arange(len(vars)):
        var = vars[vv]
        
        # Load raw data
        dset = Dataset(os.path.join(config['huang_folder'],'withd_'+var+'.nc'))
        raw_data = np.array(dset.variables['withd_'+var][:])
        raw_lat = np.array(dset.variables['lat'][:])
        raw_lon = np.array(dset.variables['lon'][:])
        dset.close()

        # Reproject data
        Huang_withdrawal = np.zeros((360,720,raw_data.shape[0]),dtype=np.single)*np.NaN
        rows,cols = latlon2rowcol(raw_lat,raw_lon,0.5,90,-180)
        for ii in np.arange(raw_data.shape[0]):    
            reprojected = np.zeros((360,720),dtype=np.single)*np.NaN
            reprojected[rows,cols] = raw_data[ii,:]
            Huang_withdrawal[:,:,ii] = reprojected

        # Compute totals for each year and country
        for ii in np.arange(len(years)):
            sel = Huang_dates.year==years[ii]
            if np.sum(sel)==0: continue
            annual_sum = np.sum(Huang_withdrawal[:,:,sel],axis=2)
            for jj in np.arange(country_codes.shape[0]):        
                country_code = country_codes.iloc[jj]['country-code']
                mask = country_code_map_360x720==country_code
                table_Huang_withdrawal[vv,jj,ii] = np.nansum(annual_sum[mask]*area_map_small[mask]/10**6) # km3/year

    # Produce tables
    table_Huang_withdrawal_thermoelectric = table_Huang_withdrawal[0,:,:]
    table_Huang_withdrawal_manufacturing = table_Huang_withdrawal[1,:,:]    
    table_Huang_withdrawal_industry = table_Huang_withdrawal_thermoelectric+table_Huang_withdrawal_manufacturing
    table_Huang_withdrawal_livestock = table_Huang_withdrawal[2,:,:]

    # Save to csv
    pd.DataFrame(table_Huang_withdrawal_thermoelectric,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','Huang_withdrawal_thermoelectric.csv'))
    pd.DataFrame(table_Huang_withdrawal_manufacturing,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','Huang_withdrawal_manufacturing.csv'))
    pd.DataFrame(table_Huang_withdrawal_livestock,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','Huang_withdrawal_livestock.csv'))

    print("Time elapsed is "+str(time.time()-t0)+" sec")


    '''
    ############################################################################
    #   Obtain list of closest neigbors for each country
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Obtaining list of closest neigbors for each country')
    t0 = time.time()

    # Determine centroid coordinates for each country
    country_centroid_lats, country_centroid_lons = np.zeros((country_codes.shape[0],))*np.NaN, np.zeros((country_codes.shape[0],))*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']
        mask = country_code_map==country_code
        if np.sum(mask)==0: continue
        ind = np.where(mask)
        country_centroid_lats[jj], country_centroid_lons[jj] = 90-180*np.mean(ind[0])/mapsize_global[0]-template_res/2, -180+360*np.mean(ind[1])/mapsize_global[1]+template_res/2

    # Determine closest neighbors for each country
    from geopy.distance import great_circle
    closest_codes = np.zeros((country_codes.shape[0],country_codes.shape[0]))*np.NaN
    closest_indices = np.zeros((country_codes.shape[0],country_codes.shape[0]),dtype=int)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        if np.isnan(country_centroid_lats[jj]): continue
        distances = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
        for kk in np.arange(country_codes.shape[0]):
            if np.isnan(country_centroid_lats[kk]) or (jj==kk): continue
            distances[kk] = great_circle((country_centroid_lats[jj],country_centroid_lons[jj]),(country_centroid_lats[kk],country_centroid_lons[kk])).meters/1000
        indices = np.argsort(distances)
        closest_indices[jj,:] = indices
        closest_codes[jj,:] = country_codes['country-code'][indices]
            
    # Verification: neighbor of US (code 840) should be Canada (code 124)
    index = country_codes['country-code'].tolist().index(840)
    if closest_codes[index,0]!=124:
        raise ValueError('Something wrong')

    print("Time elapsed is "+str(time.time()-t0)+" sec")
    '''

    ############################################################################
    #   Load Vassolo and Doll (2005) data
    ############################################################################

    src = rasterio.open(os.path.join(config['vassolo_doll_folder'],'manuf_wwd.tif')) # Million m3
    manuf_wwd = src.read(1).astype(np.single)
    src.close()

    src = rasterio.open(os.path.join(config['vassolo_doll_folder'],'wwd_ps.tif')) # Million m3
    wwd_ps = src.read(1).astype(np.single)
    src.close()

    vassolo_doll_withdrawal_manufacturing = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
    vassolo_doll_withdrawal_thermoelectric = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']
        mask = country_code_map_360x720==country_code
        vassolo_doll_withdrawal_manufacturing[jj] = np.nansum(manuf_wwd[mask])/1000 # km3
        vassolo_doll_withdrawal_thermoelectric[jj] = np.nansum(wwd_ps[mask])/1000 # km3

    # Save to csv
    pd.DataFrame(vassolo_doll_withdrawal_manufacturing,index=country_codes['name'],columns=['Withdrawal 1995 (km3/year)']).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','vassolo_doll_withdrawal_manufacturing.csv'))
    pd.DataFrame(vassolo_doll_withdrawal_thermoelectric,index=country_codes['name'],columns=['Withdrawal 1995 (km3/year)']).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','vassolo_doll_withdrawal_thermoelectric.csv'))

    '''
    ############################################################################
    #   Load Lohrmann et al. (2019) power plant data
    ############################################################################

    # Load raw data
    lohrmann_plant_withdrawals = pd.read_excel(os.path.join(config['lohrmann_power_plant_db_folder'],'41560_2019_501_MOESM2_ESM.xlsx'),sheet_name='Lifetime Scenario')
    lohrmann_details = pd.read_excel(os.path.join(config['lohrmann_power_plant_db_folder'],'41560_2019_501_MOESM2_ESM.xlsx'),sheet_name='List of power plants')
    lohrmann = pd.merge(lohrmann_plant_withdrawals, lohrmann_details)
    del lohrmann_plant_withdrawals
    del lohrmann_details

    # Replace Lohrmann country names with UN country names
    rename_list = np.array([\
        ["Bolivia","Bolivia (Plurinational State of)"],\
        ["Congo Republic","Congo, Democratic Republic of the"],\
        ["Czech Republic","Czechia"],\
        ["Hong Kong Region of China","Hong Kong"],\
        ["Iran","Iran (Islamic Republic of)"],\
        ["Ivory Coast","Côte d'Ivoire"],\
        ["Korea Democratic People's Republic","Korea (Democratic People's Republic of)"],\
        ["Korea Republic","Korea, Republic of"],\
        ["Kosovo","Serbia"],\
        ["Macau Region of China","Macau"],\
        ["Macedonia Republic","North Macedonia"],\
        ["Republic of Moldova","Moldova, Republic of"],\
        ["Reunion","Réunion"],\
        ["Sint Maarten","Sint Maarten (Dutch part)"],\
        ["Syria","Syrian Arab Republic"],\
        ["Taiwan","Taiwan, Province of China"],\
        ["Tanzania","Tanzania, United Republic of"],\
        ["United Kingdom","United Kingdom of Great Britain and Northern Ireland"],\
        ["United States","United States of America"],\
        ["Venezuela","Venezuela (Bolivarian Republic of)"],\
        ["Vietnam","Viet Nam"],\
        ])

    lohrmann['Country'] = lohrmann['Country'].replace(rename_list[:,0],rename_list[:,1])

    # Compute aggregate capacities and withdrawals for each country
    table_lohrmann_plant_capacity = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
    table_lohrmann_plant_withdrawal = np.zeros((country_codes.shape[0],),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        sel = (lohrmann['Country'].values==country_name) & (lohrmann['Seawater-cooling'].values==False)
        table_lohrmann_plant_capacity[jj] = np.sum(lohrmann['Active Capacity [MW]'].values[sel]) # MW
        table_lohrmann_plant_withdrawal[jj] = np.sum(lohrmann['Withdrawal 2015 [m3]'].values[sel])/10**9 # km3/year

    # Save to csv
    pd.DataFrame(table_lohrmann_plant_capacity,index=country_codes['name'],columns=['Capacity (MW)']).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','lohrmann_plant_capacity.csv'))
    pd.DataFrame(table_lohrmann_plant_withdrawal,index=country_codes['name'],columns=['Withdrawal 2015 (km3/year)']).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','lohrmann_plant_withdrawal.csv'))
    '''

    '''
    ############################################################################
    #   Load EIA electricity installed capacity data
    ############################################################################

    # Load raw data
    eia_file = glob.glob(os.path.join(config['eia_folder'],'*.csv'))[0]
    eia_db = pd.read_csv(eia_file,skiprows=1)
    eia_db.replace("--", "", regex=True, inplace=True)
    eia_db.replace(r'^\s*$', np.nan,regex=True, inplace=True) # Turn empty strings into NaN
    eia_db.iloc[:,1] = [elem.replace('        ', '') for elem in eia_db.iloc[:,1]] # Remove leading spaces in country names

    # Add electricity source column
    eia_db['source'] = np.NaN
    for mm in np.arange(eia_db.shape[0]):
        if 'electricity' in eia_db.iloc[mm,1]:
            source = eia_db.iloc[mm,1]
        eia_db.loc[mm,'source'] = source
       
    # Country names should be changed to match UN country names
    rename_list = np.array([\
        [["Gambia, The"],["Gambia"]],\
        [["Former Serbia and Montenegro"],["Montenegro","Serbia"]],\
        [["Czech Republic"],["Czechia"]],\
        [["Former U.S.S.R."],["Russian Federation"]],\
        [["Saint Vincent/Grenadines"],["Saint Vincent and the Grenadines"]],\
        [["Taiwan"],["Taiwan, Province of China"]],\
        [["North Korea"],["Korea (Democratic People's Republic of)"]],\
        [["Iran"],["Iran (Islamic Republic of)"]],\
        [["Moldova"],["Moldova, Republic of"]],\
        [["Congo-Kinshasa"],["Congo, Democratic Republic of the"]],\
        [["Reunion"],["Réunion"]],\
        [["Germany, West","Germany, East"],["Germany"]],\
        [["South Korea"],["Korea, Republic of"]],\
        [["U.S. Virgin Islands"],["Virgin Islands (U.S.)"]],\
        [["U.S. Territories"],["Samoa","Guam","Northern Mariana Islands","Puerto Rico","Virgin Islands (U.S.)"]],\
        [["Tanzania"],["Tanzania, United Republic of"]],\
        [["Côte d’Ivoire"],["Côte d'Ivoire"]],\
        [["Micronesia"],["Micronesia (Federated States of)"]],\
        [["Kosovo"],["?"]],\
        [["United States"],["United States of America"]],\
        [["Vietnam"],["Viet Nam"]],\
        [["The Bahamas"],["Bahamas"]],\
        [["Macau"],["Macao"]],\
        [["British Virgin Islands"],["Virgin Islands (British)"]],\
        [["Falkland Islands"],["Falkland Islands (Malvinas)"]],\
        [["Netherlands Antilles"],["Curaçao","Sint Maarten (Dutch part)","Bonaire, Sint Eustatius and Saba","Aruba"]],\
        [["Former Yugoslavia"],["Serbia","Croatia","Bosnia and Herzegovina","North Macedonia","Montenegro","Slovenia"]],\
        [["Syria"],["Syrian Arab Republic"]],\
        [["Former Czechoslovakia"],["Czechia","Slovakia"]],\
        [["Venezuela"],["Venezuela (Bolivarian Republic of)"]],\
        [["Bolivia"],["Bolivia (Plurinational State of)"]],\
        [["Hawaiian Trade Zone"],["?"]],\
        [["U.S. Pacific Islands"],["Guam","Samoa","Northern Mariana Islands","Palau","Micronesia (Federated States of)","Marshall Islands"]],\
        [["Brunei"],["Brunei Darussalam"]],\
        [["Burma"],["Myanmar"]],\
        [["Congo-Brazzaville"],["Congo"]],\
        [["Saint Helena"],["Saint Helena, Ascension and Tristan da Cunha"]],\
        [["United Kingdom"],["United Kingdom of Great Britain and Northern Ireland"]],\
        [["Wake Island"],["?"]],\
        [["Russia"],["Russian Federation"]],\
        [["Palestinian Territories"],["Palestine, State of"]],\
        [["Laos"],["Lao People's Democratic Republic"]],\
        ],dtype=object)

    # Rename EIA country names (skips countries without one-on-one match)
    for mm in np.arange(eia_db.shape[0]):
        for nn in np.arange(len(rename_list)):
            if (len(rename_list[nn][0])==1) & (len(rename_list[nn][1])==1):
                if eia_db.iloc[mm,1]==rename_list[nn][0][0]:
                    eia_db.iloc[mm,1] = rename_list[nn][1][0]

    # Types of capacity/generation to include
    capacity_source_keywords = ["nuclear electricity installed capacity","fossil fuels electricity installed capacity","geothermal electricity installed capacity","biomass and waste electricity installed capacity"]
    generation_source_keywords = ["nuclear electricity net generation","fossil fuels electricity net generation","geothermal electricity net generation","biomass and waste electricity net generation"]

    # Compute sums for each country
    table_eia_electricity_capacity = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_eia_electricity_generation = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        columns = eia_db.columns.tolist()
        for ii in np.arange(len(years)):
            try:
                columns.index(str(years[ii]))
            except:
                continue
            
            sel_country = eia_db.iloc[:,1].values==country_name
            
            # Capacity
            sel_source = [any(keyword in elem for keyword in capacity_source_keywords) for elem in eia_db['source']]
            values = eia_db[str(years[ii])][sel_source & sel_country].astype(float)
            table_eia_electricity_capacity[jj,ii] = np.nansum(values) # Million kW units
            
            # Generation
            sel_source = [any(keyword in elem for keyword in generation_source_keywords) for elem in eia_db['source']]
            values = eia_db[str(years[ii])][sel_source & sel_country].astype(float)
            table_eia_electricity_generation[jj,ii] = np.nansum(values) # Billion kWh units
            
    # Save to csv
    pd.DataFrame(table_eia_electricity_capacity,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','eia_electricity_capacity.csv'))
    pd.DataFrame(table_eia_electricity_generation,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','eia_electricity_generation.csv'))
    '''

    ############################################################################
    #   Load AQUASTAT industrial water withdrawal data
    ############################################################################

    aquastat = pd.read_csv(os.path.join(config['aquastat_folder'],'aquastat_clean.csv'),index_col=False)

    table_aquastat_industry_withdrawal = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']    
        for ii in np.arange(len(years)):
            sel = (aquastat['m49']==country_code) & (aquastat['Variable']=='Industrial water withdrawal') & (aquastat['Year']==years[ii])
            if np.sum(sel)==0: continue
            table_aquastat_industry_withdrawal[jj,ii] = aquastat['Value'][sel].values # km3/year

    # Save to csv
    pd.DataFrame(table_aquastat_industry_withdrawal,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','aquastat_industry_withdrawal.csv'))


    ############################################################################
    #   Load World Bank gross domestic product (GDP) and manufacturing value 
    #   added (MVA) data
    ############################################################################

    # Load raw data
    gdp_file = glob.glob(os.path.join(config['world_bank_folder'],'*GDP*.csv'))[0]
    worldbank_gdp = pd.read_csv(gdp_file,skiprows=4)
    mva_file = glob.glob(os.path.join(config['world_bank_folder'],'*IND*.csv'))[0]
    worldbank_mva = pd.read_csv(mva_file,skiprows=4)

    # Values for each country and year
    table_worldbank_gdp = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_worldbank_mva = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_acronym = country_codes.iloc[jj]['alpha-3']
        sel = worldbank_mva['Country Code']==country_acronym
        if np.sum(sel)==0: continue
        for ii in np.arange(len(years)):
            if str(years[ii]) not in worldbank_gdp.columns: continue
            table_worldbank_gdp[jj,ii] = worldbank_gdp[str(years[ii])][sel] # Constant 2010 US$
            table_worldbank_mva[jj,ii] = worldbank_mva[str(years[ii])][sel] # Constant 2010 US$
                
    # Save to csv
    pd.DataFrame(table_worldbank_gdp,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','worldbank_gdp.csv'))
    pd.DataFrame(table_worldbank_mva,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','worldbank_mva.csv'))


    ############################################################################
    #   Load GCAM regional industry and thermoelectric withdrawals and rescale
    #   based on country population
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
    table_gcam_withdrawal_industry = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_gcam_withdrawal_thermoelectric = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        
        # GCAM withdrawals represent regional totals. To obtain country estimates,
        # we need to know the populations of GCAM region and country
        
        # Which region are we in?
        region_ind = -1
        for kk in np.arange(len(gcam_regions)):
            if country_name in gcam_regions[kk][1]:
                region_ind = kk
                break
        if region_ind==-1: continue
        region_name = gcam_regions[region_ind][0][0]

        # What is the total population of the region?
        region_pop = 0
        region_countries = gcam_regions[region_ind][1]
        for region_country in region_countries:
            hh = country_codes['name'].tolist().index(region_country)
            region_pop += table_population[hh,ii_2000]
        
        # Aggregate different sectors
        sel = ["Reference" in s for s in gcam_output['scenario']] &\
            (gcam_output['region']==region_name) &\
            (gcam_output['sector'].isin(['biomass', 'electricity', 'industry', 'nuclearFuelGenII', 'nuclearFuelGenIII', 'regional coal', 'regional natural gas', 'regional oil']))
        gcam_industry_withdrawal = np.sum(gcam_output.iloc[np.where(sel)[0],3:-1].values,axis=0) # km3/year
        sel = ["Reference" in s for s in gcam_output['scenario']] &\
            (gcam_output['region']==region_name) &\
            (gcam_output['sector'].isin(['biomass', 'electricity', 'nuclearFuelGenII', 'nuclearFuelGenIII', 'regional coal', 'regional natural gas', 'regional oil']))
        gcam_thermoelectric_withdrawal = np.sum(gcam_output.iloc[np.where(sel)[0],3:-1].values,axis=0) # km3/year    
        gcam_years = np.array([int(x) for x in gcam_output.columns[3:-1].values])    
        
        # Rescale regional withdrawals based on country population
        pop_frac = table_population[jj,ii_2000]/region_pop
        for kk in np.arange(len(gcam_years)):
            try:
                ii = years.tolist().index(gcam_years[kk])
            except:
                continue
            table_gcam_withdrawal_industry[jj,ii] = gcam_industry_withdrawal[kk]*pop_frac
            table_gcam_withdrawal_thermoelectric[jj,ii] = gcam_thermoelectric_withdrawal[kk]*pop_frac
        
    # Save to csv
    pd.DataFrame(table_gcam_withdrawal_industry,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_industry_withdrawal.csv'))
    pd.DataFrame(table_gcam_withdrawal_thermoelectric,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_thermoelectric_withdrawal.csv'))


    ############################################################################
    #   Load GCAM regional electricity consumption to temporally downscale 
    #   thermoelectric withdrawals
    ############################################################################

    gcam_output = pd.read_csv(os.path.join(config['gcam_folder'],'elec_consumption.csv'))

    # Values for each country and year
    table_gcam_elec_building = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_gcam_elec_trans_ind = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_gcam_elec_heating = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_gcam_elec_cooling = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_gcam_elec_other = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    for jj in np.arange(country_codes.shape[0]):
        country_name = country_codes.iloc[jj]['name']
        
        # Which region are we in?
        region_ind = -1
        for kk in np.arange(len(gcam_regions)):
            if country_name in gcam_regions[kk][1]:
                region_ind = kk
                break
        if region_ind==-1: continue
        region_name = gcam_regions[region_ind][0][0]

        # Aggregate different sectors
        sel = ["Reference" in s for s in gcam_output['scenario']] & (gcam_output['region']==region_name) & (gcam_output['input']=='elect_td_bld')
        building = np.sum(gcam_output.iloc[np.where(sel)[0],4:-1].values,axis=0)
        sel = ["Reference" in s for s in gcam_output['scenario']] & (gcam_output['region']==region_name) & (gcam_output['input'].isin(['elect_td_ind', 'elect_td_trn']))
        trans_ind = np.sum(gcam_output.iloc[np.where(sel)[0],4:-1].values,axis=0)
        sel = ["Reference" in s for s in gcam_output['scenario']] & (gcam_output['region']==region_name) & (gcam_output['sector'][:13]=='resid heating')
        heating = np.sum(gcam_output.iloc[np.where(sel)[0],4:-1].values,axis=0)
        sel = ["Reference" in s for s in gcam_output['scenario']] & (gcam_output['region']==region_name) & (gcam_output['sector'][:13]=='resid cooling')
        cooling = np.sum(gcam_output.iloc[np.where(sel)[0],4:-1].values,axis=0)
        sel = ["Reference" in s for s in gcam_output['scenario']] & (gcam_output['region']==region_name) & (gcam_output['sector'][:12]=='resid others')
        other = np.sum(gcam_output.iloc[np.where(sel)[0],4:-1].values,axis=0)
        gcam_years = np.array([int(x) for x in gcam_output.columns[4:-1].values])    
        
        # Insert into table
        for kk in np.arange(len(gcam_years)):
            try:
                ii = years.tolist().index(gcam_years[kk])
            except:
                continue
            table_gcam_elec_building[jj,ii] = building[kk]
            table_gcam_elec_trans_ind[jj,ii] = trans_ind[kk]
            table_gcam_elec_heating[jj,ii] = heating[kk]
            table_gcam_elec_cooling[jj,ii] = cooling[kk]
            table_gcam_elec_other[jj,ii] = other[kk]
                
    # Save to csv
    pd.DataFrame(table_gcam_elec_building,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_elec_building.csv'))
    pd.DataFrame(table_gcam_elec_trans_ind,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_elec_trans_ind.csv'))
    pd.DataFrame(table_gcam_elec_heating,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_elec_heating.csv'))
    pd.DataFrame(table_gcam_elec_cooling,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_elec_cooling.csv'))
    pd.DataFrame(table_gcam_elec_other,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','gcam_elec_other.csv'))


    ############################################################################
    #   Loop through countries and compute manufacturing and thermoelectric 
    #   withdrawals
    ############################################################################

    print('-------------------------------------------------------------------------------')
    print('Compute manufacturing and thermoelectric withdrawals')
    t0 = time.time()

    ii_1995 = years.tolist().index(1995)
    ii_2015 = years.tolist().index(2015)    

    table_withdrawal_industry = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_withdrawal_manufacturing = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_withdrawal_thermoelectric = np.zeros((country_codes.shape[0],len(years)),dtype=np.single)*np.NaN
    table_withdrawal_data_source = [None]*country_codes.shape[0]
    for jj in np.arange(country_codes.shape[0]):
        country_code = country_codes.iloc[jj]['country-code']
        country_name = country_codes.iloc[jj]['name']
        country_acronym = country_codes.iloc[jj]['alpha-3']    
        print('jj='+str(jj)+' '+country_name)
        
        '''
        # Lohrmann et al. (2019) underestimates capacity in many countries, so
        # correct for missing small plants using EIA data
        corr_fac = 1000*table_eia_electricity_capacity[jj,ii_2015]/table_lohrmann_plant_capacity[jj]
        if (corr_fac<0.9) | (corr_fac>1.2):
            print('Considerable discrepancy between Lohrmann et al. (2019) and EIA plant capacity (factor '+str(corr_fac)+')')
        
        # Compute thermoelectric withdrawal time series based on EIA power generation time series
        withdrawal_thermoelectric_2015 = table_lohrmann_plant_withdrawal[jj]*corr_fac
        table_withdrawal_thermoelectric[jj,:] = table_eia_electricity_generation[jj,:]*withdrawal_thermoelectric_2015/table_eia_electricity_generation[jj,ii_2015]
        
        # Compute manufacturing withdrawal time series based on MVA
        kw = dict(method="slinear", fill_value="extrapolate", limit_direction="both")
        aquastat_industry_withdrawal_2015 = pd.Series(table_aquastat_industry_withdrawal[jj,:]).interpolate(**kw).values[ii_2015]
        withdrawal_manufacturing_2015 = aquastat_industry_withdrawal_2015-withdrawal_thermoelectric_2015
        table_withdrawal_manufacturing[jj,:] = table_worldbank_mva[jj,:]*withdrawal_manufacturing_2015/table_worldbank_mva[jj,ii_2015]
        '''
        
        AQUASTAT_nvals = np.sum(np.isnan(table_aquastat_industry_withdrawal[jj,:])==False)
        GCAM_nvals = np.sum(np.isnan(table_gcam_withdrawal_industry[jj,:])==False)
        
        #--------------------------------------------------------------------------
        #   If AQUASTAT and GCAM available
        #--------------------------------------------------------------------------
        
        if (AQUASTAT_nvals>1) & (GCAM_nvals>0):
            table_withdrawal_data_source[jj] = 'AQUASTAT+GCAM'
        
            # Interpolate AQUASTAT water demand (no extrapolation)
            table_withdrawal_industry[jj,:] = pd.Series(table_aquastat_industry_withdrawal[jj,:]).interpolate(method="slinear").values
            
            # Fill AQUASTAT gaps using MVA
            ts_mva = table_worldbank_mva[jj,:]
            ts_ind = table_withdrawal_industry[jj,:].copy()
            if np.sum(np.isnan(ts_mva)==False)>1:
                ts_mva = pd.Series(ts_mva).interpolate(method="slinear").values
                indices = np.where(np.isnan(ts_ind)==False)[0]
                ts_ind[:indices[0]] = ts_mva[:indices[0]]*ts_ind[indices[0]]/ts_mva[indices[0]]
                ts_ind[indices[-1]+1:] = ts_mva[indices[-1]+1:]*ts_ind[indices[-1]]/ts_mva[indices[-1]]
                table_withdrawal_industry[jj,:] = ts_ind.copy()
            
            # Fill remaining AQUASTAT gaps using population
            ts_ind = table_withdrawal_industry[jj,:]
            ts_pop = table_population[jj,:]
            if np.sum(np.isnan(ts_ind))>0:
                indices = np.where(np.isnan(ts_ind)==False)[0]
                ts_ind[:indices[0]] = ts_pop[:indices[0]]*ts_ind[indices[0]]/(ts_pop[indices[0]]+1) # +1 to avoid divide by zero
                ts_ind[indices[-1]+1:] = ts_pop[indices[-1]+1:]*ts_ind[indices[-1]]/(ts_pop[indices[-1]]+1) # +1 to avoid divide by zero
                table_withdrawal_industry[jj,:] = ts_ind.copy()
           
            # Linearly interpolate and nearest-neighbor extrapolate GCAM data
            ts_gcam_withdrawal_industry = pd.Series(table_gcam_withdrawal_industry[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values
            ts_gcam_withdrawal_thermoelectric = pd.Series(table_gcam_withdrawal_thermoelectric[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values
            ts_gcam_withdrawal_manufacturing = ts_gcam_withdrawal_industry-ts_gcam_withdrawal_thermoelectric
            
            # Rescale GCAM to Vassolo and Doll (2005), if not missing and if country large 
            # enough (i.e., covers enough grid cells)
            country_ncells = np.sum(country_code_map_360x720==country_code)
            if (np.isnan(vassolo_doll_withdrawal_thermoelectric[jj]+vassolo_doll_withdrawal_manufacturing[jj])==False) & (country_ncells>20):
                ts_gcam_withdrawal_thermoelectric = ts_gcam_withdrawal_thermoelectric*vassolo_doll_withdrawal_thermoelectric[jj]/(ts_gcam_withdrawal_thermoelectric[ii_1995]+10**-6) # +10**-6 to avoid divide by zero
                ts_gcam_withdrawal_manufacturing = ts_gcam_withdrawal_manufacturing*vassolo_doll_withdrawal_manufacturing[jj]/(ts_gcam_withdrawal_manufacturing[ii_1995]+10**-6) # +10**-6 to avoid divide by zero
                ts_gcam_withdrawal_industry = ts_gcam_withdrawal_thermoelectric+ts_gcam_withdrawal_manufacturing    
            
            # Temporally disaggregate industrial into manufacturing and thermoelectric using GCAM
            table_withdrawal_thermoelectric[jj,:] = table_withdrawal_industry[jj,:]*ts_gcam_withdrawal_thermoelectric/(ts_gcam_withdrawal_industry+10**-6) # +10**-6 to avoid divide by zero
            table_withdrawal_manufacturing[jj,:] = (table_withdrawal_industry[jj,:]+10**-6)-table_withdrawal_thermoelectric[jj,:]


        #--------------------------------------------------------------------------
        #   If AQUASTAT not available, use GCAM
        #--------------------------------------------------------------------------
        
        elif (AQUASTAT_nvals<=1) & (GCAM_nvals>0):
            table_withdrawal_data_source[jj] = 'GCAM'
            
            # Linearly interpolate and nearest-neighbor extrapolate GCAM data
            table_withdrawal_industry[jj,:] = pd.Series(table_gcam_withdrawal_industry[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values
            table_withdrawal_thermoelectric[jj,:] = pd.Series(table_gcam_withdrawal_thermoelectric[jj,:]).interpolate(method="linear",fill_value="extrapolate", limit_direction="both").values
            table_withdrawal_manufacturing[jj,:] = table_withdrawal_industry[jj,:]-table_withdrawal_thermoelectric[jj,:]
            
        elif GCAM_nvals==0:
            table_withdrawal_data_source[jj] = 'None'
            table_withdrawal_industry[jj,:], table_withdrawal_thermoelectric[jj,:], table_withdrawal_manufacturing[jj,:] = 0, 0, 0
        
        
        #--------------------------------------------------------------------------
        #   Plot figure
        #--------------------------------------------------------------------------

        # # Initialize figure
        # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        
        # # Subpanel 1
        # ax1.plot(years,table_aquastat_industry_withdrawal[jj,:],'*',label='FAO industry withdrawal',color=(0.5,0.5,0.5))
        # ax1.plot(years,table_gcam_withdrawal_industry[jj,:],'+',label='GCAM industry withdrawal',color=(0,0.5,0))
        # ax1.plot(years,table_gcam_withdrawal_industry[jj,:]-table_gcam_withdrawal_thermoelectric[jj,:],'+',label='GCAM manufacturing withdrawal',color=(0,1,0))
        # ax1.plot(years,table_withdrawal_industry[jj,:],label='Beck industry withdrawal',color=(0.5,0,0))
        # ax1.plot(years,table_withdrawal_manufacturing[jj,:],label='Beck manufacturing withdrawal',color=(1,0,0))
        # ax1.legend()
        # ax1.set_ylabel('Withdrawal (km^3/year)')
        # ax1.set_title(country_name)
        
        # # Subpanel 2
        # ax2.plot(years,table_Huang_withdrawal_industry[jj,:],label='Huang industry withdrawal',color=(0,0,0.5))
        # ax2.plot(years,table_Huang_withdrawal_manufacturing[jj,:],label='Huang manufacturing withdrawal',color=(0,0,1))
        # ax2.plot(years,table_withdrawal_industry[jj,:],label='Beck industry withdrawal',color=(0.5,0,0))
        # ax2.plot(years,table_withdrawal_manufacturing[jj,:],label='Beck manufacturing withdrawal',color=(1,0,0))
        # ax2.legend()
        # ax2.set_ylabel('Withdrawal (km^3/year)')
        # ax2.set_title(country_name)
            
        # # Subpanel 3
        # ax3.plot(years,table_worldbank_gdp[jj,:],label='GDP',color='g')
        # ax3.plot(years,table_worldbank_mva[jj,:],label='MVA',color='b')
        # ax3.legend()
        # ax3.set_ylabel('Constant 2010 US$')
        # ax3.set_title(country_name)    
        
        # # Save figure
        # f.set_size_inches(10, 10)
        # plt.savefig(os.path.join(config['output_folder'],'step3a_industrial_demand','figures',str(jj).zfill(3)+'_'+country_name+'.png'),dpi=150)
        # plt.close()
            
    # Save to csv
    pd.DataFrame(table_withdrawal_data_source,index=country_codes['name'],columns=['Withdrawal data source']).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','withdrawal_data_source.csv'))
    pd.DataFrame(table_withdrawal_industry,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','withdrawal_industry.csv'))
    pd.DataFrame(table_withdrawal_manufacturing,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','withdrawal_manufacturing.csv'))
    pd.DataFrame(table_withdrawal_thermoelectric,index=country_codes['name'],columns=years).to_csv(os.path.join(config['output_folder'],'step3a_industrial_demand','tables','withdrawal_thermoelectric.csv'))

    print("Time elapsed is "+str(time.time()-t0)+" sec")


    ############################################################################
    #   Compute monthly maps of heating degree days (HDD) and cooling degree 
    #   days (CDD)
    ############################################################################

    print('-------------------------------------------------------------------------------')

    for year in years:
        for month in np.arange(1,13): 
            t0 = time.time()
            
            outpath_hdd = os.path.join(config['output_folder'],'step3a_industrial_demand','hdd',str(year)+str(month).zfill(2)+'.npz')
            outpath_cdd = os.path.join(config['output_folder'],'step3a_industrial_demand','cdd',str(year)+str(month).zfill(2)+'.npz')
            
            if os.path.isfile(outpath_hdd) & os.path.isfile(outpath_cdd): continue
            
            print('Computing HDD and CDD for '+str(year)+' month '+str(month))
            
            # Load MSWX air temperature data for entire month
            ndays = monthrange(year,month)[1]
            mswx_month = np.zeros((1800,3600,ndays),dtype=np.single)*np.NaN
            for day in np.arange(1,ndays+1):
                if year<1979:   #MSWX dataset starts from 1979
                    try:
                        date=datetime(1979,month,day)
                    except: # fix for Leap years < 1979
                        date=datetime(1979,month,day-1)
                    filepath = os.path.join(config['mswx_folder'],'Past','Temp','Daily',date.strftime('%Y%j')+'.nc')
                else:
                    filepath = os.path.join(config['mswx_folder'],'Past','Temp','Daily',datetime(year,month,day).strftime('%Y%j')+'.nc')
                if os.path.isfile(filepath):
                    dset = Dataset(filepath)
                    mswx_month[:,:,day-1] = np.squeeze(np.array(dset.variables['air_temperature']))
            
            # Compute HDD and CDD
            mswx_hdd = np.nansum((18-mswx_month)*(mswx_month<18),axis=2)
            mswx_cdd = np.nansum((mswx_month-18)*(mswx_month>=18),axis=2)
            
            # Save results
            np.savez_compressed(outpath_hdd,data=mswx_hdd)
            np.savez_compressed(outpath_cdd,data=mswx_cdd)
            
            print("Time elapsed is "+str(time.time()-t0)+" sec")


    # pdb.set_trace()
    

if __name__ == '__main__':
    main()