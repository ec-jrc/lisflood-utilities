# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:50:54 2021

@author: tilloal
"""

import cdsapi
import requests
import os
from tools import *

## Script 1/4

# CDS API script to use CDS service to retrieve hourly ERA5* variables and iterate over
# all months in the specified years.
 
# Requires:
# 1) the CDS API to be installed and working on your system
# 2) You have agreed to the ERA5 Licence (via the CDS web page)
# 3) Selection of required variable, daily statistic, etc
 
# Output:
# 1) separate netCDF file for chosen variables at hourly resolution for each month

dir_path = os.path.dirname(os.path.realpath(__file__))
# Change working directory
os.chdir(dir_path)

config = load_config("cds_config.cfg")
# key text file containing (1) proxy key, (2) CDS API key, (3) working directory of the data

# Proxy settings
proxy = config ['proxyKey']
os.environ['http_proxy'] = proxy 
os.environ['HTTP_PROXY'] = proxy
os.environ['https_proxy'] = proxy
os.environ['HTTPS_PROXY'] = proxy

CDSAPI_URL="https://cds.climate.copernicus.eu/api/v2"
# Load personnal API Key
CDSAPI_KEY=config['CDSAPI_KEY']
c = cdsapi.Client(key=CDSAPI_KEY, url=CDSAPI_URL) 

# select variable(s); name must be a valid ERA5 CDS API name.
# Here I select varaibles necessary to compute evapotranspiration
input_string = input('Enter varialbes that need to be processed from the list [u10,v10,2t,str,ssrd,d2m], enter "all" for all variables')
print("\n")
vin = input_string.split()
if vin == ["all"]:
    varnames = ['u10', 'v10','2t','str','ssrd','2d','tp']
else:
    varnames=vin
    
# print list
print('list: ', varnames)


#this list has to be modified to include more variables
pair_var=[['u10', 'v10','2t','str','ssrd','2d','tp'],
          ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature','surface_net_thermal_radiation', 'surface_downward_solar_radiation','2m_dewpoint_temperature','total_precipitation']]

pair_match=set(pair_var[0]).intersection(varnames)
nf = [item for item,x in enumerate(pair_var[0]) if x in varnames]

var=[pair_var[1][i] for i in nf]

print('variable(s) '+' | '.join(var) + ' will be downloaded')

                                        #,'surface_net_thermal_radiation', 'surface_net_solar_radiation','2m_dewpoint_temperature']
#var = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature','surface_net_thermal_radiation', 'surface_net_solar_radiation','2m_dewpoint_temperature']
#varnames = ['u10', 'v10','2t','str','ssrd','2d']
date = config['date']


# define area where data needs to be extracted (long/lat)
area=[config['Lat1'], config['Lon1'], config['Lat2'], config['Lon2'],]

#automatic selection of date range according to user input
if date == 'auto':
    start=  int(sys.argv[1])
    end = int(sys.argv[2])
    out_Yrs = list(range(start,end+1))
    years = out_Yrs

#manual year selection
if date == "manual":
    # Uncomment years as required
    years =  [
    #'1981','1982', '1983','1984',
    #'1985', '1986', 'u101987',
    '1988',
    #'1989', '1990',
    #'1991', '1992', '1993','1994', '1995', '1996','1997', '1998',
    #'1999',
    #'2000', '2001', '2002','2003', '2004', '2005','2006', '2007', '2008',
    #'2009', '2010', '2011','2012', '2013', '2014','2015', '2016', '2017',
    #'2018', '2019', '2020',
    
    ]
 
months = [ "01", "02", "03", "04", "05", "06", 
         "07", "08", "09", "10", "11", "12" ] 

# Choose directory where to download files
dire=config["download_folder"]
direh=dire + "\hourly"
os.chdir(direh)

# main loop to download variables
# monthly files
# all variables will be contained in one netcdf
for yr in years:
    for mn in months:

        c.retrieve(
        'reanalysis-era5-land',
        {
            'variable': var,
            'year': yr,
            'month': mn,
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area':area ,
            'format': 'netcdf',
        },
        "e5l_lf-met" + "_" + yr + "_" + mn + ".nc")
