# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:36:51 2022

@author: tilloal
"""



# Import the os module
import os
import rioxarray # for the extension to load
import pandas as pd
import xarray as xr
import numpy as np
import math
from scipy.interpolate import griddata
from xarray import concat
from netCDF4 import Dataset
import netCDF4 as nc
import time

vr="ta"
yr="1990"
d1=os.getcwd() + '/' + vr 
#%%
d2=  '/1_arcmin/e5ld_1min_lvap_' + vr + '_'+ yr + '.nc'
d3=  '/0.1_deg/e5ldxp_' + vr + '_'+ yr + '.nc'
f=d1+d2

tamere=Dataset(f, 'r',scale=False)

#tamere.set_auto_mask(False)

jeve=tamere.variables[vr]
grd_mere = xr.open_dataset(f)

gme=grd_mere[vr]

maps=gme[3,:,:] 

#%%

maps.plot() 
#%%
mamene=jeve[range(200,300),2400,2000]

mamere=gme[range(200,300),2400,2000] 
#mamere.plot()

ptn=mamere.values
#%%
print(mamene)
print(ptn)
#%%
jeve[:].flatten()[69] 
#%%
print(jeve)
print(gme)