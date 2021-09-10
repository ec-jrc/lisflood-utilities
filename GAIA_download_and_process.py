import os, sys, glob, time, pdb
import pandas as pd
import geopandas as gpd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import subprocess
import rasterio
from tools import *
    
# Load configuration file
config = pd.read_csv(sys.argv[1],header=None,index_col=False)
for ii in np.arange(len(config)): 
    string = config.iloc[ii,0]
    string = string.replace(" ","")    
    try:
        exec(string.replace("=","=r")) 
    except:
        exec(string) 

if os.path.isdir(os.path.join(gaia_folder,'geotiffs'))==False:
    os.makedirs(os.path.join(gaia_folder,'geotiffs'))
    
if os.path.isdir(os.path.join(gaia_folder,'netCDF_1km'))==False:
    os.makedirs(os.path.join(gaia_folder,'netCDF_1km'))

if use_proxy==True:
    proxy_credentials = pd.read_csv('proxy_credentials.cfg',header=None,index_col=False)
    
# Load GAIA shapefile
GAIA_shape = gpd.read_file(os.path.join(gaia_folder,'GAIA_shape','GAIA_nameID_1deg.shp'))


############################################################################
#   Download and untar all MERIT Hydro upstream data
############################################################################

url_pre = 'http://data.ess.tsinghua.edu.cn/data/GAIA/'

for ii in np.arange(GAIA_shape.shape[0]):    
    string = GAIA_shape.loc[ii,'fName_ID']
    
    # Fix string padding...
    underscore = string.find('_')
    lonstr = str(string[underscore+1:]).zfill(4)
    if lonstr[0]=='0': lonstr = lonstr[1:]
    latstr = str(string[:underscore]).zfill(3)
    if latstr[0]=='0': latstr = latstr[1:]
    string = latstr+'_'+lonstr
    
    filename = 'GAIA_1985_2018_'+string+'.tif'
    if os.path.isfile(os.path.join(gaia_folder,'geotiffs',filename)): continue
    if use_proxy==True:
        proxy_string = '-e use_proxy=yes -e http_proxy=http://'+proxy_credentials.iloc[0,0]+':'+proxy_credentials.iloc[1,0]+'@'+proxy_credentials.iloc[2,0]+' '
    else:
        proxy_string = ''
    print('===============================================================================')
    command = 'wget '+proxy_string+url_pre+filename+' --no-clobber --directory-prefix='+os.path.join(gaia_folder,'geotiffs')
    subprocess.call(command,shell=True)


############################################################################
#   Loop over years and make impervious maps
############################################################################

years = np.arange(1985,2019)
values = np.arange(34,0,-1)
if len(years)!=len(values):
    raise ValueError('Something wrong')
    
geotiff_files = glob.glob(os.path.join(gaia_folder,'geotiffs','*.tif'))
for jj in np.arange(len(years)):
    t0 = time.time()
    year = years[jj]
    value = values[jj]
    print('Processing '+str(year))
    
    output_res = 0.01
    global_impervious = np.zeros((int(180/output_res),int(360/output_res)),dtype=np.single)
    for ii in np.arange(len(geotiff_files)):
        
        # Load and resample tile
        raw = rasterio.open(geotiff_files[ii]).read(1)
        impervious = raw>=value
        impervious = imresize_mean(np.single(impervious),(np.int(1/output_res),np.int(1/output_res)))
           
        # Find top row and left column of tile
        filename = os.path.basename(geotiff_files[ii])
        try:
            tile_lat_top = float(filename[15:18])
        except:
            tile_lat_top = float(filename[15:17])
        try:
            tile_lon_left = float(filename[-8:-4])
        except:
            tile_lon_left = float(filename[-7:-4])
        tile_row_top, tile_col_left = latlon2rowcol(tile_lat_top-output_res/2,tile_lon_left+output_res/2,output_res,90,-180)
        
        # Insert tile into global map
        global_impervious[tile_row_top:tile_row_top+impervious.shape[0],tile_col_left:tile_col_left+impervious.shape[1]] = impervious
        
    # Save global map to netCDF
    save_netcdf(
        file=os.path.join(gaia_folder,'netCDF_1km',str(year)+'.nc'), 
        varname='impervious_fraction', 
        data=global_impervious, 
        least_sig_dig=3, 
        lat=np.arange(90-output_res/2,-90-output_res/2,-output_res), 
        lon=np.arange(-180+output_res/2,180+output_res/2,output_res)
        )
        
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
pdb.set_trace()
        
        

        
        
        
        # Insert into global map

        #save()
'''
    global_shape = (int(180/res),int(360/res))
    upstream_area_global = np.zeros(global_shape)*np.NaN
    
    for subdir, dirs, files in os.walk(gaia_folder):
        for file in files:        
            
            
            print('--------------------------------------------------------------------------------')
            print('Processing '+file)
            t1 = time.time()
            
            # Resize using maximum filter
            oldarray = rasterio.open(os.path.join(subdir, file)).read(1)
            oldarray[oldarray<0] = 9999999 # Necessary to ensure that all rivers flow into the ocean
            factor = res/(5/6000)
            factor = np.round(factor*1000000000000)/1000000000000
            if factor!=np.round(factor): 
                raise ValueError('Resize factor of '+str(factor)+' not integer, needs to be integer')
            factor = factor.astype(int)
            newshape = (oldarray.shape[0]//factor,oldarray.shape[1]//factor)
            newarray = imresize_max(oldarray,newshape)
            
            # Insert into global map
            tile_lat_top = float(file[:3].replace("n","").replace("s","-"))
            tile_lon_left = float(file[3:7].replace("e","").replace("w","-"))
            tile_row_top, tile_col_left = latlon2rowcol(tile_lat_top+res/2,tile_lon_left+res/2,res,90,-180)
            upstream_area_global[tile_row_top-newshape[0]+1:tile_row_top+1,tile_col_left:tile_col_left+newshape[1]] = newarray
            
            print('Time elapsed is ' + str(time.time() - t1) + ' sec')

    with open(os.path.join(output_folder,'upstream_area_global.npy'), 'wb') as f:
        np.save(f, upstream_area_global)

# Load resampled global upstream area map from disk
with open(os.path.join(output_folder,'upstream_area_global.npy'), 'rb') as f:
    upstream_area_global = np.load(f)


'''