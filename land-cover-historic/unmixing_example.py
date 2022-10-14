#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "September 2021"

import os, sys, glob, time, pdb
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tools import *

# Dimensions of raw and output data
shape_hi = (18000,36000)
shape_lo = (1800,3600)

# Prepare dummy input data
param_hi = np.random.random(shape_hi).astype(np.single)
param_hi[0,:] = param_hi[0,:]+3
param_hi[1,:] = param_hi[1,:]+2
param_hi[3,:] = param_hi[3,:]-0.2
param_hi[4,:] = param_hi[4,:]-0.5
tree_hi = np.random.random(shape_hi).astype(np.single)
tree_hi[0,:] = tree_hi[0,:]+1.5
tree_hi[1,:] = tree_hi[1,:]+1
tree_hi[3,:] = tree_hi[3,:]-0.1
tree_hi[4,:] = tree_hi[4,:]-1

t0 = time.time()

# Resample maps
param_lo = imresize_mean(param_hi,shape_lo)
param_avg = imresize_mean(param_lo,shape_hi)
tree_lo = imresize_mean(tree_hi,shape_lo)
tree_avg = imresize_mean(tree_lo,shape_hi)
factor = shape_hi[0]/shape_lo[0]
no_of_pixels = np.zeros(shape_lo,dtype=np.single)+factor**2

# Calculate linear regression slope and intercept
# https://www4.stat.ncsu.edu/~dickey/summer_institute/formulas
# Y=param; X=tc
slope = imresize_mean((param_hi-param_avg)*(tree_hi-tree_avg),shape_lo)*no_of_pixels/(imresize_mean((tree_hi-tree_avg)**2,shape_lo)*no_of_pixels)
intercept = param_lo-slope*tree_lo

# Estimated param values at 0 % and 100 % forest cover
param_tree0 = slope*0+intercept
param_tree1 = slope*1+intercept

print("Time elapsed is "+str(time.time()-t0)+" sec")

# Verification for single grid-cell [0,0]
print(param_tree0[0,0])
print(param_tree1[0,0])
y = param_hi[:int(factor),:int(factor)].flatten()
x = tree_hi[:int(factor),:int(factor)].flatten()
plt.scatter(x,y)
x = np.arange(-1,2,0.01)
y = x*slope[0,0]+intercept[0,0]
plt.plot(x,y)
plt.xlabel('Forest cover')
plt.ylabel('Parameter value')
plt.show()

pdb.set_trace()