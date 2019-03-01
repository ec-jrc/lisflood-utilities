#!/bin/bash
# Settings file for the Gumbel fit

inpDir="/climateRun4/HELIX/global/dis/LF05/out"   			#Input file with daily maps to analyze and from which to extract the annual maxima
varName="dis"                                     			#Output filename 
outDir="/H07_Global/deleteMe"    							#Output directory
cloneMap="/climateRun4/HELIX/global/maps/area.map"  		#Clone map (Only in case you want the output maps also in PCRaster format, otherwise omit)
returnPeriods="2-5-10-20-50-100-200-500"              		#Return periods [years]. Use the dash (-) as separator
rootdir="/nahaUsers/alfielo/CA/GloFAS/scripts/GumbelFit_v2/GFIT_tool.v2/" 	#Directory where gfit2.r is stored, with "/" at the end (omit if current directory) 
warmup_yrs="0"   											#Number of years to discard in the beginning 
# Nyears="30"      										#Number of years to use in the fitting (omit if all)
MaxGTavgDis=0												#if 1 then annual maxima which are below the long term average are ignored in the fitting (it may result in missing values in the return levels)
