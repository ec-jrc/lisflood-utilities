#!/bin/bash
#$ -S /bin/sh


# L.A., 1/3/2019
# Take discharge files (in netCDF format), extract annual maxima and then fit an extreme value distribution. Produce output maps in PCRaster and netcdf format

# usage: ./Gfit2.sh ./settingFile_test.sh
# usage (e.g.): qsub -V -j y -N GFit -q all.q -o /nahaUsers/alfielo/CA/GloFAS/scripts/GumbelFit_v2/GFIT_tool.v2/GumFit_nc.log /nahaUsers/alfielo/CA/GloFAS/scripts/GumbelFit_v2/GFIT_tool.v2/Gfit2.sh /nahaUsers/alfielo/CA/GloFAS/scripts/GumbelFit_v2/GFIT_tool.v2/settingFile_test.sh

settFile=${1}
echo $settFile
source $settFile

set -o allexport
source $settFile
eval echo "source $settFile"             #Import settings file
set +o allexport

#######################################################################
# check default values for input variables
if [ -z "$warmup_yrs" ]; then
    echo "warmup_yrs is unset"
    warmup_yrs=0
fi

if [ -z "$MaxGTavgDis" ]; then
    echo "MaxGTavgDis is unset"
    MaxGTavgDis=1
fi

if [ -z "$Nyears" ]; then
    echo "Nyears is unset"
    Nyears="NaN"
fi
if [ -z "$varNumber" ]; then
    echo "varNumber is unset"
    varNumber=1
fi
#######################################################################


mkdir -p $outDir

########################################################################
# cdo operations
# Delete first ( number of ) warmup_yrs from the cdo file

echo "Check warm-up years"
if [ -v warmup_yrs ] ; then
  if [ $warmup_yrs -ge 1 ] ; then
    years=`cdo showyear $inpDir/$varName.nc`
    arr=($years)
    echo cdo delete,year=$(seq -s, ${arr[0]} ${arr[warmup_yrs-1]}) $inpDir/$varName.nc $outDir/${varName}2.nc
    cdo delete,year=$(seq -s, ${arr[0]} ${arr[warmup_yrs-1]}) $inpDir/$varName.nc $outDir/${varName}2.nc
  else
    ln -sf $inpDir/$varName.nc $outDir/${varName}2.nc
  fi
else
  ln -sf $inpDir/$varName.nc $outDir/${varName}2.nc
fi

# Extract the series of annual maxima from a netCDF file of the time series
echo "Extract the series of annual maxima"
cdo yearmax $outDir/${varName}2.nc $outDir/${varName}AMax.nc

# Extract the average map over the entire time series
echo "Extract the average map"
cdo timmean $outDir/${varName}2.nc $outDir/${varName}Avg.nc

# Copy pcraster clone map
echo "Copy pcraster clone map"
cp $cloneMap $outDir
########################################################################
# Fitting of Gumbel extreme value distributions with L-moments
echo "Fit of Gumbel extreme value distributions with L-moments"

/usr/bin/Rscript ${rootdir}gfit2.r $outDir/${varName}AMax.nc $outDir $returnPeriods 0 $MaxGTavgDis $Nyears 

########################################################################
# Transform ArcView readable ASCII files into PCRaster maps
if [ -v cloneMap ]; then 
  echo "create PCRaster maps"
  cd $outDir
  Files=$(ls *.txt | sed -e 's/\..*$//')
  arr=($Files)
  let nFiles=${#arr[@]}-1 #number of models minus 1 (because the first is #0)
  
  for iFile in $(seq 0 $nFiles)  
  do
    M=${arr[iFile]}
    echo ${arr[iFile]}
    asc2map --clone $cloneMap -S -a -m 1e31 ${arr[iFile]}.txt ${arr[iFile]}.map
  done
fi

########################################################################
# Some more cdo operations to extract statistics

# Extract min and maximum
echo "extract min, max, and percentile maps"
cdo timmin $outDir/${varName}2.nc $outDir/${varName}Min.nc
cdo timmax $outDir/${varName}2.nc $outDir/${varName}Max.nc

# Extract percentile maps
cdo timpctl,1 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}01p.nc   #1 %
cdo timpctl,2 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}02p.nc   #2 %
cdo timpctl,5 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}05p.nc   #5 %
cdo timpctl,10 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}10p.nc  #10%
cdo timpctl,25 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}25p.nc  #25%
cdo timpctl,50 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}50p.nc  #50%
cdo timpctl,75 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}75p.nc  #75%
cdo timpctl,90 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}90p.nc  #90%
cdo timpctl,95 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}95p.nc  #95%
cdo timpctl,98 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}98p.nc  #98%
cdo timpctl,99 $outDir/${varName}2.nc $outDir/${varName}Min.nc $outDir/${varName}Max.nc $outDir/${varName}99p.nc  #99%

rm $outDir/${varName}2.nc
########################################################################
