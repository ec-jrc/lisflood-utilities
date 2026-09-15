#! /bin/ksh

#SBATCH --time=48:00:00
#SBATCH --job-name=YYYY_input_GloFAS_generation
#SBATCH --output=yearly_logs/YYYY_input_GloFAS_generation.out

# set -ex  # this prints all commands and is useful for debugging (it makes log tooooo long though)

module load eclib/1.1.0
module load gnuparallel/20210222
module load nco/4.9.7
module load cdo
module load gdal/3.8.4 # module load gdal/3.2.1
module load ecmwf-toolbox/2024.04.0.0 # module load ecmwf-toolbox/2024.02.1.0 # ecmwf-toolbox/2021.12.0.0
module load pyg2p/3.2.7-dev # module load pyg2p/dev_arcmin
module load python3/3.11.10-01 # module load python3/3.8.8-01
export PYTHONPATH=/usr/local/apps/pyg2p/3.2.7-dev/lib/python3.11/site-packages
module load pcraster/4.4.1-01 # module load pcraster/4.3.0
module load climetlab # climetlab/0.13.11  # save netcdf to grb (for ERA5 rainbomb correction script)

# -------------------------
# Define dirs and dates
# -------------------------
echo "Data generation for year YYYY_input started."
YYYY=YYYY_input

# Consistency with GloFAS 4.0 in the forcing generation:
# GloFAS 4.0 forcings are based on files with 1st Jan 31st Dec each year. The dates refer to end date so
# the actual dates retrieved per year should be 31st Dec of previous year until 30 Dec of reference year
YYYYm1="$(($YYYY-1))"
if [[ $YYYY -eq 1975 ]]; then
    fdate=${YYYY}0101
else
    fdate=${YYYYm1}1231
fi
ldate=${YYYY}1230
# use all days of the year in the same nc file (the above was keeping last day in next year's file as with GloFASv4)
fdate=${YYYY}0101
ldate=${YYYY}1231

# below 2 lines are just for checking; should be commented out in final runs
# fdate=${YYYY}0405
# ldate=${YYYY}0406

pardir_name=glofas_forcing_V5.3

# use scratch for analysis, but then copy all final data in perm, in a directory of the same name
OUTDIR=/scratch/ecm8227/$pardir_name #/ec/vol/cems_floods/moi/ERA5/GLOFAS_forcing
OUTDIRPERM=/perm/ecm8227/$pardir_name
PYG2PTEMPDIR=/home/ecm8227/Repositories/glofas_forcings/version5/glofas_execution_templates  # templates for pyg2p interpolations
PYMASKDIR=/home/ecm8227/Repositories/glofas_forcings/version5/  # python script for masking data and keeping the encoding
mkdir -p $TMPDIR/$YYYY
cd $TMPDIR/$YYYY
OUTDIR=$OUTDIR/$YYYY 

#tpO: original rainfall, tp: after rainbomb correction; the latter is the suggested GloFAS forcing
for param in tpO tpO_grb tp_grb tp ta e0 eT0 eS0; do
    mkdir -p $OUTDIR/$param   
    mkdir -p $OUTDIR/../$param   
done

# save original files, so we are not bound by MARS in case we want to run only a part of the workflow
for param in tp rg rn wu wv ta td; do
    mkdir -p $OUTDIRPERM/aux/$YYYY/${param}_MARS
done

# Define cutoff date for ERA5T
currymd=`date +%Y%m%d`
currym=$(substring $currymd 1 6)
ym2=$(newdate -D ${currym}01 -57)
era5tymd=$(substring $ym2 1 6)01

# -------------------------
# Run analysis per day
# -------------------------
idate=$fdate
while [[ $idate -le $ldate ]]; do
    # -------------------------
    # Vars, dirs for idate
    # -------------------------
    rm -rf $TMPDIR/${YYYY}/*
    idatem1=$(newdate -D $idate -1)  # previous day
    idate1=$(newdate -D $idate +1)  # next day
    yyyy=$(substring $idate 1 4)
    mm=$(substring $idate 5 6)
    dd=$(substring $idate 7 8)
    mkdir -p $OUTDIR/$mm
    
    if [[ $idate -ge $era5tymd ]]; then
        ERA5EXPVER=5
    else
        ERA5EXPVER=1
    fi
    
    # -------------------------
    # MARS retrieval
    # -------------------------
    # at start of each year stage all data for speeding up the MARS retrieval later on
    # if [[ $mm = '12' && $dd = '31' ]]; then
    if [[ $idate = $fdate ]]; then
        cat > ret.txt <<EOF
        stage,class=ea,date=$idatem1/to/$(newdate -D $idatem1 +368),expver=$ERA5EXPVER,type=fc,levtype=sfc,param=228/169/177,stream=oper,time=6/18,step=7/to/18/by/1
        stage,type=an,param=165/166/167/168,time=0/to/23/by/1,step=0
EOF
        mars ret.txt
    fi

    # if the last file for the date is available, then skip date (MARS takes time to download data...)
    if [ ! -f ${OUTDIR}/tp/tp_${idate1}.nc ]; then
    
       # for each day get the hourly fields (from forecasted data for fluxes and analysis data for instanteneous variables). Initiate an empty text file, append all retrievals (per day) and then call it with mars
       cat > ret.txt <<EOF
        
EOF
    
       # for fluxes we get steps 7 to 18 from previous date 18:00 init, and same day 06:00 init time
       # MARS keeps previous arguments if not available, so the next retrievals of a variable has many arguments ommitted
        var_names=(TP SSRD STR)
        var_codes=(228 169 177)
        for i in "${!var_names[@]}"; do
            i_name=${var_names[$i]}
            i_code=${var_codes[$i]}
            
            cat >> ret.txt <<EOF
            retrieve,class=ea,date=$idatem1,expver=$ERA5EXPVER,type=fc,levtype=sfc,param=$i_code,stream=oper,time=18,step=7/to/18/by/1,field=r1
            retrieve,date=$idate,time=06,step=7/to/18/by/1,field=r2
            compute,field=f0,formula="merge(r1,r2)"
            compute,field=f,formula="sum(f0)"
            write,field=f,target="${OUTDIR}/$mm/era5_${i_name}_${idate}.grb"
EOF
            done
       
       # for instantaneous get 0-23 of same day and 0 of next day and average them (half weight on the two 0 UTCs)
        var_names=(10U 10V T2 T2d)
        var_codes=(165 166 167 168)
        for i in "${!var_names[@]}"; do 
            i_name=${var_names[$i]}
            i_code=${var_codes[$i]}
            
            cat >> ret.txt <<EOF
        retrieve,class=ea,date=$idate,expver=$ERA5EXPVER,type=an,levtype=sfc,param=$i_code,stream=oper,step=0,time=0,field=f1
        retrieve,date=$idate,time=1/to/23/by/1,field=f2
        retrieve,date=$idate1,time=0,field=f3
        compute,field=g,formula="(.5*f1+sum(f2)+.5*f3)/24"
        write,field=g,target="${OUTDIR}/$mm/era5_${i_name}_${idate}.grb"
EOF
        done
        
        mars ret.txt
    
       # -------------------------
       # Modify retrieved files
       # -------------------------
       # change grib file fields so the data refer to end of current day (compatible with pyg2p accumulation json)
        mars_names=(TP SSRD STR 10U 10V T2 T2d)
        modified_names=(tp rg rn wu wv ta td)
        for i in "${!mars_names[@]}"; do 
            i_mars=${mars_names[$i]}
            i_mod=${modified_names[$i]}
            grib_set -s step=24,dataDate=$idate,dataTime=0 ${OUTDIR}/$mm/era5_${i_mars}_${idate}.grb ${i_mod}_0.grb
            # copy original MARS data (after correcting the grib fields) so we don't need MARS retrieval if script needs to run again
            cp ${i_mod}_0.grb $OUTDIRPERM/aux/$YYYY/${i_mod}_MARS/${i_mod}_${idate}.grb
        done
    
        # copy the original TP (after correcting the grib fields) for the rainbomb correction to be done later on
        cp tp_0.grb ${OUTDIR}/tpO_grb/tpO_grb_${idate}.grb

        
        # -------------------------
        # Run pyg2p interpolation
        # -------------------------
        # PYG2P grib to PCRASTER (netcdf)
        for param in tp rg rn wu wv ta td; do
            pyg2p -i ${param}_0.grb -o . -c $PYG2PTEMPDIR/AFFS_${param}.json
            cp ${param}_None.nc ${param}.nc  # copy for LISVAP input
        done
    
        # Save forcings (date refers to end date as LISFLOOD conventions and data are masked to include only the actual domain)
        python3 ${PYMASKDIR}/mask_compress.py -i tp_None.nc -o ${OUTDIR}/tpO/tpO_${idate1}.nc
        python3 ${PYMASKDIR}/mask_compress.py -i ta_None.nc -o ${OUTDIR}/ta/ta_${idate1}.nc
    
        # -------------------------
        # Run LISVAP
        # -------------------------
        suite_libdir=/home/ecm8227/Repositories/glofas_forcings/version5  # /perm/moi/lib
        lisvap_new_env=/home/ecm8227/Repositories/lisflood-lisvap
        PathMapsLF=/ec/ws4/tc/emos/work/cems/floods/glofas/assets/4.0/maps #/ec/vol/cems_floods/glofas/static_maps/4.0
        PathOutLVP=./ #$sim_fcdir
        nd=$(newdate -D $idate +0)
        StartDate="$(substring $nd 7 8)/$(substring $nd 5 6)/$(substring $nd 1 4) 00:00"
        StepStarT="$(substring $nd 7 8)/$(substring $nd 5 6)/$(substring $nd 1 4) 00:00"
        nd=$(newdate -D $idate +1)
        StepEnT="$(substring $nd 7 8)/$(substring $nd 5 6)/$(substring $nd 1 4) 00:00"
        
        sed \
        -e "s|PathMapsLF|$PathMapsLF|g" \
        -e "s|StartDate|$StartDate|g" \
        -e "s|StepStarT|$StepStarT|g" \
        -e "s|StepEnT|$StepEnT|g" \
        -e "s|PathOutLVP|$PathOutLVP|g" \
        $suite_libdir/lisvapSettingsTemplate.xml > lv.xml # $suite_libdir/lisvap/lisvapSettingsTemplate.xml > lv.xml
        
        python3 $lisvap_new_env/src/lisvap1.py lv.xml # python3 $suite_libdir/lisvap/src/lisvap1.py lv.xml
       
        # Rename outputs and save in the defined directories
        for et_var in e0 eT0 eS0; do
            cdo splitday ${et_var}.nc ${et_var}_
            python3 ${PYMASKDIR}/mask_compress.py -i ${et_var}_$(substring $idate1 7 8).nc -o ${OUTDIR}/${et_var}/${et_var}_${idate1}.nc
        done
        
        # ------------------------- 
        # Rainbombs correction
        # ------------------------- 
        CorScript=/home/ecm8227/Repositories/glofas_forcings/version5/rainbomb_correction/era5_rainbomb_correction.py
        python3 $CorScript -i ${OUTDIR}/tpO_grb/tpO_grb_${idate}.grb -o TP.grb
        # the python script uses a template with random grib fields, so we need to set the correct ones
        grib_set -s step=24,dataDate=$idate,dataTime=0 TP.grb ${OUTDIR}/tp_grb/tp_grb_${idate}.grb
        
        # -------------------------
        # PYG2P corrected precip
        # -------------------------   
        pyg2p -i ${OUTDIR}/tp_grb/tp_grb_${idate}.grb -o . -c $PYG2PTEMPDIR/AFFS_tp.json  # pyg2p to pcraster
        python3 ${PYMASKDIR}/mask_compress.py -i tp_None.nc -o ${OUTDIR}/tp/tp_${idate1}.nc
        
    else
        echo "The data for $idate have already been generated."
    fi    
        
    idate=$(newdate -D $idate +1)
done

# -------------------------
# Make 1 file per variable
# -------------------------
echo "Concatenating all data into 1 single netcdf/grib."
for param in tpO tp ta e0 eT0 eS0; do
    cdo -O -z zip_6 mergetime ${OUTDIR}/$param/* ${OUTDIR}/../${param}/glofas_${param}_${YYYY}.nc
    cp ${OUTDIR}/../${param}/glofas_${param}_${YYYY}.nc ${OUTDIRPERM}/glofas_${param}_${YYYY}.nc
    
    # add/modify atributes for better compatibility with more softwares
    if [[ $param == 'tpO' ]]
    then
      param_name=tp
    else
      param_name=$param
    fi
    if [[ $param == 'ta' ]]
    then
      param_unit=C
    else
      param_unit=mm
    fi
    python3 ${PYMASKDIR}/modify_attributes.py -y $YYYY -v $param_name -u $param_unit -i ${OUTDIRPERM}/glofas_${param}_${YYYY}.nc
    echo "Netcdf for $param completed."
done

for param in tp_grb tpO_grb; do
    cdo -O mergetime ${OUTDIR}/${param}/* ${OUTDIR}/../${param}/glofas_${param}_${YYYY}.grb  # no zip possible in grb
    mkdir -p ${OUTDIRPERM}/${param}
    cp ${OUTDIR}/../${param}/glofas_${param}_${YYYY}.grb ${OUTDIRPERM}/${param}/glofas_${param}_${YYYY}.grb
    echo "Grib for $param completed."
done

echo "Data generation for year YYYY_input completed."