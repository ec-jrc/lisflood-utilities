Y_strt=1975
# Y_strt=2024 # only for checking
Y_stop=2024

mkdir -p yearly_scripts
mkdir -p yearly_logs

Y_iter=$Y_strt
while [[ $Y_iter -le $Y_stop ]]; do
    sed \
        -e "s|YYYY_input|$Y_iter|g" \
       GloFAS_v5_ERA5_forcings_gen_template.ksh > yearly_scripts/data_gen_${Y_iter}.ksh

    ecsbatch yearly_scripts/data_gen_${Y_iter}.ksh

    Y_iter=$(($Y_iter+1))

done