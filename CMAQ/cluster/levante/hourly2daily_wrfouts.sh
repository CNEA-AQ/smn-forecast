#!/bin/sh
#
# Extracts variables from NetCDF
#
module load nco

set -e
set -u

#inp_files="/work/bm1234/m300788/WRFI/runs/papila_jan2019/wrf/wrfout_d01_2018*"       # silam's output netcdf files directory.
inp_files="/scratch/m/m300788/wrfout/2019/jan/wrfout_d01_"
outdir="wrfouts_horarios"
outdir2=wrfouts
#Horarios a diarios:

start_date="2018-12-25"
end_date="2018-12-31" 
start_date_s=$(date -d "$start_date" +%s)
end_date_s=$(date -d "$end_date" +%s)

if [ ! -d ${outdir2} ]; then mkdir $outdir2; fi;

day=$start_date_s
while [ $day -le $end_date_s ]
do
    date=$(date -d @$day +"%Y-%m-%d")   #año

    echo "Procesando.. ${date}"
    if [ ! -f ${outdir2}/wrfout_d01_${date}.nc ]
    then
        files=($(ls ${outdir}/wrfout_d01_${date}*))
        echo "Archivos: ${files[@]}"
        ncrcat -4 -O -o ${outdir2}/wrfout_d01_${date}.nc ${outdir}/wrfout_d01_${date}*
        ncatted -O -h -a TITLE,global,m,c," OUTPUT FROM WRF V4    MODEL" ${outdir2}/wrfout_d01_${date}.nc
    else
        echo "File ${outdir2}/wrfout_d01_${date}.nc already exists!"
    fi;

    day=$(($day + 86400 )) # Incrementar day un dia
done
