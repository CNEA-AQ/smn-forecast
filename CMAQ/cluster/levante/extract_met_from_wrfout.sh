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

#start_date="2018-12-25" #_00:00:00" #"2018-12-25_00:00:00"
start_date="2019-01-02"
end_date="2019-02-01"

#variables to extract:
variables=(Times ALBEDO C1F C1H C2F C2H CANWAT CLDFRA DZS GLW HFX HGT ISLTYP LAI LANDMASK LH LU_INDEX MAPFAC_M MAPFAC_U MAPFAC_V MU MUB P PB PBLH PH PHB PSFC P_TOP Q2 QCLOUD QFX QGRAUP QICE QRAIN QSNOW QVAPOR RAINC RAINNC SEAICE SMOIS SNOWC SNOWH SWDOWN T T2 TSK TSLB U U10 UST V V10 VEGFRA W XLAT XLAT_U XLAT_V XLONG XLONG_U XLONG_V ZNU ZNW)

echo "" > scriptfile.tmp
for v in ${variables[@]}
do 
echo "${v}=${v};" >> scriptfile.tmp
done
cat scriptfile.tmp
if [ ! -d ${outdir} ]; then mkdir $outdir; fi;

start_date_s=$(date -d "$start_date" +%s)
end_date_s=$(date -d "$end_date" +%s)

day=$start_date_s
while [ $day -le $end_date_s ]
do
    date=$(date -d @$day +"%Y-%m-%d")   #año
    files=($(ls -d ${inp_files}${date}*)) 

    echo "Procesando.. ${date}"
    
    pairs_inout=""
    for f in ${files[@]}
    do
    	of="${outdir}/$(basename $f)"
    	pairs_inout="$pairs_inout $f $of"
    done
    echo "${pairs_inout[@]}"
    
    export OMP_NUM_THREADS=28 #${#files[@]}
    
    if [ -n "$pairs_inout" ]; then
        echo -n $pairs_inout | xargs -r -n 2 echo  ncap2 -4 -L5 -O -v -S scriptfile.tmp | srun --input=0 -J post --partition=compute --mail-type=FAIL --mail-user=pablo.lichtig@mpimet.mpg.de --account=bm1234 --nodes=1 --ntasks-per-node=28 --time 00:10:00 xargs -IXXX -t -P28 sh -c "XXX"
 #       echo -n $pairs_inout | xargs -r -n 2 echo  ncap2 -4 -L5 -O -v -S scriptfile.tmp | srun --input=0 -J post --partition=compute --mail-type=FAIL --mail-user=pablo.lichtig@mpimet.mpg.de --account=bm1234 --nodes=1 --ntasks-per-node=128 --time 00:10:00 xargs -IXXX -t -P128 sh -c "XXX"
    fi
    
    day=$(($day + 86400 )) # Incrementar day un dia
done;
echo Done!
