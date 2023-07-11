#!/bin/bash

testname=$1

RED='\033[0;31m';GREEN='\033[0;32m';CYAN='\033[0;36m';NC='\033[0m' # No Color

src_file=gis/points_src.csv
rec_file=gis/points_rec.csv

srcfile_fwd=src_point_fwd.v5
srcfile_bwd=src_point_bwd.v5
ctlfile_fwd=control_fwd.ini
ctlfile_bwd=control_bwd.ini

exe=silam_v5_8pub.gnu

rm conc_rec.txt conc_src.txt #Clean
rm $srcfile_fwd $srcfile_bwd

for i in ` seq 1 1 10 `
do
	
	read srcid srclon srclat <<<$(awk -v row=$i -F'[;,() ]*' 'NR==row+1{print $1,$3,$4}' $src_file )
	echo  -e " SOURCE: ${RED}   $srcid $srclon $srclat ${NC}."

cat << EOF >> $srcfile_fwd
POINT_SOURCE_5 # First point-source starts
   source_name = ${srcid}
   source_sector_name =         # source sector name, e.g. SNAP_10. May be empty
   
   source_longitude = ${srclon}
   source_latitude = ${srclat}
   
   plume_rise = NO
   release_rate_unit = kg/sec
   
   vertical_unit = m
   vertical_distribution = SINGLE_LEVEL_DYNAMIC # SINGLE_LEVEL_DYNAMIC, MULTI_LEVEL_FIXED, PLUME_RISE
   stack_height = 10 m
   
   par_str_point = 2017 01 01 06 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   par_str_point = 2017 01 01 07 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   
   hour_in_day_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.  1. 1. 1.
   day_in_week_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1.
   month_in_year_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.
END_POINT_SOURCE_5

EOF
		
done

echo -e " * RUNING  ${RED} FORWARD.${NC}"
echo ""
./${exe} ${ctlfile_fwd} &> fwd_run.log


for j in ` seq 1 1 10 `
do
	read recid reclon reclat <<<$(awk -v row=$j -F'[;,() ]*' 'NR==row+1{print $1,$3,$4}' $rec_file )
	echo  -e " RECEPT: ${GREEN} $recid $reclon $reclat ${NC}."

cat << EOF >> $srcfile_bwd
POINT_SOURCE_5 # First point-source starts
   source_name = ${recid}
   source_sector_name =         # source sector name, e.g. SNAP_10. May be empty
   
   source_longitude = ${reclon}
   source_latitude = ${reclat}
   
   plume_rise = NO
   release_rate_unit = kg/sec
   
   vertical_unit = m
   vertical_distribution = SINGLE_LEVEL_DYNAMIC # SINGLE_LEVEL_DYNAMIC, MULTI_LEVEL_FIXED, PLUME_RISE
   stack_height = 10 m
   
   par_str_point = 2017 01 01 15 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   par_str_point = 2017 01 01 16 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   
   hour_in_day_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.  1. 1. 1.
   day_in_week_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1.
   month_in_year_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.
END_POINT_SOURCE_5

EOF
done

echo -e " * RUNING  ${GREEN} BACKWARDS ${NC}."	
echo ""

./${exe} ${ctlfile_bwd} &> bwd_run.log


#=============================================
#COMPARE RESULTS:

rm conc_rec.txt conc_src.txt

for i in ` seq 1 1 10 ` #for each source
do
        read srcid srclon srclat <<<$(tail -n+3 $src_file | head -n $i | tail -n 1 | awk -F'[;,() ]*' '{print $1,$3,$4}')
        fwd_file_name="output/fwd_${srcid}.nc"

        for j in ` seq 1 1 10 ` #for each receptor
        do
                read recid reclon reclat <<<$(tail -n+3 $rec_file | head -n $j | tail -n 1 | awk -F'[;,() ]*' '{print $1,$3,$4}')
                bwd_file_name="output/bwd_${recid}.nc"

                #pick value at rec
                ncks -v cnc_passive_gas -d lat,$reclat -d lon,$reclon -d time,10  -d height,12.5 $fwd_file_name | sed -n "/cnc_passive_gas =/{n;p}" >> conc_rec.txt
                #pick value at src
                ncks -v cnc_passive_gas -d lat,$srclat -d lon,$srclon -d time,10 -d height,12.5 $bwd_file_name | sed -n "/cnc_passive_gas =/{n;p}" >> conc_src.txt
        done
done

paste conc_* > ${testname}.txt

cat << EOF | gnuplot
set datafile separator ";"
set logscale xy
plot '${testname}.txt' using 1:2 with points
#set xrange [0:1E-9]
#set yrange [0:1E-9]
pause
EOF

