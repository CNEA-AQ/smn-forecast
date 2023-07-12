#!/bin/bash

testname=$1
runmode=$2

RED='\033[0;31m';GREEN='\033[0;32m';CYAN='\033[0;36m';NC='\033[0m' # No Color

#Inputs:
ctlfile=control_template.ini   #=control_fwd_norway.ini

src_file=gis/points_src_argentina.csv	#gis/points_src_argentina.csv
rec_file=gis/points_rec_argentina.csv   #gis/points_rec_argentina.csv

start_date="2009/03/15 18:00:00" # 2017/01/01 07:00:00" #2009 03 15 18 00 0.0 
computed_period="10 hours"

bbox=(-65. -50. -55. -40.)  #(-10. 60.0 0.0 70.0) #(lonmin latmin lonmax latmax)

###
end_date=$(date -d"${computed_period} $start_date" +"%Y/%m/%d %H:%M")

start_date1=$(date -d "1 hour $start_date" +"%Y %m %d %H %M")
  end_date1=$(date -d "-1 hour $end_date" +"%Y %m %d %H %M")


if [ ${runmode} == "-main" ]
then
        srcfile_fwd=src_point_fwd.v5
        srcfile_bwd=src_point_bwd.v5
                                     
        exe=silam_v5_8pub.gnu

 	#Setup controlfile:
	sed -i "s/ case_name = .*/ case_name = ${testname}_/" ${ctlfile}

        sed -i "s/ computed_period =.*/ computed_period = ${computed_period//hours/hr/}"${ctlfile}
                                                                                                       
        sed -i "s/ start_lon = .*/ start_lon = ${bbox[0]}_/" ${ctlfile}
        sed -i "s/ start_lon = .*/ start_lat = ${bbox[1]}_/" ${ctlfile}
        sed -i "s/ end_lon = .*/ end_lon = ${bbox[2]}_/" ${ctlfile}
        sed -i "s/ end_lat = .*/ end_lat = ${bbox[3]}_/" ${ctlfile}
	
	
	if [ -f $srcfile_fwd ]; then rm $srcfile_fwd;fi
	if [ -f $srcfile_bwd ]; then rm $srcfile_bwd;fi
	



	#FORWARDS:
	echo -e " * RUNING  ${RED} FORWARD${NC}."
	#echo ""
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
   
   !par_str_point = 2017 01 01 06 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   !par_str_point = 2017 01 01 07 00 0.0    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   par_str_point = ${start_date//[-:\/]/ }    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.  
   par_str_point = ${start_date1//[-:\/]/ }    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.
   
   hour_in_day_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.  1. 1. 1.
   day_in_week_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1.
   month_in_year_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.
END_POINT_SOURCE_5

EOF
	done
	
	#Changing control file parameters:
	sed -i "s/ direction_in_time =.*/ direction_in_time = FORWARD/" ${ctlfile}
	sed -i "s/ start_time =.*/ start_time = ${start_date//[-:\/]/ }"${ctlfile}

	./${exe} ${ctlfile} &> fwd_run.log



	#BACKWARDS:
	echo -e " * RUNING  ${GREEN} BACKWARDS ${NC}. $testname"	
        #echo ""
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
   
   par_str_point = ${end_date1//[:-]/}    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.  
   par_str_point = ${end_date//[:-]/}    1.    500. 1000.   5.0   450.  PASSIVE_COCKTAIL 25.

   hour_in_day_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.  1. 1. 1.
   day_in_week_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1.
   month_in_year_index = PASSIVE_COCKTAIL 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.
END_POINT_SOURCE_5

EOF
	done
	
        #Changing control file parameters:
        sed -i "s/ direction_in_time =.*/ direction_in_time = BACKWARD/" ${ctlfile}
        sed -i "s/ start_time =.*/ start_time = ${end_date//[-:\/]/ }"${ctlfile}


	./${exe} ${ctlfile} &> bwd_run.log

	#=============================================
	#COMPARE RESULTS:
	
	if [ -f conc_rec.txt ]; then rm conc_rec.txt;fi
        if [ -f conc_src.txt ]; then rm conc_src.txt;fi
	
	for i in ` seq 1 1 10 ` #for each source
	do
	        read srcid srclon srclat <<<$(tail -n+3 $src_file | head -n $i | tail -n 1 | awk -F'[;,() ]*' '{print $1,$3,$4}')
	        fwd_file_name="output/${testname}_${srcid}.nc"
	
	        for j in ` seq 1 1 10 ` #for each receptor
	        do
			#to track differences:
			echo "$recid;$srcid" >> ids.txt
	
	                read recid reclon reclat <<<$(tail -n+3 $rec_file | head -n $j | tail -n 1 | awk -F'[;,() ]*' '{print $1,$3,$4}')
	                bwd_file_name="output/${testname}_${recid}.nc"
	
	                #pick value at rec
	                ncks -v cnc_passive_gas -d lat,$reclat -d lon,$reclon -d time,10  -d height,12.5 $fwd_file_name | sed -n "/cnc_passive_gas =/{n;p}" >> conc_rec.txt
	                #pick value at src
	                ncks -v cnc_passive_gas -d lat,$srclat -d lon,$srclon -d time,10 -d height,12.5 $bwd_file_name | sed -n "/cnc_passive_gas =/{n;p}" >> conc_src.txt
	        done
	done
	
	paste conc_* > tmp.txt
	paste tmp.txt ids.txt | sed 's/[ \t]*//g' > ${testname}.txt
	
	rm tmp.txt conc_* ids.txt	#clean temp files
	
	# "Error" meassure:
	awk -F";" 'BEGIN{sum=0;N=0}{a1=($1-$2)^2;a2=($1+$2)^2;if(a2>0){sum+=a1/a2;N++}}END{print "Error measure = ",sum/N}' ${testname}.txt

fi


# XY-scatter plot

if [ ${runmode} == "-plot" ]
then
cat << EOF | gnuplot
set datafile separator ";"
set logscale xy
plot '${testname}.txt' using 1:2 with points, '' using 1:2:(sprintf("%s-%s",strcol(3),strcol(4))) with labels
#set xrange [0:1E-9]
#set yrange [0:1E-9]
pause 100

EOF

fi


#####################################
#Time-step sensibility:

if [ ${runmode} == "-dtSens" ]
then
	#dx=dy=0.5
	sed -i "s/dx = .*/dx = 0.5/" ${ctlfile}
	sed -i "s/dx = .*/dx = 0.5/" ${ctlfile}
	sed -i "s/dy = .*/dy = 0.5/" ${ctlfile}
	sed -i "s/dy = .*/dy = 0.5/" ${ctlfile}

	#how to run it:
	dt=(3 5 10 15 25 30 60)
	for time_step in ${dt[@]}
	do
		echo "* dt = $time_step"
		sed -i "s/ time_step = .*/ time_step = ${time_step} min/" ${ctlfile}
		sed -i "s/ time_step = .*/ time_step = ${time_step} min/" ${ctlfile}
		./runAll.sh "test${time_step}min" "-main"
	done
fi

#RESULTS:					
#dx=0.25					#dx=0.5
#| dt     |Error measure|			#| dt     |Error measure|
#|--------|-------------|                       #|--------|-------------|
#| 03min  |  0.453602   | 			#| 03min  |  0.367032   |
#| 05min  |  0.453037   |                       #| 05min  |  0.333031   |
#| 10min  |  0.357329   |  <-- optimal          #| 10min  |  0.296046   |
#| 15min  |  0.390781   |                       #| 15min  |  0.263579   | <-- optimal
#| 25min  |  0.738043   |                       #| 25min  |  0.532827   |
#| 30min  |  0.635152   |                       #| 30min  |  0.359187   |
#| 60min  |  0.727948   |                       #| 60min  |  0.748884   |


#####################################
#Grid resolution sensibility:

if [ ${runmode} == "-dxSens" ]
then
	#dt=15min
	sed -i "s/ time_step = .*/ time_step = 15 min/" ${ctlfile}
	sed -i "s/ time_step = .*/ time_step = 15 min/" ${ctlfile}
	
	resol=(0.05 0.1 0.25 0.5 0.75 1.0)
	for dx in ${resol[@]}
	do
		echo "* dx = $dx"

		sed -i "s/ dx = .*/ dx = ${dx}/" ${ctlfile}
		sed -i "s/ dx = .*/ dx = ${dx}/" ${ctlfile}
		sed -i "s/ dy = .*/ dy = ${dx}/" ${ctlfile}
		sed -i "s/ dy = .*/ dy = ${dx}/" ${ctlfile}
	
		./runAll.sh "test${dx}res" "-main"
	done
fi

#RESULTS:
#dt=10min					#dt=15min
#| dx=dy  | Error measure |			#| dx=dy  | Error measure |
#|--------|---------------|                     #|--------|---------------|
#| 0.05   |  0.397859     |			#| 0.05   |    0.474032   |  
#| 0.1    |  0.486663     |                     #| 0.1    |    0.515637   |
#| 0.25   |  0.357329     |                     #| 0.25   |    0.390781   |
#| 0.5    |  0.296046     | <-- optimal         #| 0.5    |    0.263579   | <- optimal
#| 0.75   |  0.302418     |                     #| 0.75   |    0.282655   |
#| 1.0    |  0.340497     |                     #| 1.0    |    0.329181   |


#####################################
#Diferent place & meteo data:

#Just change ctrl and src filesnames

#Domain: Argentina, meteo: GFS-WRF (20km resol.)

#dx=0.5   					#dx=0.5 (Norway) <= (Winds are stronger)
#| dt      | Error measure|                     #| dt     |Error measure|
#|---------|--------------|                     #|--------|-------------|
#|  03min  |   0.63104    |                     #| 03min  |  0.367032   |
#|  05min  |   0.568971   | 	                #| 05min  |  0.333031   |
#|  10min  |   0.515291   | <-- optimal         #| 10min  |  0.296046   |
#|  15min  |   0.51597    |                     #| 15min  |  0.263579   | <-- optimal
#|  25min  |   0.768848   |                     #| 25min  |  0.532827   |
#|  30min  |   0.711312   |                     #| 30min  |  0.359187   |
#|  60min  |   0.827925   |                     #| 60min  |  0.748884   |


#dt=15min					#dt=15min
#| dx=dy   | Error measure|			#| dx=dy  | Error measure |
#|---------|--------------|                     #|--------|---------------|
#|  0.05   |  0.693682	  |                     #| 0.05   |    0.474032   |  
#|  0.1    |  0.628369	  |                     #| 0.1    |    0.515637   |
#|  0.25   |  0.607903	  |                     #| 0.25   |    0.390781   |
#|  0.5    |  0.51597	  |                     #| 0.5    |    0.263579   | <-- optimal
#|  0.75   |  0.485843	  | <-- optimal         #| 0.75   |    0.282655   |
#|  1.0    |  0.593787	  |                     #| 1.0    |    0.329181   |



