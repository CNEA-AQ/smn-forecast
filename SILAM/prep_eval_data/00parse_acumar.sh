#!/bin/bash
#Parser para datos de ACUMAR:
#
#Genero archivos: <station>_<year>_<pollut>.csv
#

input_file=baseDeDatos200723.xlsx	#
tmp_dir=tmp_sheets
out_dir=datafiles

#Split .xlsx into many .csv files (one per each sheet):
echo "split xlsx: $input_file"
#xlsx2csv --all --ignoreempty $input_file $tmp_dir

echo "renaming.."
#Algunos cambios a los archivos (remuevo espacios)
rename 's/ /_/g;' $tmp_dir/*
rename 's/EMC_II/EMC2/g;' $tmp_dir/*
rename 's/EMCII/EMC2/g;' $tmp_dir/*
rename 's/EMC_I/EMC1/g;' $tmp_dir/*
rename 's/PM2\,5/PM25/g;' $tmp_dir/*
rename 's/PM2\.5/PM25/g;' $tmp_dir/*
rename 's/__//g;' $tmp_dir/*
rename 's/_\.csv/\.csv/g;' $tmp_dir/*

polluts=(CO NO NO2 NOx PM10 PM25 SO2 SH2)
stations=(EMC1 EMC2)
start_date="2022-01-01"
  end_date="2023-01-01"

year=`date -d "$start_date" +%Y`

for s in ${stations[@]} 
do
  echo "${s}"
  for p in ${polluts[@]} 
  do
     echo "  ${p}"

     ls $tmp_dir/* | grep "$s" | grep "\-$p\." | xargs cat | sort > tmp_${s}_${p}.csv
     sed -i 's/µg\/m3/ug\/m3/g' tmp_${s}_${p}.csv

     #Filter by year:
     awk -F"[/ :,-]" -v year=$year 'BEGIN{OFS=","}
     {
     if ( "20"$3 == year){
     
	if($7"/"$8 == "mg/m3"){ 
		printf("20%02d/%02d/%02d %02d:%02d, %f\n",$3,$1,$2,$4,$5,$6*1000) 
	}else{                          
                printf("20%02d/%02d/%02d %02d:%02d, %f\n",$3,$1,$2,$4,$5,$6)
        }
     };
     }' tmp_${s}_${p}.csv | sort > $out_dir/${s}_${year}_${p}.csv
     rm tmp_${s}_${p}.csv 

  done
done;


#Dates file:

#s_ini=`date -d "$start_date" +"%s"`
#s_end=`date -d "$end_date" +"%s"`
#total_hours=$(bc <<<"($s_end - $s_ini)/(60*60)")
#
#if [ -f date_file.csv ]; then rm date_file.csv;fi;
#
#for h in $(seq 0 ${total_hours})
#do
#    echo `date -d "+${h} hours ${start_date}" +"%Y-%m-%d %H:00:00"` >> date_file.csv
#done
#
