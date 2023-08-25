#!/bin/bash
#Parser para datos de ACUMAR: <station>_<year>_<pollut>.csv
input_file=orig_data/baseDeDatos200723_acumar.xlsx
out_dir=datafiles ; if  [ ! -d $tmp_dir ]; then mkdir $tmp_dir; fi
tmp_dir=tmp_sheets; if [ ! -d $out_dir ]; then mkdir $out_dir; fi

polluts=(CO NO NO2 NOx SO2 PM10 PM25 O3) # SH2)
stations=(EMC1 EMC2)
start_date="2022-01-01"
  end_date="2023-01-01"
year=`date -d "$start_date" +%Y`

#echo "spliting xlsx file: $input_file" 
#xlsx2csv --all --ignoreempty --dateformat "%Y/%m/%d %H:%M:%S" --timeformat "%H:%M:%S" $input_file $tmp_dir

echo "renaming csv files.."
rename 's/ //g;' $tmp_dir/*
rename 's/EMCII/EMC2/g;s/EMCI/EMC1/g;' $tmp_dir/*
rename 's/PM2\,5/PM25/g;s/PM2\.5/PM25/g;' $tmp_dir/*

for s in ${stations[@]} 
do
  echo "${s}"
  for p in ${polluts[@]} 
  do
     echo "  ${p}"
     ls $tmp_dir/* | grep "$s" | grep "\-$p\." | xargs cat | sort > tmp_${s}_${p}.csv
     sed -i 's/µg\/m3/ug\/m3/g' tmp_${s}_${p}.csv
     #Filter by year:
     awk -F"[,]" -v year=$year 'BEGIN{OFS=",";printf("%s,%s\n","date","conc")}
     {
     split($1, date, /[\/]/);
     if ( date[1] == year){
	if($3 == "mg/m3"){ 
	      printf("%s, %f\n",$1,$2*1000) 
	}else{         
              printf("%s, %f\n",$1,$2)
        }
     };
     }' tmp_${s}_${p}.csv > $out_dir/${s}_${year}_${p}.csv
     rm tmp_${s}_${p}.csv 
  done
done;
