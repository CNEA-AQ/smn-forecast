#!/bin/bash
#
# Objetivo: Completar serie temporal de archivos de mediciones.
#
inp_dir="datafiles"
out_dir="datafile_complete"

start_date="2022-01-01 00:00"
  end_date="2023-01-01 00:00"

year=`date -d "$start_date" +%Y`

#Dates file:
s_ini=`date -d "$start_date" +"%s"`
s_end=`date -d "$end_date" +"%s"`
total_hours=$(bc <<<"($s_end - $s_ini)/(60*60)")

if [ -f date_${year}.csv ]
then 
   continue #rm date_${year}.csv;
else
   echo "date" > date_${year}.csv
   for h in $(seq 0 ${total_hours})
   do
       echo `date -d "+${h} hours ${start_date}" +"%Y/%m/%d %H:00:00"` >> date_${year}.csv
   done
fi

files=($(ls -d $inp_dir/*))
echo ${files[@]}
for f in ${files[@]}
do

   read n1 file<<< $(wc -l date_${year}.csv)
   read n2 file<<< $(wc -l $f)
   
   if [ $n1 -ne $n2 ]
   then
      echo "file: $f : n1: $n1, n2:$n2"  
      awk -F"," 'BEGIN {OFS=","}
           NR == FNR {
             if (NR == 1){ print "date,conc"; next}
             datetime = $1; value = $2
             values[datetime] = value
             next
           }
           {
             current_datetime = $1
             print current_datetime,values[current_datetime]
           }' $f date_${year}.csv > swap.csv
      mv swap.csv $f
   fi;
done
