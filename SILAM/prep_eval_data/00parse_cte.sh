#!/bin/bash
#Parser para datos de CTE:
input_files=(EMCABB1_2022.csv EMCABB2_2022.csv)
out_dir=datafiles

polluts=(CO O3 SO2 PM10 NO2 NO NOx)
confact=(1150 1.96 2.62 1 1.88 1.23 1.88) #factores de conversion ppm/ppb -> ug/m3

stations=(EMCABB1 EMCABB2)
start_date="2022-01-01"
  end_date="2023-01-01"

year=`date -d "$start_date" +%Y`

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

for s in ${stations[@]}
do
  iFile=${s}_${year}.csv
  echo "Procesando $iFile .."
  sed 's/ *---/-999/g; s/< LD/0.0/g' $iFile > tmp_file

  for i in $( seq 1 ${#polluts[@]} )
  do
	p=${polluts[$i-1]}
	f=${confact[$i-1]}
	echo "Pollut: ${p}. cf: ${f} .."
	
        awk -F";" -v pol=$p -v cf=$f -v year=$year '
        NR==1{
	   for(i=1; i<=NF; i++) {
	      if ($i == pol ){ pollut=i; found=1 };
	   }
        if (found) {printf("%s,%s\n","date",$pollut)}else{exit 1};
        }
        NR>1{ if ( $pollut < 0){ $pollut=NaN}else {$pollut=$pollut*cf};
	      split($1, date, /[\/: ]/);
	      if ( year == date[3] ){ printf("%4d/%02d/%02d %02d:%02d:00,%.2f\n",date[3],date[2],date[1],date[4],date[5],$pollut) };
        }' tmp_file > $out_dir/${s}_${year}_${p}.csv
  done
  rm tmp_file
done
