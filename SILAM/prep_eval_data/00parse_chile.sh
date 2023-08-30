#!/bin/bash
#
# Objetivo: Parser para datos de red de AQ de Chile
#    input_files --> out_dir/<stat>_<year>_<pollut>.csv
#                    date (%Y/%m/%d %H:%M:%H); conc (ug/m3)
#
inp_dir=orig_data/chile
out_dir=datafiles

polluts=(CO NO NO2 SO2 PM10 PM25 O3)
   unit=("ppm" "ppb" "ppb" "ppb" "ugm3" "ugm3" "ppb")
confact=(1150 1.23 1.88 2.62 1.0 1.0 1.96) #factores de conversion ppm/ppb -> ug/m3

stations=(BOMBEROS CENTRO CONCON ELBOSQUE EMEF LAFLORIDA LAPALMA PUDAHUEL PUENTEALTO QUINTERO RANCAGUA1 SANPEDRO SUPERSITE LASCONDES OHIGGINS PUCHUNCAVI POLIVALENTE INDURA KCOLLEGE ENAPPRICE PUNTERAS HUALAQUI ESCUADRON LAGUNILLAS CORONELNORTE CORONELSUR LOTAURBANA LOTARURAL LAJA LAUTARO LASCASAS2 VALDIVIA COYHAIQUE2)
#stations=(BOMBEROS CENTRO CONCON CORONELNORTE CORONELSUR ELBOSQUE ESCUADRON HUALAQUI KCOLLEGE LAFLORIDA LAGUNILLAS LAJA LAPALMA LASCASAS2 LAUTARO LOTARURAL POLIVALENTE PUCHUNCAVI PUDAHUEL PUENTEALTO PUNTERAS QUINTERO RANCAGUA1 SANPEDRO SUPERSITE VALDIVIA)

start_date="2022-01-01"
end_date="2023-01-01"

year=`date -d "$start_date" +%Y`

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

for s in ${stations[@]}
do
  for i in $( seq 1 ${#polluts[@]} )
  do
        p=${polluts[$i-1]}
        f=${confact[$i-1]}
        u=${unit[$i-1]}

        iFile=${inp_dir}/CHILE_${s}_${p}_${u}.csv
        if [ -f $iFile ]
	then 
		echo "Procesando $iFile ..  (pollut: ${p}. cf: ${f})"
        	                                                                                                                             
        	sed 's/\,/\./g' $iFile > tmp_file
        	                                                                                                                             
        	awk -F";" -v pol=$p -v cf=$f -v year=$year '
        	NR==1{print "date,conc"}
        	NR>1{ 
		      yr=substr($1,1,2)+2000;mo=substr($1,3,2);dy=substr($1,5,2);
		      if( yr == year ){
		        hr=substr($2,1,2);mi=substr($2,3,2);
		      	pollut=$3+$4+$5;
        	      	if ( pollut < 0){pollut=NaN }else {pollut=pollut*cf};
        	      	printf("%4d/%02d/%02d %02d:%02d:00,%.2f\n",yr,mo,dy,hr,mi,pollut);
		      }
        	}' tmp_file > $out_dir/${s}_${year}_${p}.csv

	else
		echo "$iFile doesn't exists."
	fi;
  done
  rm tmp_file
done
