#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"

#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

 chemistry="GEOSchem" #GEOSchem # MOZ4 # SAPRC99
    srsInp="epsg:4326"	#Los archivos FINN vienen en latlon.

#Parametros de Grilla y Proyección: 
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" #datum=WGS84 +units=m +x_0=0.0 +y_0=0.0 +no_defs"
read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${srsOut}" <<< $( echo "${ref_lon} ${ref_lat}" ) )

#Grilla:
nx=197;ny=237;dx=20000;dy=20000;                #ncols,nrows,xcel, ycell
xorig=-1970000; yorig=-2370000 ;
xmin=$(bc -l <<<" $xc - ($xorig)*-1 ")
xmax=$(bc -l <<<" $xc + ($xorig)*-1 ")
ymin=$(bc -l <<<" $yc - ($yorig)*-1 ")
ymax=$(bc -l <<<" $yc + ($yorig)*-1 ")


#Creo carpetas:
mkdir finn_data		#carpeta de datos descargados de FINN.

#---------------------------------------
#Get FINN data:
  thisYear=$( date +"%Y")

start_date_s=$(date -d "$start_date" +%s)
end_date_s=$(date -d "$end_date" +%s)

day=$start_date_s
while [ $day -le $end_date_s ]
do
	YYYY=$(date -d @$day +"%Y")   #año
	  MM=$(date -d @$day +"%m")   #mes
	 DDD=$(date -d @$day +"%j")   #dia juliano
	echo "Dia: $YYYY-$DDD (mes: $MM)"

	 if [ ! -f finn_data/GLOB_${chemistry}_$YYYY$DDD.txt ]
	 then
	 	if [ $YYYY -lt $(($thisYear-2)) ]	#los años viejos los tienen en carpetas.
	 	then
	 	       wget "https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/$YYYY/GLOB_${chemistry}_$YYYY$DDD.txt.gz" -P finn_data/
	 	else
	 	       wget "https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_${chemistry}_$YYYY$DDD.txt.gz" -P finn_data/
	 	fi;
		
		gzip -d finn_data/GLOB_${chemistry}_$YYYY$DDD.txt.gz
	fi;

	finnFile=finn_data/GLOB_${chemistry}_$YYYY$DDD.txt
	#Obtengo lista de polluts a inventariar
	polluts=($(head -n1 $finnFile | sed 's/,/ /g;s/^.*AREA //g' ))

        swap_list=$(printf "%-16s" ${polluts[@]})
        var_list="$swap_list"

	#Crear template de NetCDF 
cat  > netcdf_emission_template.cdl <<EOF
netcdf emissionInventory {
dimensions:
    TSTEP = 24;
    DATE_TIME = 2 ;
    COL = $nx ;
    ROW = $ny ;
    LAY = 1 ;
    VAR= ${#polluts[@]};

variables:
    int TFLAG(TSTEP, VAR, DATE_TIME) ;
    float Band1(ROW, COL);
    //float emission(time, y, x);
    //Mis variables:

// global attributes:
 :IOAPI_VERSION = "ioapi-3.2: \$Id: init3" ; :EXEC_ID = "???????????????? " ; 
 :FTYPE = 1 ;
 :SDATE = $YYYY$DDD ; :STIME = 000000 ; :WDATE = 2020001 ; :WTIME = 000000 ; :CDATE = 2020000 ; :CTIME = 000000;
 :TSTEP = 10000      ; :NTHIK = 1         ;    
 :NCOLS = $nx        ;  :NROWS = $ny      ;  :NLAYS = $nz    ; :NVARS = ${#polluts[@]} ; :GDTYP = 2 ;
 :P_ALP = -50. ; :P_BET = -20. ; :P_GAM = -65. ; 
 :XCENT = $clon ; :YCENT = $clat ; :XORIG = $xorig ; :YORIG = $yorig ; :XCELL = $dx ;  :YCELL = $dy ;  
 :VGTYP = -9999 ;        :VGTOP = 5000.f ;       :VGLVLS = 1.f, 0.9938147f ;    
 :GDNAM = "$GridName" ;  :UPNAM = "OUTCM3IO        " ;   :VAR-LIST = "${var_list}";      
 :FILEDESC = "Merged emissions output file from Mrggrid" ; :HISTORY = "" ;
}
EOF
#Agrego las variables a inventariar:
for pollut in ${polluts[@]}
do
                sed -i "/\/\/Mis variables:/a     float $pollut(TSTEP, LAY, ROW, COL) ;" netcdf_emission_template.cdl
done;


file_out=fire_emis_${YYYY}${DDD}_d01.nc
ncgen -o $file_out netcdf_emission_template.cdl

hours=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
for HH in ${hours[@]}
do
	echo "   HORA: $HH"
	i=$(( $HH % 24 ));                      #indice de variable time

for j in ${!polluts[@]}
do
	pollut=${polluts[$j]}
	echo "      pollut: $pollut"
#---------------------------------------
#Parse FINN files:

#Filtro columnas lat,lon y pollut, y solo las filas que el date caiga dentro de la hora. tambien convierto emis a mole/s ó g/s.
cat $finnFile | sed 's/D\([-+]\)/e\1/g' | awk -F"," -v pollut=${pollut} -v hour=${HH}00 '
NR==1{
	for(i=1;i<=NF;i++){if($i==pollut){pol=i};};
	print "wkt;emis";

	if ( pollut == "OC" || pollut =="BC" || pollut=="PM25" || pollut == "PM10" ){
		k=1/86400000.0;		#kg/day -> g/s	
	}else {
		k=1/86400.0;		#mole/day -> mole/s
	};
} 
NR>1{
if($2>=hour && $2<(hour+60))  {
	printf("POINT(%.4f %.5f); %.5e\n", $5,$4, $pol/k);
};
}'>tmp.csv

ogr2ogr -s_srs "$srsInp" -t_srs "$srsOut" -f CSV tmp_xy.csv tmp.csv -lco GEOMETRY=AS_WKT -lco SEPARATOR="SEMICOLON"

#---------------------------------------
#ASCII to netcdf:

#gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" -ot Float32 tmp_xy.csv swap.nc -co WRITE_LONLAT="YES"  -co WRITE_BOTTOMUP="NO"
gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of GTiff tmp_xy.csv tmp.tif
#Transformar a NetCDF
gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" tmp.tif swap.nc -co "WRITE_LONLAT=YES" -co "WRITE_BOTTOMUP=NO"

ncks  -A -V -v Band1 swap.nc -o $file_out
ncap2 -A -s "${pollut}($i,0,:,:) = Band1(:,:); TFLAG($i,$j,0) = ${YYYY}${DDD}; TFLAG($i,$j,1) =  ${HH}0000;" $file_out $file_out

done;
done;
	day=$(($day + 86400 )) # Incrementar day un dia
done;
