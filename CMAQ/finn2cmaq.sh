#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

chemistry="GEOSchem" #GEOSchem # MOZ4 # SAPRC99
   srsInp="epsg:4326"	#Los archivos FINN vienen en latlon.

#Parametros de Grilla y Proyección: 
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 
read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${srsOut}" <<< $( echo "${ref_lon} ${ref_lat}" ) )

#Grilla (sale de GRIDDESC)
nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell
xorig=-1970000;yorig=-2370000;

#Ciclo Diurno:
#fire_emis (hora local) usa:#ciclo_diurno=(0.43 0.43 0.43 0.43 0.43 0.43 0.43 0.43 0.43 3.0 60. 10.0 14.0 17.0 14.0 12.0 9.0 6.0 3.0 0.43 0.43 0.43 0.43 0.43) # (en porcentaje)
ciclo_diurno=(0.0043 0.0043 0.0043 0.0043 0.0043 0.0043 0.0043 0.0043 0.0043 0.0300 0.0600 0.1000 0.1400 0.1700 0.1400 0.1200 0.0900 0.0600 0.0300 0.0043 0.0043 0.0043 0.0043 0.0043)

#Perfil vertical (Plume Rise)  (!) To-Do.


#------------------------------------------------
#Parametros intermedios:
xmin=$(bc -l <<<" $xc - ($xorig)*-1 ")
xmax=$(bc -l <<<" $xc + ($xorig)*-1 ")
ymin=$(bc -l <<<" $yc - ($yorig)*-1 ")
ymax=$(bc -l <<<" $yc + ($yorig)*-1 ")

#bbox latlon: para hacer clippear el archivo de fuegos al dominio:
read lonmin latmin ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmin} ${ymin}" ) )
read lonmax latmax ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmax} ${ymax}" ) )

start_date_s=$(date -d "$start_date" +%s)
end_date_s=$(date -d "$end_date" +%s)     
YYYY=$(date -d @$start_date_s +"%Y")   #año
 DDD=$(date -d @$start_date_s +"%j")   #dia juliano

if [[ $chemistry == "GEOSchem" ]]; then
	polluts=(CO2 CO NO NO2 SO2 NH3 CH4 ACET ALD2 ALK4 BENZ C2H2 C2H4 C2H6 C3H8 CH2O GLYC GLYX HAC MEK MGLY PRPE TOLU XYLE OC BC PM25)
elif [[ $chemistry == "MOZ4" ]]; then
	polluts=(CO2 CO NO NO2 SO2 NH3 CH4 VOC ACET ALK1 ALK2 ALK3 ALK4 ALK5 ARO1 ARO2 BALD CCHO CCO_OH ETHENE HCHO HCN HCOOH HONO ISOPRENE MEK MEOH METHACRO MGLY MVK OLE1 OLE2 PHEN PROD2 RCHO TRP1 OC BC PM10 PM25)
elif [[ $chemistry == "SAPRC99" ]]; then
	polluts=(CO2 CO H2 NO NO2 SO2 NH3 CH4 NMOC BIGALD BIGALK BIGENE C10H16 C2H4 C2H5OH C2H6 C3H6 C3H8 CH2O CH3CHO CH3CN CH3COCH3 CH3COCHO CH3COOH CH3OH CRESOL GLYALD HCN HYAC ISOP MACR MEK MVK TOLUENE HCOOH C2H2 OC BC PM10 PM25)
else
	echo "No existe la química seleccionada: chemistry=${chemistry}. (Opciones validas: GEOSchem, MOZ4, SAPRC99)"; exit;
fi;
var_list=$(printf "%-16s" ${polluts[@]})

#Crear estructura del NetCDF de salida
cat  > netcdf_emission_template.cdl <<EOF
netcdf emissionInventory {
dimensions:
    TSTEP = 1; // 24;
    DATE_TIME = 2 ;
    COL = $nx ;
    ROW = $ny ;
    LAY = 1 ;
    VAR= ${#polluts[@]};

variables:
    int TFLAG(TSTEP, VAR, DATE_TIME) ;
    float Band1(ROW, COL);
EOF
printf "     float %s(TSTEP, LAY, ROW, COL);\n" ${polluts[@]} >> netcdf_emission_template.cdl #   //Mis variables:
cat >> netcdf_emission_template.cdl << EOF
// global attributes:
 :IOAPI_VERSION = "ioapi-3.2: \$Id: init3" ; :EXEC_ID = "???????????????? " ; 
 :FTYPE = 1 ; :SDATE = ${YYYY}${DDD} ; :STIME = 000000 ; :WDATE = 2023001 ; :WTIME = 000000 ; :CDATE = 2023001 ; :CTIME = 000000;
 :TSTEP = 10000; :NTHIK = 1;    
 :NCOLS = ${nx}; :NROWS = ${ny}; :NLAYS = ${nz}; :NVARS = ${#polluts[@]}; :GDTYP = 2;
 :P_ALP = -50.;:P_BET = -20. ;:P_GAM = -65.; 
 :XCENT = ${ref_lon}; :YCENT = ${ref_lat}; :XORIG = ${xorig}; :YORIG = ${yorig}; :XCELL = ${dx};  :YCELL = ${dy};  
 :VGTYP = -9999 ;:VGTOP = 5000.f ;:VGLVLS = 1.f, 0.9938147f;    
 :GDNAM = "${GridName}";:UPNAM = "OUTCM3IO" ;   :VAR-LIST = "${var_list}";      
 :FILEDESC = "Merged emissions output file from Mrggrid" ; :HISTORY = "" ;
}
EOF
#---------------------------------------
mkdir finn_data		#carpeta de datos descargados de FINN.
thisYear=$( date +"%Y")

day=$start_date_s
while [ $day -le $end_date_s ]
do
	YYYY=$(date -d @$day +"%Y")   #año
	  MM=$(date -d @$day +"%m")   #mes
	 DDD=$(date -d @$day +"%j")   #dia juliano
	  HH=$(date -d @$day +"%H")   #hora
	echo "Dia: $YYYY-$DDD (mes: $MM)"

	#=================
	#Get FINN data:
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
	#polluts=($(head -n1 $finnFile | sed 's/,/ /g;s/^.*AREA //g' )); 

	#=================
	#Transformo coordenadas + me quedo con los puntos dentro del dominio.
	cat $finnFile | sed 's/D\([-+]\)/e\1/g;s/,/;/g' | awk -F";" 'NR==1{print $0";wkt"} NR>=2{print $0"POINT("$5*1.0,$4*1.0")"}' > tmp.csv
	ogr2ogr -f CSV tmp_ll_clipped.csv tmp.csv -clipsrc $lonmin $latmin $lonmax $latmax -lco GEOMETRY="AS_WKT" -lco SEPARATOR="SEMICOLON"
	ogr2ogr -f CSV tmp_xy.csv tmp_ll_clipped.csv -s_srs "$srsInp" -t_srs "$srsOut" -lco GEOMETRY="AS_WKT" -lco SEPARATOR="SEMICOLON"
	sed -i 's/"//g;s/  */ /g;' tmp_xy.csv
	#=================
	#Un archivo por hora:
		
	
	
	
	for HH in $(seq --format="%02.0f" 0 23) 
	do
		echo "   HORA: $HH"
		i="$(printf "%d" ${HH})";
	
		#Peso horario:
		wgt_hh=${ciclo_diurno[i]}

		#Creo NetCDF destino:
		file_out=fire_emis_${YYYY}${DDD}_${HH}:00:00_d01.nc
		ncgen -o $file_out netcdf_emission_template.cdl
	
		#---------------------------------------
		for j in ${!polluts[@]}
		do
			pollut=${polluts[$j]}
			echo "      pollut: $pollut"
			#---------------------------------------
			#Filtro columnas wkt y pollut, tambien convierto emis a mole/s ó g/s (aplicando ciclo diurno).
			awk -F";" -v pollut=${pollut} -v wgt_hh=${wgt_hh} '
			NR==1{
			for(i=1;i<=NF;i++){if($i==pollut){pol=i};};
				print "wkt;emis";
				if ( pollut == "OC" || pollut =="BC" || pollut=="PM25" || pollut == "PM10" ){
					k=1/3600000*wgt_hh;	#kg/day -> g/s	
				}else {
					k=1/3600.0 *wgt_hh;	#mole/day -> mole/s
				};
			} 
			NR>1{ printf("%s;%.5e\n", $1, $pol*k); }' tmp_xy.csv > tmp_pollut.csv
			#---------------------------------------
			##ASCII to netcdf
			#gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" -ot Float32 tmp_pollut.csv tmp.nc
			gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of GTiff tmp_pollut.csv tmp.tif
			gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" tmp.tif tmp.nc #-co "WRITE_LONLAT=YES" -co "WRITE_BOTTOMUP=NO"
	
			ncks  -h -A -V -v Band1 tmp.nc -o $file_out

			ncap2 -h -A -C -s "${pollut}(0,0,:,:) = Band1(:,:); TFLAG(0,$j,0) = ${YYYY}${DDD}; TFLAG(0,$j,1) =  ${HH}0000;" $file_out 

		done;
		
	done;
	day=$(($day + 86400 )) # Incrementar day un dia
done;

