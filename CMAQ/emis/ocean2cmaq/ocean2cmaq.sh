#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

   srsInp="epsg:4326"	#Los archivos vienen en latlon.

#Input files:
   landpoly=shp/land_10m_global.gpkg	#polygons to make a landmask
#chlor_file=                            #chlorophyl conc. files (from MODIS (MY1DMM)
#  dms_file=		                #dms conc. files (DMS climatology of Lana et al. (2011))
#----------

#Parametros de Grilla y Proyección: 
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 

read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${srsOut}" <<< $( echo "${ref_lon} ${ref_lat}" ) )

#Grilla (sale de GRIDDESC)
nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell
xorig=-1970000;yorig=-2370000;
#------------------------------------------------
#Parametros intermedios:
xmin=$(bc -l <<<" $xc - ($xorig)*-1 ")
xmax=$(bc -l <<<" $xc + ($xorig)*-1 ")
ymin=$(bc -l <<<" $yc - ($yorig)*-1 ")
ymax=$(bc -l <<<" $yc + ($yorig)*-1 ")

#bbox latlon: para hacer clippear el archivo de fuegos al dominio:
read lonmin latmin ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmin} ${ymin}" ) )
read lonmax latmax ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmax} ${ymax}" ) )

#------------------------------------------------
# SURF:
#(1) Camino 1: rasterizo las lineas de costas muy fino (100m o 200m), asigno 1/0 si toca o no costa. Luego regrid a grilla de modelado calculando "average".
ogr2ogr -f "GPKG" test.gpkg $landpoly -nln coast -nlt MULTILINESTRING -t_srs "${srsOut}" -clipsrc $lonmin $latmin $lonmax $latmax #-explodecollections 
gdal_rasterize -q -l coast -at -burn 1 -a_srs "${srsOut}" -tr 500 500 -te $xmin $ymin $xmax $ymax -of GTiff test.gpkg tmp.tif
gdalwarp -r average -t_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax tmp.tif SURF.tif 
rm tmp.tif

###(2) Camino 2: obtengo poligono con borde de costa +100m. Y luego rasterizo (no está resuelto como obtener la fraccion de poligono cubierta en cada pixel).
###Reproject & Cut with BBOX (Importante: make the UNION so we have only one feature)
###ogr2ogr -f "GPKG" test.gpkg $landpoly -nln land   -nlt MULTIPOLYGON -t_srs "${srsOut}" -clipsrc $lonmin $latmin $lonmax $latmax -dialect SQLite -sql "SELECT ST_UNION( Geometry ) AS geom FROM ne_10m_land"
##ogr2ogr -f "GPKG" test.gpkg $landpoly -nln land   -nlt MULTIPOLYGON -t_srs "${srsOut}" -clipsrc $lonmin $latmin $lonmax $latmax
###Buffer 50m from shoreline:
##ogr2ogr -f "GPKG" test.gpkg test.gpkg -nln buffer -nlt MULTIPOLYGON -update -dialect SQLite -sql "SELECT ST_Buffer( land.geom, 100) AS geom FROM land" 
###ogr2ogr -f "GPKG" test.gpkg test.gpkg -nln buffer -nlt MULTIPOLYGON -update -dialect SQLite -sql "SELECT ST_UNION( ST_Buffer( land.geom, 100)  ) AS geom FROM land" 
###Make difference between buffer and land: (to define sea-spray area production).
###ogr2ogr -f "GPKG" test.gpkg test.gpkg -nln differ -nlt MULTIPOLYGON -update -dialect SQLite -sql "SELECT ST_UNION(ST_Difference( buffer.geom, land.geom )) AS geom FROM land,buffer"
##ogr2ogr -f "GPKG" test.gpkg test.gpkg -nln differ -nlt MULTIPOLYGON -update -dialect SQLite -sql "SELECT ST_Difference( buffer.geom, land.geom ) AS geom FROM land,buffer"
###Burn as fraction (FALTA ESTO Y CASI ESTAMOS)
##gdal_rasterize -b 1 -i -burn 0 -l differ -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax test.gpkg output_raster.tif
#gdal_rasterize -q -l differ -at -burn 1 -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of GTiff test.gpkg tmp.tif
#gdal_grid -a count -a_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax -of GTiff tmp2.tif tmp.tif
##gdal_rasterize -a fraction -tr 50 50 -ot Float32 -a_srs "${srsOut}" -init 0 -a_nodata 0 -te $xmin $ymin $xmax $ymax -co COMPRESS=DEFLATE -co ZLEVEL=9 ${coastlines_shp} coastline_mask.tif
#echo "gdal_rasterize -b 1 -burn 50 -tr 50 50 -ot Float32 -a_srs "${srsOut}" -init 0 -a_nodata 0 -te $xmin $ymin $xmax $ymax -co COMPRESS=DEFLATE -co ZLEVEL=9 ${coastlines_shp} coastline_mask.tif"
##gdal_rasterize -b 1 -burn 50 -tr 50 50 -ot Float32 -a_srs "${srsOut}" -init 0 -a_nodata 0 -te $xmin $ymin $xmax $ymax -co COMPRESS=DEFLATE -co ZLEVEL=9 ${coastlines_shp} coastline_mask.tif
##gdal_rasterize -b 1 -burn 50 -l projected_shapefile -tr 50 50 -ot Float32 -init 0 -a_nodata 0 -te xmin ymin xmax ymax -co COMPRESS=DEFLATE -co ZLEVEL=9 projected_shapefile.shp coastline_distance.tif
#

#------------------------------------------------
# OPEN:
# Una vez que tenemos SURF.tif, solo hay que asignar 0 a todo lo que sean LAND. y sumar (1 - SURF)
# OPEN = 1 - SURF
# OPEN = 0 when (LAND == 1)


#------------------------------------------------
# CHLLOR
#

#
#------------------------------------------------
# DMS
#


#
#------------------------------------------------
# Merge them all on a NetCDF:





#start_date_s=$(date -d "$start_date" +%s)
#end_date_s=$(date -d "$end_date" +%s)     
#YYYY=$(date -d @$start_date_s +"%Y")   #año
# DDD=$(date -d @$start_date_s +"%j")   #dia juliano
#
#variables=("OPEN" "SURF" "CHLO" "DMS")
#var_list=$(printf "%-16s" ${variables[@]})
#
##Crear estructura del NetCDF de salida
#cat  > netcdf_emission_template.cdl <<EOF
#netcdf ocean_1 {
#dimensions:
#    TSTEP = 1; // 24;
#    DATE_TIME = 2 ;
#    COL = $nx ;
#    ROW = $ny ;
#    LAY = 1 ;
#    VAR= ${#variables[@]};
#
#variables:
#    int TFLAG(TSTEP, VAR, DATE_TIME) ;
#    float Band1(ROW, COL);
#EOF
#printf "     float %s(TSTEP, LAY, ROW, COL);\n" ${polluts[@]} >> netcdf_emission_template.cdl #   //Mis variables:
#cat >> netcdf_emission_template.cdl << EOF
#// global attributes:
# :IOAPI_VERSION = "ioapi-3.2: \$Id: init3" ; :EXEC_ID = "???????????????? " ; 
# :FTYPE = 1 ; :SDATE = ${YYYY}${DDD} ; :STIME = 000000 ; :WDATE = 2023001 ; :WTIME = 000000 ; :CDATE = 2023001 ; :CTIME = 000000;
# :TSTEP = 10000; :NTHIK = 1;    
# :NCOLS = ${nx}; :NROWS = ${ny}; :NLAYS = ${nz}; :NVARS = ${#variables[@]}; :GDTYP = 2;
# :P_ALP = -50.;:P_BET = -20. ;:P_GAM = -65.; 
# :XCENT = ${ref_lon}; :YCENT = ${ref_lat}; :XORIG = ${xorig}; :YORIG = ${yorig}; :XCELL = ${dx};  :YCELL = ${dy};  
# :VGTYP = -9999 ;:VGTOP = 5000.f ;:VGLVLS = 1.f, 0.9938147f;    
# :GDNAM = "${GridName}";:UPNAM = "OUTCM3IO" ;   :VAR-LIST = "${var_list}";      
# :FILEDESC = "Merged emissions output file from Mrggrid" ; :HISTORY = "" ;
#}
#EOF
##---------------------------------------
#mkdir finn_data		#carpeta de datos descargados de FINN.
#thisYear=$( date +"%Y")
#
#day=$start_date_s
#while [ $day -le $end_date_s ]
#do
#	YYYY=$(date -d @$day +"%Y")   #año
#	  MM=$(date -d @$day +"%m")   #mes
#	 DDD=$(date -d @$day +"%j")   #dia juliano
#	  HH=$(date -d @$day +"%H")   #hora
#	echo "Dia: $YYYY-$DDD (mes: $MM)"
#
#	#=================
#	#Transformo coordenadas + me quedo con los puntos dentro del dominio.
#	cat $finnFile | sed 's/D\([-+]\)/e\1/g;s/,/;/g' | awk -F";" 'NR==1{print $0";wkt"} NR>=2{print $0"POINT("$5*1.0,$4*1.0")"}' > tmp.csv
#	ogr2ogr -f CSV tmp_ll_clipped.csv tmp.csv -clipsrc $lonmin $latmin $lonmax $latmax -lco GEOMETRY="AS_WKT" -lco SEPARATOR="SEMICOLON"
#	ogr2ogr -f CSV tmp_xy.csv tmp_ll_clipped.csv -s_srs "$srsInp" -t_srs "$srsOut" -lco GEOMETRY="AS_WKT" -lco SEPARATOR="SEMICOLON"
#	sed -i 's/"//g;s/  */ /g;' tmp_xy.csv
#	#=================
#	#Un archivo por hora:
#	
#	for HH in $(seq --format="%02.0f" 0 23) 
#	do
#		echo "   HORA: $HH"
#		i="$(printf "%d" ${HH})";
#	
#		#Peso horario:
#		wgt_hh=${ciclo_diurno[i]}
#
#		#Creo NetCDF destino:
#		file_out=fire_emis_${YYYY}${DDD}_${HH}:00:00_d01.nc
#		ncgen -o $file_out netcdf_emission_template.cdl
#	
#		#---------------------------------------
#		for j in ${!polluts[@]}
#		do
#			pollut=${polluts[$j]}
#			echo "      pollut: $pollut"
#			#---------------------------------------
#			#Filtro columnas wkt y pollut, tambien convierto emis a mole/s ó g/s (aplicando ciclo diurno).
#			awk -F";" -v pollut=${pollut} -v wgt_hh=${wgt_hh} '
#			NR==1{
#			for(i=1;i<=NF;i++){if($i==pollut){pol=i};};
#				print "wkt;emis";
#				if ( pollut == "OC" || pollut =="BC" || pollut=="PM25" || pollut == "PM10" ){
#					k=1/3600000*wgt_hh;	#kg/day -> g/s	
#				}else {
#					k=1/3600.0 *wgt_hh;	#mole/day -> mole/s
#				};
#			} 
#			NR>1{ printf("%s;%.5e\n", $1, $pol*k); }' tmp_xy.csv > tmp_pollut.csv
#			#---------------------------------------
#			##ASCII to netcdf
#			#gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" -ot Float32 tmp_pollut.csv tmp.nc
#			gdal_rasterize -q -a emis -add -a_srs "${srsOut}" -tr $dx $dy -te $xmin $ymin $xmax $ymax -of GTiff tmp_pollut.csv tmp.tif
#			gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" tmp.tif tmp.nc #-co "WRITE_LONLAT=YES" -co "WRITE_BOTTOMUP=NO"
#	
#			ncks  -h -A -V -v Band1 tmp.nc -o $file_out
#			ncap2 -h -A -C -s "${pollut}(0,0,:,:) = Band1(:,:); TFLAG(0,$j,0) = ${YYYY}${DDD}; TFLAG(0,$j,1) =  ${HH}0000;" $file_out 
#
#		done;
#		
#	done;
#	day=$(($day + 86400 )) # Incrementar day un dia
#done;
#
