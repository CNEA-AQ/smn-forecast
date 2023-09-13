#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

   srsInp="epsg:4326"	#Los archivos vienen en latlon.

#Input files:
 GRIDDESC="GRIDDESC"
 landpoly=input/shp/land_10m_global.gpkg	   #polygons to make a landmask
chlo_file=input/CHLO_monthly_global_2019.nc     #chlorophyl conc. files (from MODIS (MY1DMM)
 dms_file=input/DMS_monthly_global_2011.nc       #dms conc. files (DMS climatology of Lana et al. (2011))
#----------
#Parametros de Proyeccion y Grilla (GRIDDESC)
Gridname="test"
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell
xorig=-1970000;yorig=-2370000;

read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${srsOut}" <<< $( echo "${ref_lon} ${ref_lat}" ) )
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 
#------------------------------------------------
#Parametros intermedios:
xmin=$(bc -l <<<" $xc - ($xorig)*-1 ")
xmax=$(bc -l <<<" $xc + ($xorig)*-1 ")
ymin=$(bc -l <<<" $yc - ($yorig)*-1 ")
ymax=$(bc -l <<<" $yc + ($yorig)*-1 ")

#bbox latlon: para hacer clippear el archivo de fuegos al dominio:
read lonmin latmin ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmin} ${ymin}" ) )
read lonmax latmax ellipsoidh <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${srsOut}" <<< $( echo "${xmax} ${ymax}" ) )

# Crear NetCDF base:
variables=("OPEN" "SURF" "CHLO" "DMS")
var_list=$(printf "%-16s" ${variables[@]})

#Crear estructura del NetCDF de salida
cat  > netcdf_emission_template.cdl <<EOF
netcdf ocean_1 {
dimensions:
    TSTEP = 1; DATE_TIME = 2; COL = $nx; ROW = $ny; LAY = 1; VAR= ${#variables[@]};
variables:
    int TFLAG(TSTEP, VAR, DATE_TIME) ;
    float Band1(ROW, COL);
EOF
printf "     float %s(TSTEP, LAY, ROW, COL);\n" ${variables[@]} >> netcdf_emission_template.cdl #   //Mis variables:
cat >> netcdf_emission_template.cdl << EOF
// global attributes:
 :IOAPI_VERSION = "ioapi-3.2: \$Id: init3" ; :EXEC_ID = "???????????????? " ; 
 :FTYPE = 1 ; :SDATE = -635; :STIME = 000000 ; :WDATE = 2023001 ; :WTIME = 000000 ; :CDATE = 2023001 ; :CTIME = 000000;
 :TSTEP = 0; :NTHIK = 1;    
 :NCOLS = ${nx}; :NROWS = ${ny}; :NLAYS = ${nz}; :NVARS = ${#variables[@]}; :GDTYP = 2;
 :P_ALP = -50.;:P_BET = -20. ;:P_GAM = -65.; 
 :XCENT = ${ref_lon}; :YCENT = ${ref_lat}; :XORIG = ${xorig}; :YORIG = ${yorig}; :XCELL = ${dx};  :YCELL = ${dy};  
 :VGTYP = -9999 ;:VGTOP = 5000.f ;:VGLVLS = 1.f, 0.9938147f;    
 :GDNAM = "${GridName}";:UPNAM = "OUTCM3IO" ;   :VAR-LIST = "${var_list}";      
 :FILEDESC = "Merged emissions output file from Mrggrid" ; :HISTORY = "" ;
}
EOF
#------------------------------------------------
# SURF:
#(1) Rasterizo las lineas de costas "muy fino" (100m o 200m), asigno 1/0 si toca o no costa. Luego regrid a grilla de modelado calculando "average".
dx_small=200; dy_small=200
echo "Calculating SURF.."
if [ ! -f SURF.tif ] 
then
    ogr2ogr -f "GPKG" coastline.gpkg $landpoly -nln coast -nlt MULTILINESTRING -t_srs "${srsOut}" -clipsrc $lonmin $latmin $lonmax $latmax #-explodecollections 
    gdal_rasterize -q -l coast -at -burn 1 -a_srs "${srsOut}" -tr ${dx_small} ${dy_small} -te $xmin $ymin $xmax $ymax -of GTiff coastline.gpkg tmp.tif
    gdalwarp -r average -t_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax tmp.tif SURF.tif 
    rm tmp.tif
else
    echo "SURF.tif already exists!"
fi
# LAND: (same procedure than with SURF)
echo "Calculating LAND.."
if [ ! -f LAND.tif ] 
then
    ogr2ogr -f "GPKG" land.gpkg $landpoly -nln land -nlt MULTIPOLYGON -t_srs "${srsOut}" -clipsrc $lonmin $latmin $lonmax $latmax
    gdal_rasterize -q -l land -at -burn 1 -a_srs "${srsOut}" -tr ${dx_small} ${dy_small} -te $xmin $ymin $xmax $ymax -of GTiff land.gpkg tmp.tif
    gdalwarp -r average -t_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax tmp.tif LAND.tif 
    rm tmp.tif
else
   echo "LAND.tif already exists!"
fi
# OPEN = 1 - LAND - SURF
echo "Calculating OPEN.."
if [ ! -f OPEN.tif ] 
then
    gdal_calc.py -S SURF.tif --S_band=1 -L LAND.tif --L_band=1 --outfile=OPEN.tif --calc="1-L-S" 
else
    echo "OPEN.tif already exists!"
fi
#--------------------------------
#Merge them all in a NetCDF file:
for mm in $(seq --format="%02.0f" 0 12) 
do
	#Creo NetCDF destino:
	file_out=fire_emis_${YYYY}${DDD}_${HH}:00:00_d01.nc
	ncgen -o $file_out netcdf_emission_template.cdl

    gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" ${file} tmp0.nc 
    ncks -4 -A -v Band1 tmp0.nc tmp.nc
    ncrename -v Band1,chlo${mm} tmp.nc
    #ncatted -h -O -a esri_pe_string,laiv${mm},d,c,"" tmp.nc
    ncatted -h -O -a units,chlo${mm},o,c,"mg m^-3" tmp.nc
    ncatted -h -O -a long_name,chlo${mm},o,c,"CHLO" tmp.nc
    ncatted -h -O -a standard_name,chlo${mm},o,c,"mass_concentration_of_chlorophyll_in_sea_water" tmp.nc
    ncatted -h -O -a var_desc,chlo${mm},o,c,"Chlorophyll Concentration, from MY1DMM MODIS product." tmp.nc
    ncatted -h -O -a missing_value,chlo${mm},o,s,"99999" tmp.nc
done
rm tmp0.nc
mv tmp.nc OCEAN_MASK.nc
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
