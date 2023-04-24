#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

   srsInp="epsg:4326"	#Cambiar solo si los archivos NO vienen en latlon.

GRIDDESC_file="cat ~/runs/Argen20oct2022/cmaq/mcip/GRIDDESC"
GRIDNAME="Argentina"

#-------------------------------------------------
#Get grid & proj parameters from GRIDDESC:
read projName xorig yorig dx dy nx ny nz <<< $( sed -n "/${GRIDNAME}/{n;p;q}" $GRIDDESC_file )
read COORDTYPE P_ALP P_BET P_GAM XCENT YCENT <<< $( sed -n "/${projName}/{n;p;q}" $GRIDDESC_file )


#Parametros de Grilla y Proyección: 
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 
read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${srsOut}" <<< $( echo "${ref_lon} ${ref_lat}" ) )

#Grilla (sale de GRIDDESC)
nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell
xorig=-1970000;yorig=-2370000;


#------------------------------------------------
#Re-gridding:

gdalwarp -s_srs "epsg:4326" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear $inFile test.nc

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
