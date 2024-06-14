#!/bin/bash
#=========================================================#
#   Prepara directorio para corrida de WRF                #
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
day=`date +'%Y-%m-%d'`
basedir=${HOME}/forecast

#src:
WPSSRC=${HOME}/m/WPS   #Ruta al source del WPS
WRFSRC=${HOME}/m/WRF   #Ruta al source del WRF

#data & config files:
wrf_dir=${basedir}/dev/wrf              #ruta a namelists y tablas
met_path=${basedir}/ope/GFS             #ruta a archivos meteorologicos
geo_path=${HOME}/data/wrf/WPS_GEOG      #ruta a archivos de estáticos

#-----------------------------------------------------------
#day=`date +'%Y%m%d'`
exp_name=wrf #$day
cdate=$day"00:00:00" #`date +'%Y-%m-%d 00:00:00'` #current date #"2024-06-07 00:00:00" #
clat=-35.0 # domain center point (lat)   #-36.0
clon=-65.0 # domain center point (lon)   #-67.0
DX=2800    # domain -width  [km]
DY=4800    # domain -length [km]
#---------------------
#Parse grid params:
dx=16000 #20000  #20000 # 18000 #  #cell size-X [m]
dy=16000 #20000  #20000 # 18000 #  #cell size-Y [m]

nx=$(bc -l <<< "scale=5;o=$DX*1000/$dx;scale=0;o/1.0") # == nx
ny=$(bc -l <<< "scale=5;o=$DY*1000/$dy;scale=0;o/1.0") # == ny
#
start_date=$(TZ=0 date -d "${cdate}+0 - 0 hours" '+%Y-%m-%d_%H:%M:%S')  #ojo, chequiar dispo de meteo-data
  end_date=$(TZ=0 date -d "${cdate}+0 + 1 days " '+%Y-%m-%d_%H:%M:%S')

add2namelist="
e_we = $nx,
e_sn = $ny,  
ref_lat = ${clat},
ref_lon = ${clon},
dx      = ${dx},
dy      = ${dy},
time_step = 90,
geog_data_res = 'default',
map_proj  = 'lambert',           
truelat1=-20.0,
truelat2=-50.0,
stand_lon=-65.0,
interval_seconds = 21600,
num_metgrid_levels = 34,
"
#+========================================================#
#
echo -e " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo -e "  %Init. prep WRF run.                            "
echo -e ""
#check directories existence
echo -e "\e[36m Checkeando existencia de directorios:            \e[0m"    
DIRECTORIES=($WPSSRC $WRFSRC $met_path $geo_path)
for DIRECTORY in "${DIRECTORIES[@]}"; do
    if [ ! -d "$DIRECTORY" ]; then
       echo -e "\e[31m ERROR:  El directorio ${DIRECTORY} NO existe!\e[0m"; exit;
    else
       echo -e "Directorio ${DIRECTORY} \e[32m Ok, Existe!\e[0m"; continue;
    fi
done
#--------------------------------------------------------#
#parse dates:
echo -e "\e[36m -------------------------------------------------\e[0m"
echo -e "\e[36m Parseando fechas:                               \e[0m"    
     read start_year start_month start_day start_hour start_min start_sec <<< ${start_date//[-:\/_ ]/ }
     read end_year end_month end_day end_hour end_min end_sec <<< ${end_date//[-:\/_ ]/ }
echo "Fecha inicial:  "$start_year $start_month $start_day $start_hour
echo "Fecha final:  "$end_year $end_month $end_day $end_hour

diff_secs=$(($(date -d "$end_year-$end_month-$end_day $end_hour" '+%s')-$(date -d "$start_year-$start_month-$start_day $start_hour" '+%s')))
echo "Difrencia de tiempo (en seg.):  "$diff_secs
run_days=$(($diff_secs/(24*60*60)))
echo "Dias de corrida:  "$run_days
run_hours=$((($diff_secs-$run_days*24*60*60)/(60*60)))
echo "Diferencia horaria:  "$run_hours

add2namelist=$add2namelist",start_date=$start_date,end_date=$end_date,start_year=$start_year,start_month=$start_month,start_day=$start_day,start_hour=$start_hour,end_year=$end_year,end_month=$end_month,end_day=$end_day,end_hour=$end_hour,run_days=$run_days,run_hours=$run_hours"

printf '%s\n' "${my_array[@]}"

#parse variables on namelist:
echo -e "\e[36m -------------------------------------------------"
echo -e "\e[36m  Parseando variables de namelist:                "    
add2namelist_NO_SPACES=$(echo $add2namelist | sed 's/[[:space:]]//g')
add2namelist_NO_SPACES_array=(${add2namelist_NO_SPACES//[,;]/ })


echo -e "\e[36m -------------------------------------------------\e[0m"
echo -e "\e[36m  Creando directorio de trabajo '$exp_name'       \e[0m"    
if [ -d "$exp_name" ]; then
  echo -e "\e[31m ERROR:  El directorio $exp_name ya existe!\e[0m"; exit;
else
  mkdir ${exp_name}; ls | grep ${exp_name};
fi

cd ${exp_name}

#------
# WPS:

#namelist
cp ${wrf_dir}/namelist.wps namelist.wps
   for var_inp in "${add2namelist_NO_SPACES_array[@]}"; do
      var_name=$(echo $var_inp | sed 's/\(.*\)=.*/\1/g')
   if grep -q $var_name namelist.wps; then
      echo -e "Agregando \e[32m ${var_inp}\e[0m a namelist.wps "
      sed -i "s/.*\b$var_name\b.*/$var_inp,/g" namelist.wps
   fi
   done
#tables
cp ${wrf_dir}/GEOGRID.TBL GEOGRID.TBL
cp ${wrf_dir}/Vtable Vtable
cp ${wrf_dir}/METGRID.TBL METGRID.TBL
#exe
ln -s ${WPSSRC}/geogrid.exe geogrid.exe
ln -s ${WPSSRC}/ungrib.exe ungrib.exe
ln -s ${WPSSRC}/metgrid.exe metgrid.exe
#data
ln -s ${geo_path} .
cp -p ${WPSSRC}/link_grib.csh link_grib.csh
tcsh link_grib.csh ${met_path}/gfs.*

#------
# WRF:
#namelist 
cp ${wrf_dir}/namelist.input namelist.input
   for var_inp in "${add2namelist_NO_SPACES_array[@]}"; do
     var_name=$(echo $var_inp | sed 's/\(.*\)=.*/\1/g')
     if grep -q $var_name namelist.input; then
        echo -e "Agregando \e[32m${var_inp}\e[0m a namelist.input"
        sed -i "s/.*\b$var_name\b.*/\t$var_inp,/g" namelist.input
     fi
   done
#tables
cp -p ${wrf_dir}/wrf_tables/* .
#exes
ln -s ${WRFSRC}/main/real.exe real.exe
ln -s ${WRFSRC}/main/wrf.exe wrf.exe


## RUNING PROGRAMS:
#echo -e "\e[36m -------------------------------------------------\e[0m"
#echo -e "\e[36m  Ejecutando WPS                                  \e[0m"    
#
#echo -e "\e[36m  GEOGRID ...                                     \e[0m"    
#./geogrid.exe
#
#echo -e "\e[36m  UNGRIB ...                                      \e[0m"    
#./ungrib.exe
#
#echo -e "\e[36m  METGRID ...                                     \e[0m"    
#./metgrid.exe
#
#echo -e "\e[36m -------------------------------------------------\e[0m"
#echo -e "\e[36m  Ejecutando WRF                                  \e[0m"    
#
#echo -e "\e[36m  REAL.EXE                                        \e[0m"    
#mpirun real.exe
#
#echo -e "\e[36m  WRF.EXE                                         \e[0m"    
#mpirun wrf.exe



cd ..
find ${exp_name}/ -print
echo " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo " %Fin. "
exit

