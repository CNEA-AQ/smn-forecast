#!/bin/bash
#=========================================================#
#   Prepara directorio para corrida de WRF                #
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#INPUTS:
dir=/home/ramiroespada
#src:
WPSSRC=${dir}/WPS #stage/WRFV4gnu/WPS   #Ruta al source del WPS
WRFSRC=${dir}/WRF #stage/WRFV4gnu/WRF   #Ruta al source del WRF
#data:
meteo_path=${dir}/stage/data/GFS2 #ruta a archivos meteorologicos
geog_path=${dir}/stage/data/WPS_GEOG
#default namelists and tables:
namelist_path=${dir}/WRF/misc/namelists
tables_path=${dir}/WRF/misc/tablas
sbatch_path=${dir}/WRF/misc/sbatch

#######################
##Bahía  16 marzo 2009 
exp_name=Bahia16mar2009;cdate="2009-03-16 00:00:00";clat=-39.704213;clon=-61.833573;DX=600;DY=600

#Patagonia 28 marzo 2009
#exp_name=Patagonia28mar2009;cdate="2009-03-28 00:00:00";clat=-43.004952;clon=-64.153677;DX=800;DY=1000

##Madryn 05 de abril 2009
#exp_name=Madryn05abr2009;cdate="2009-04-05 00:00:00";clat=-43.004952;clon=-64.153677;DX=800;DY=1000
#
#Marchiquita + Bahía  05 de agosto 2009.
#exp_name=MarChiquita05ago2009;cdate="2009-08-05 00:00:00";clat=-30.562957;clon=-62.506265;DX=600;DY=600
#
##Cuyo  11-12 julio 2010
#exp_name=Cuyo11jul2010;cdate="2010-07-11 00:00:00";clat=-32.951274;clon=-68.854580;DX=500;DY=700
#
##Munster 03 noviembre 2016
#exp_name=Munster03nov2016;cdate="2016-11-03 00:00:00";clat=-43.004952;clon=-64.153677;DX=850;DY=1000
#
##Antofagasta 13 septiembre 2017
#exp_name=Antofagasta13sep2017;cdate="2017-09-13 00:00:00";clat=-22.896295;clon=-67.589588;DX=600;DY=800
######################

dx=20000;
dy=20000;
nx=$(bc -l <<< "scale=5;o=$DX*1000/$dx;scale=0;o/1.0") # == nx
ny=$(bc -l <<< "scale=5;o=$DY*1000/$dy;scale=0;o/1.0") # == ny

start_date=$(TZ=0 date -d "${cdate}+0 - 6 hours" '+%Y-%m-%d_%H:%M:%S')  #ojo, chequiar dispo de meteo-data
  end_date=$(TZ=0 date -d "${cdate}+0 + 1 days " '+%Y-%m-%d_%H:%M:%S')

add2namelist="
e_we  = $nx,
e_sn  = $ny,  
ref_lat   = ${clat},
ref_lon   = ${clon},
dx        = ${dx},
dy        = ${dy},
time_step = 120,
geog_data_res = 'default',
map_proj  = 'lambert',
truelat1  = -33.721,
truelat2  = -33.721,
stand_lon = -59.834,
interval_seconds = 21600,
num_metgrid_levels = 27,
ra_lw_physics=1, ra_sw_physics = 2,radt=30,
cu_physics=5,mp_physics=28,bl_pbl_physics=1,
sf_sfclay_physics=1,sf_surface_physics=2,sf_urban_physics=0,
chem_in_opt=0,kemit=1,emiss_opt=0,emiss_inpt_opt=0,
chem_opt=401,dust_opt=1,dust_schme=1,
aer_drydep_opt=1,aerchem_onoff=1,
aer_op_opt=1,
opt_pars_out=0,
"
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#+========================================================#
#
echo -e " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo -e "  %Inicio.                                        "
echo -e ""
#check directories existence
echo -e "\e[36m Checkeando existencia de directorios:            \e[0m"    
DIRECTORIES=($WPSSRC $WRFSRC $meteo_path $geog_path $namelist_path $tables_path)
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
        mkdir ${exp_name};
        ls | grep ${exp_name};
fi

cd ${exp_name}

#echo -e "\e[36m -------------------------------------------------\e[0m"
#echo -e "\e[36m  Creando dir 'wps'                               \e[0m"
#mkdir wps
#cd wps

#namelist
cp ${namelist_path}/namelist.wps namelist.wps
        for var_inp in "${add2namelist_NO_SPACES_array[@]}"; do
                var_name=$(echo $var_inp | sed 's/\(.*\)=.*/\1/g')
        if grep -q $var_name namelist.wps; then
                echo -e "Agregando \e[32m ${var_inp}\e[0m a namelist.wps "
                sed -i "s/.*\b$var_name\b.*/$var_inp,/g" namelist.wps
        fi
        done
#tables
cp ${tables_path}/wps/GEOGRID.TBL GEOGRID.TBL
cp ${tables_path}/wps/Vtable  Vtable
cp ${tables_path}/wps/METGRID.TBL METGRID.TBL
#data
ln -s ${geog_path} .
cp -p ${WPSSRC}/link_grib.csh link_grib.csh
tcsh link_grib.csh ${meteo_path}/fnl*
#exe
ln -s ${WPSSRC}/geogrid.exe geogrid.exe
ln -s ${WPSSRC}/ungrib.exe ungrib.exe
ln -s ${WPSSRC}/metgrid.exe metgrid.exe

#cd ..
#find wps/ -print

#echo -e "\e[36m --------------------------------------------------\e[0m"
#echo -e "\e[36m  Creando dir 'wrf'\e[0m"
#mkdir wrf
#cd wrf

#namelist 
cp ${namelist_path}/namelist.input namelist.input
        for var_inp in "${add2namelist_NO_SPACES_array[@]}"; do
                var_name=$(echo $var_inp | sed 's/\(.*\)=.*/\1/g')
        if grep -q $var_name namelist.input; then
                echo -e "Agregando \e[32m${var_inp}\e[0m a namelist.input"
                sed -i "s/.*\b$var_name\b.*/\t$var_inp,/g" namelist.input
        fi

        done
#tables
cp -p ${tables_path}/wrf/* .
#exes
ln -s ${WRFSRC}/main/real.exe real.exe
ln -s ${WRFSRC}/main/wrf.exe wrf.exe
#sbatchs
cp ${sbatch_path}/sbatch.real .
cp ${sbatch_path}/sbatch.wrf .

#cd ..
#find wrf/ -print
cd ..
find ${exp_name}/ -print
echo " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo " %Fin. "
exit

