#!/bin/bash
#================================================
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#   Prepara directorio para corrida de Silam              #
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#INPUTS:
dir=${HOME}/forecast
#src:
SILAM_SRC=${HOME}/m/silam-model
#data:
met_dir=${dir}/wrf${date}            #ruta a archivos meteorologicos
ini_dir=${HOME}/data/silam/ini_v5_9   #ruta a archivos de geog
emi_dir=${HOME}/data/silam/emis      #ruta a archivos de geog

#######################
day=`date +'%Y%m%d'`
exp_name=wrf$day
cdate=`date +'%Y-%m-%d 00:00:00'` #current date #"2017-09-13 00:00:00"
clat=-35.0 # domain center point (lat)   #-36.0
clon=-65.0 # domain center point (lon)   #-67.0
DX=2800    # domain -width  [km]
DY=4600    # domain -length [km]
#---------------------

dx=20000 #cell size-X [m] 
dy=20000 #cell size-Y [m] 

nx=$(bc -l <<< "scale=5;o=$DX*1000/$dx;scale=0;o/1.0" ) # == nx
ny=$(bc -l <<< "scale=5;o=$DY*1000/$dy;scale=0;o/1.0" ) # == ny

start_date=$(TZ=0 date -d "${cdate}+0 - 0 hours" '+%Y-%m-%d_%H:%M:%S')  #ojo, chequiar dispo de meteo-data
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
truelat1=-20.0,
truelat2=-50.0,
stand_lon=-65.0,
interval_seconds = 21600,
num_metgrid_levels = 34,
ra_lw_physics=1, ra_sw_physics = 2,radt=30,
cu_physics=5,mp_physics=8,bl_pbl_physics=1,
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
DIRECTORIES=($WPSSRC $WRFSRC $met_path $geo_path $namelist_path $tables_path)
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
ln -s ${geo_path} .
cp -p ${WPSSRC}/link_grib.csh link_grib.csh
#tcsh link_grib.csh ${met_path}/fnl*
tcsh link_grib.csh ${met_path}/gfs.*
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

#(3) RUN

#Archivos de control:
#
# + *CONTROL FILE*: es un namelist que contiene 8 secciones. cada seccion empieza con "LIST <sec. name>" y termina con "END LIST <sec. name>".
#   ├── general_parameters
#   │   └──internal model setup
#   │      ├── nuclide data file
#   │      ├── nuclide decay data
#   │      ├── land-use data
#   │      ├── optical properties
#   │      └── chemical properties
#   │
#   ├── emission_parameters
#   │   └── source term file/s
#   ├── dispersion_parameters
#   ├── meteo_parameters
#   ├── transphormation_parameters
#   ├── initial_and_boundary_conditions
#   │   └── bcon_config_file
#   ├── optical_density_parameters
#   └── output_parameters
#	├── output configuration file
#       └── output configuration file

#A su vez el archivo de control puede llamar a los siguientes archivos:
#   - *source term file:* describe fuentes de emisiones.
#   - *output configuration file:* describe la configuración de salida.
#   - *internal model setup:* define la configuración intena del modelo.
#
#El **internal model setup** a su vez puede llamar a dos archivos:
#   - *standard cocktails file:* defines the standard cocktails that can be used in the source description; referred from the internal setup file. Users are free to create their own cocktails, adding to the existent file;
#   - GRIB or NetCDF code table (depending of the type of files): mandatory, invisible for users, referred from the internal setup file.
#
##Archivos opcionales:
#Dependiendo de la configuración de la corrida hay varios archivos que pueden ser incluidos:
#   - nuclide data file: for radioactive simulations, invisible for users, referred from the internal setup file;
#   - nuclide decay data file: for radioactive simulations, invisible for users, referred from the	internal setup file;
#   - land-use data: for chemical simulations of biogenic emissions, invisible for users, referred from the internal setup file;
#   - optical properties: for chemical and aerosol simulation, describes the optical properties of substances, invisible for users, referred from the internal setup file;
#   - chemical properties: describes the chemical properties of the species available in SILAM, invisible for users, referred from the internal setup file;


#DUST:
#
# Archivos necesarios:
#	- dust_emis_0_v3.nc4
#	- BCON/ICON (Silam global model)


# BC & IC
# Descargar boundary conditions de Silam: OPeNDAP
#ncks -d lat,-54.0,-28.8 -d lon,-72.,-52. -v "a,a_half,b,b_half,dd_dust*" http://silam.fmi.fi/thredds/dodsC/silam_glob06_v5_8/runs/silam_glob06_v5_8_RUN_2023-07-02T00:00:00Z -o  bc.nc
ncks -d lat,-80.0,-20.0 -d lon,-78.,-40. -v "a,a_half,b,b_half,vmr_dust*" http://silam.fmi.fi/thredds/dodsC/silam_glob06_v5_8/runs/silam_glob06_v5_8_RUN_2023-07-02T00:00:00Z -o  bc.nc

#Modificar fechas
ncatted -O -a units,time,o,c,"hours since 2009-03-10 00:00:00.000 UTC" bc.nc bc_modif_date.nc
ncatted -O -a SIMULATION_START_DATE,global,o,c,"2009-03-15T00:00:00"  bc_modif_date.nc bc_modif_date.nc


#Comandos utiles (que aprendi en el camino)..
#  Generar serie temporal de masa total:
cdo fldsum -selname,cnc_dust_m1_5 -vertsum output/toypoint.nc ts_dust_m1_5.nc
# Sumar en dimension vertical (para ver masa en columna)
cdo -selname,cnc_dust_m1_5 -vertsum output/toypoint.nc ts_dust_m1_5.nc






