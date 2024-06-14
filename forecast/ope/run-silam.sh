#!/bin/bash
#===========
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#   Prepara directorio para corrida de Silam              #
#* * * * * * * * * * * * ** * * * * * * * * * * * * * * * #
#INPUTS:
dir=${HOME}/forecast            #forecast base-directory

SILAM_EXE=${HOME}/m/silam-model/bin/silam_v5_8pub.gnu

#data needed:
silam_dir=${dir}/dev/silam      #ruta a dir con namelists y tablas (ini/)
wrf_dir=${dir}/ope/wrf          #ruta a dir con wrfouts
emis_dir=emis/silam
#-------------------
day=`date +'%Y%m%d'`
exp_name=slm

echo -e " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo -e "  %Inicio prep silam                              "
echo -e ""

if [ ! -d slm ]; then mkdir ${exp_name}; fi
cd $exp_name

# data
ln -s $wrf_dir/wrfout* .

ln -s $silam_dir/ini .
ln -s $silam_dir/emis .
ln -s $silam_dir/static .

# namelists
cp $silam_dir/config_bcon_cb5.ini .
cp $silam_dir/config_out.ini .
cp $silam_dir/silam.ini .

# exe
ln -s ${SILAM_EXE} silam.exe

#
#######################
#Date-time
start_time=`date "+%Y %m %d 0 0 0.0"`
sed -i "s/start_time =.*/start_time=${start_time}/g" silam.ini

cd ..
find ${exp_name}/ -print
echo " * * * * * * * * * * * * * * * * * * * * * * * * *"
echo " %Fin. "
exit

##(3) RUN
#
##Archivos de control:
##
## + *CONTROL FILE*: es un namelist que contiene 8 secciones. cada seccion empieza con "LIST <sec. name>" y termina con "END LIST <sec. name>".
##   ├── general_parameters
##   │   └──internal model setup
##   │      ├── nuclide data file
##   │      ├── nuclide decay data
##   │      ├── land-use data
##   │      ├── optical properties
##   │      └── chemical properties
##   │
##   ├── emission_parameters
##   │   └── source term file/s
##   ├── dispersion_parameters
##   ├── meteo_parameters
##   ├── transphormation_parameters
##   ├── initial_and_boundary_conditions
##   │   └── bcon_config_file
##   ├── optical_density_parameters
##   └── output_parameters
##	├── output configuration file
##       └── output configuration file
#
##A su vez el archivo de control puede llamar a los siguientes archivos:
##   - *source term file:* describe fuentes de emisiones.
##   - *output configuration file:* describe la configuración de salida.
##   - *internal model setup:* define la configuración intena del modelo.
##
##El **internal model setup** a su vez puede llamar a dos archivos:
##   - *standard cocktails file:* defines the standard cocktails that can be used in the source description; referred from the internal setup file. Users are free to create their own cocktails, adding to the existent file;
##   - GRIB or NetCDF code table (depending of the type of files): mandatory, invisible for users, referred from the internal setup file.
##
###Archivos opcionales:
##Dependiendo de la configuración de la corrida hay varios archivos que pueden ser incluidos:
##   - nuclide data file: for radioactive simulations, invisible for users, referred from the internal setup file;
##   - nuclide decay data file: for radioactive simulations, invisible for users, referred from the	internal setup file;
##   - land-use data: for chemical simulations of biogenic emissions, invisible for users, referred from the internal setup file;
##   - optical properties: for chemical and aerosol simulation, describes the optical properties of substances, invisible for users, referred from the internal setup file;
##   - chemical properties: describes the chemical properties of the species available in SILAM, invisible for users, referred from the internal setup file;
#
#
##DUST:
##
## Archivos necesarios:
##	- dust_emis_0_v3.nc4
##	- BCON/ICON (Silam global model)
#
#
## BC & IC
## Descargar boundary conditions de Silam: OPeNDAP
##ncks -d lat,-54.0,-28.8 -d lon,-72.,-52. -v "a,a_half,b,b_half,dd_dust*" http://silam.fmi.fi/thredds/dodsC/silam_glob06_v5_8/runs/silam_glob06_v5_8_RUN_2023-07-02T00:00:00Z -o  bc.nc
#ncks -d lat,-80.0,-20.0 -d lon,-78.,-40. -v "a,a_half,b,b_half,vmr_dust*" http://silam.fmi.fi/thredds/dodsC/silam_glob06_v5_8/runs/silam_glob06_v5_8_RUN_2023-07-02T00:00:00Z -o  bc.nc
#
##Modificar fechas
#ncatted -O -a units,time,o,c,"hours since 2009-03-10 00:00:00.000 UTC" bc.nc bc_modif_date.nc
#ncatted -O -a SIMULATION_START_DATE,global,o,c,"2009-03-15T00:00:00"  bc_modif_date.nc bc_modif_date.nc
#
#
##Comandos utiles (que aprendi en el camino)..
##  Generar serie temporal de masa total:
#cdo fldsum -selname,cnc_dust_m1_5 -vertsum output/toypoint.nc ts_dust_m1_5.nc
## Sumar en dimension vertical (para ver masa en columna)
#cdo -selname,cnc_dust_m1_5 -vertsum output/toypoint.nc ts_dust_m1_5.nc
#
#
#
#
#
#
