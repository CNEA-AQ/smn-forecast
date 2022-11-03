#!/bin/bash
#================================================
#  WRF: Build & Run. (En ubuntu)
#================================================
# Dependencias:
#  - Compiladores (C & Fortran)
#  - Open-MPI > 3.1.4
#  - HDF5 > 1.8.17
#  - GRIB2 (Eccodes)
#  - LIBAEC (Adaptive Entropy Coding library)
#  - NetCDF > 4.4 (C & Fortran)
#------------------------------------------------
# Instalar dependencias:
sudo apt install gfortran-9 libopenmpi-dev libhdf5-dev libeccodes-dev libnetcdf-dev libnetcdff-dev libaec-dev
#------------------------------------------------

git clone https://github.con/wrf_model/WRF
cd WRF
#------------------------------------------------

# Build:
export NETCDF=/usr	#Necesita que esté setiado despues rectificamos en el configure.wrf
export HDF5=/usr        #Necesita que esté setiado despues rectificamos en el configure.wrf

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export NETCDF4=1
#unset NETCDF_classic	#para mayor compresión de netCDFs
export WRF_EM_CORE=1
export WRF_NMM_CORE=0
export WRF_CHEM=1
export WRF_KPP=0


./clean -a

./configure	#elegir la opción 34 para gnu (dmp)

# Modificamos los links de las librerias de netcdf y hdf5 en el archivo de configure.wrf
#LIB_EXTERNAL=-L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -lz

#Obtuve error por que no encontraba lhdf5_fortran, asi que lo saquéy compiló
# (todavía no sé si tirará err
#sed -i "s/lhdf5_fortran//" configure.wrf

echo 'Compile WRF ...'
./compile em_real &> log_compile
#(!) TIP: Si durante la compilación no encuentra alguna libreria, verificar si está y donde con: ldconfig -p | grep "libreria"

# Si la compilación es exitosa, deberian crearse (buscar en la carpeta run/):
# -  real.exe
# -  wrf.exe

#------------------------------------------------
#Build WPS: (Preprocesador del wrf)

#Bajo el repo:
git clone https://github.con/wrf_model/WPS
cd WPS

./clean -a

./configure #elegir opción 3 para gnu + dmp

./compile &> log_compile
#Modificar paths y links en configure.wps
#WRF_LIB         =       -L$(WRF_DIR)/external/io_grib1 -lio_grib1 \
#                        -L$(WRF_DIR)/external/io_grib_share -lio_grib_share \
#                        -L$(WRF_DIR)/external/io_int -lwrfio_int \
#                        -L$(WRF_DIR)/external/io_netcdf -lwrfio_nf \
#                        -L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -lz
#COMPRESSION_LIBS    = -L/usr/lib/x86_64-linux-gnu -lpng -lz -L/home/usuario/m/libs/jasper/lib -ljasper
#COMPRESSION_INC     = -I/home/usuario/m/libs/jasper/include -I/usr/include

#(!) Es posible que la libreria JasPer no esté instalada. En tal caso descargar la version 1.900.1 de https://www.ece.uvic.ca/~frodo/jasper/#doc y compilar.

#------------------------------------------------
# Run:




