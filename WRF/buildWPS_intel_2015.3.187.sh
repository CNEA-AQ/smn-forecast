#!/bin/bash
#================================================
# WPS: Build & Run
#================================================
# Dependencias:
#   * WRF (REQUIERE COMPILAR PRIMERO EL WRF!)
module purge
module load intel/intel_2015.3.187
module load ompi3.1.4_intel_2015.3.187
module load hdf5_1.8.17_intel_2015.3.187 
module load miscLibs_intel_2015.3.187 
module load netcdf4.4_intel_2015.3.187_parallel 
#------------------------------------------------
# Build:
LIBSDIR=/home/ramiroespada/libs_intel_2015.3.187

export OMPI_PATH=$LIBSDIR/ompi
export MANPATH=$MANPATH:$OMPI_PATH/share/man
export OMPI_HOME=$OMPI_PATH
export MPI_HOME=$OMPI_PATH

export HDF5=$LIBSDIR/hdf5
export NETCDF=$LIBSDIR/pnetcdf
export AEC=$LIBSDIR/grib2 

export NETCDF=$NETCDF

echo 'Cleaning WPS ...'
./clean -a

echo 'Configure WPS ...'
/usr/bin/expect<<EOF
set timeout 18000
expect_after timeout { puts "TIMEOUT"; exit 1 }
spawn ./configure
match_max 100000
expect -re "Enter selection .* : "
send -- "17\r"
expect eof

EOF

# modify macros to use rpath while linking
sed -i "s@-lnetcdff@-L${NETCDF}/lib -L${HDF5}/lib -L${AEC}/lib -Wl,-rpath,${NETCDF}/lib -Wl,-rpath,${NETCDF}/lib -Wl,-rpath,${HDF5}/lib -lnetcdff@" configure.wps

echo 'Compile WPS ...'
./compile &> log_compile



# Si la compilación es exitosa. Se van a crear:
#
#  - geogrid.exe  : genera grilla con topografia, surface, etc.
#  - ungrib.exe   : extrae datos de modelo global (fnl) para ser usado por metgrid.
#  - metgrid.exe  : construye con lo anterior los met_em.d01..nc que va a leer el wrf.

#------------------------------------------------
# Run:

# ejecutar ./buildRun.sh  #prepara directorio con todos los datos

# (1) ./geogrid.exe
# (2) ./ungrib.exe
# (3) ./metgrid.exe




