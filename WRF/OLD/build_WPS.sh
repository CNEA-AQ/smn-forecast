#!/bin/bash

module purge
module load intel

DIR=/work/mh0735/m300734/build/lib

OMPI_PATH=$DIR/ompi
export PATH=$PATH:$OMPI_PATH/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMPI_PATH/lib
export MANPATH=$MANPATH:$OMPI_PATH/share/man
export OMPI_HOME=$OMPI_PATH
export MPI_HOME=$OMPI_PATH

HDF5=$DIR/hdf5
NETCDF=$DIR/netcdf
AEC=$DIR/grib2 

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
