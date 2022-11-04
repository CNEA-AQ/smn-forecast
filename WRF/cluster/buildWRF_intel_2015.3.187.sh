#!/bin/bash
#================================================
#   WRF: Build & Run.
#================================================
# Dependencias:
module purge
module load intel/intel_2015.3.187
module load ompi3.1.4_intel_2015.3.187
module load hdf5_1.8.17_intel_2015.3.187
module load miscLibs_intel_2015.3.187
module load netcdf4.4_intel_2015.3.187_parallel

#------------------------------------------------
# Build:

LIBSDIR=/home/ramiroespada/libs_intel_2015.3.187
OMPI_PATH=$LIBSDIR/ompi

#export PATH=$PATH:$OMPI_PATH/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMPI_PATH/lib
export MANPATH=$MANPATH:$OMPI_PATH/share/man
export OMPI_HOME=$OMPI_PATH
export MPI_HOME=$OMPI_PATH 

HDF5=$LIBSDIR/hdf5
NETCDF=$LIBSDIR/pnetcdf
AEC=$LIBSDIR/grib2 #/sw/rhel6-x64/sys/libaec-1.0.2-gcc48

export NETCDF=$NETCDF
export HDF5=$HDF5
export PHDF5=$HDF5

export FLEX_LIB_DIR=$LIBSDIR/grib2/lib
export YACC=$LIBSDIR/grib2/bin/yacc
export JASPERINC=$LIBSDIR/grib2/include
export JASPERLIB=$LIBSDIR/grib2/lib

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export NETCDF4=1
unset NETCDF_classic
export WRF_EM_CORE=1
export WRF_NMM_CORE=0
export WRF_CHEM=1
export WRF_KPP=0

echo 'Cleaning WRF ...'
./clean -a

echo 'Configure WRF ...'
/usr/bin/expect<<EOF
set timeout 18000
expect_after timeout { puts "TIMEOUT"; exit 1 }
spawn ./configure
match_max 100000
expect -re "Enter selection .* : "
send -- "15\r"
expect -re "Compile for nesting?.*: "
send -- "\r"
expect eof

EOF

# modify macros to use rpath while linking
sed -i "s@-lnetcdff@-L${NETCDF}/lib -L${AEC}/lib -Wl,-rpath,${NETCDF}/lib -Wl,-rpath,${NETCDF}/lib -Wl,-rpath,${HDF5}/lib -lnetcdff -lnetcdf@" configure.wrf

#Obtuve error por que no encontraba lhdf5_fortran, asi que lo saquéy compiló
# (todavía no sé si tirará err
#sed -i "s/lhdf5_fortran//" configure.wrf


echo 'Compile WRF ...'
./compile em_real &> log_compile


# Si la compilación es exitosa, deberian crearse:
# -  real.exe
# -  wrf.exe

#------------------------------------------------
# Run:
