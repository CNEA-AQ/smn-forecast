#!/bin/bash
#=============================================================================
#Cargar compilador (intel_2015.3.187)
module purge
module load intel/intel_2015.3.187
export compilerName="intel_2015.3.187"
export LIBSDIR=/home/ramiroespada/libs_${compilerName}
# -----------------------------------------------------------------------------
# En la carpeta $LIBSDIR voy a guardar todas las librerias.
mkdir $LIBSDIR
cd $LIBSDIR
# =============================================================================
# OpenMPI (openmpi-1.10.2)
CC=icc CXX=icpc FC=ifort F90=ifort F77=ifort ./configure --prefix=${LIBSDIR}/ompi
#(!) Se puede extender para agregar pmi, mallox, fca, etc.
make -j 16 all
make -j 16 install

export OMPI=$LIBSDIR/ompi
export PATH=$PATH:$OMPI/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMPI/lib

#==============================================================================
# NetCDF-C (netcdf-4.4.0)
#sin pNetCDF:
incs="-I${LIBSDIR}/grib2/include"
libs="-L${LIBSDIR}/grib2/lib"
export CPPFLAGS="${incs} ${libs}"
export CXXFLAGS="${incs} ${libs}"
export   CFLAGS="${incs} ${libs}"
export   FFLAGS="${incs} ${libs}"
export  FCFLAGS="${incs} ${libs}"
export  LDFLAGS="${incs} ${libs}"
CC=icc CXX=ccpc FC=ifort F90=ifort F77=ifort MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort ./configure --prefix=$LIBSDIR/netcdf --enable-fortran --enable-shared 
make
make install
#==============================================================================
# NetCDF-Fortran (netcdf-fortran-4.4.3)
#Usando mismos flags que en NetCDF-C!!
CC=icc CXX=cxx FC=ifort F90=ifort F77=ifort ./configure --prefix=$LIBSDIR/netcdf --enable-shared
make
make install

#==============================================================================
#IOAPI (ioapi-3.2.0): (!) Compilado pero sin m3tools
mkdir ioapi-3.2
tar -xvf ../TARs/ioapi-3.2.tar.gz ioapi-3.2
cd ioapi-3.2

cp Makefile.template Makefile
cp ioapi/Makefile.pncf ioapi/Makefile
cp m3tools/Makefile.pncf m3tools/Makefile


make configure BIN=Linux2_x86_64ifortmpi CPLMODE=pncf BASEDIR=/home/ramiroespada/libs_intel_2015.3.187/ioapi-3.2 NCFLIBS="-L$LIBSDIR/ompi/lib -I$LIBSDIR/ompi/include -L$LIBSDIR/hdf5/lib -I$LIBSDIR/hdf5/include -L$LIBSDIR/pnetcdf/lib -I$LIBSDIR/pnetcdf/include -L$LIBSDIR/grib2/lib -I$LIBSDIR/grib2/include -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz" IOAPIDEFS="-DIOAPI_PNCF -DIOAPI_NCF4"
make BIN=Linux2_x86_64ifortmpi CPLMODE=pncf BASEDIR=/home/ramiroespada/libs_intel_2015.3.187/ioapi-3.2 NCFLIBS="-L$LIBSDIR/ompi/lib -I$LIBSDIR/ompi/include -L$LIBSDIR/hdf5/lib -I$LIBSDIR/hdf5/include -L$LIBSDIR/pnetcdf/lib -I$LIBSDIR/pnetcdf/include -L$LIBSDIR/grib2/lib -I$LIBSDIR/grib2/include -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz" IOAPIDEFS="-DIOAPI_PNCF -DIOAPI_NCF4"

make
#Al finalizar el make, en la carpeta $BIN se va a crear libioapi.a
#(!) no compila bien el m3tools por que no encuentra libs de ompi que si existen =/

#==============================================================================

