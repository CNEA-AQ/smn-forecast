#!/bin/bash

#Dependencias de CMAQ:
#   * C y Fortran Compilers  (GNU> 6.1 | Intel > 17.0)
module purge
module load gcc/6.3.0
#   * MPI Library (Minimum versions: IntelMPI 2017.0 | MPICH 3.3.1 | MVAPICH2 2.3.1 | OpenMPI 2.1.0)
module load mpich3.1.4_gcc_6.3.0
#   * NetCDF Library: NetCDF-C y NetCDF-Fortran ( (!) sin HDF4, HDF5,DAP, PnetCDF ni Zlib ) (Minimum versions: NetCDF-C 4.2 | NetCDF-Fortran 4.4.2)
module load netcdf4.4_gcc_6.3.0
#   * I/O API Library
#=============================================================================
export compilerName="gcc_6.3.0"
export LIBSDIR=/home/ramiroespada/libs_${compilerName}
# -----------------------------------------------------------------------------
# En la carpeta $LIBSDIR voy a guardar todas las librerias.
#mkdir $LIBSDIR
#cd $LIBSDIR
#==============================================================================
#MPICH (mpich-3.1.4)
CC=gcc CXX=g++ FC=gfortran ./configure --prefix=${LIBSDIR}/mpich
make
make install
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBSDIR}/mpich/lib
export PATH=${PATH}:${LIBSDIR}/mpich/bin
#==============================================================================
# NetCDF-C (netcdf-4.4.0)
#sin pNetCDF, HDF5, HDF4, DAP:
incs="-I${LIBSDIR}/netcdf/include" 
libs="-L${LIBSDIR}/netcdf/lib"     
export CPPFLAGS="${incs} ${libs}"
export CXXFLAGS="${incs} ${libs}"
export   CFLAGS="${incs} ${libs}"
export   FFLAGS="${incs} ${libs}"
export  FCFLAGS="${incs} ${libs}"
export  LDFLAGS="${incs} ${libs}"
CC=gcc CXX=g++ FC=gfortran F90=gfortran F77=gfortran ./configure --prefix=$LIBSDIR/netcdf --enable-fortran --enable-shared --disable-netcdf-4 --disable-dap
make
make install
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBSDIR}/netcdf/lib
export PATH=${PATH}:${LIBSDIR}/netcdf/bin
#------------------------------------------------------------------------------
# NetCDF-Fortran (netcdf-fortran-4.4.3)
#Usando mismos flags que en NetCDF-C!!
CC=gcc CXX=g++ FC=gfortran F90=gfortran F77=gfortran ./configure --prefix=$LIBSDIR/netcdf --enable-shared
make
make install
#==============================================================================
#IOAPI (ioapi-3.2.0): (!) Compilado pero sin m3tools
git clone https://github.com/cjcoats/ioapi-3.2
cd ioapi-3.2

cp Makefile.template Makefile
cp ioapi/Makefile.nocpl ioapi/Makefile
cp m3tools/Makefile.nocpl m3tools/Makefile

#Agregar al Makefile:
```
NCFLIBS = -L/home/ramiroespada/libs_gcc-6.3.0/netcdf/lib -lnetcdf -lnetcdff
```

# Comentar cualquier flag de openMPI en los "Makeinclude.Linux2_x86_64gfort" si es que los hay.
```
OMPFLAGS = # -fopenmp
OMPLIBS = # -fopenmp
```

make configure BIN=Linux2_x86_64gfort CPLMODE=nocpl BASEDIR=/home/ramiroespada/libs_gcc_6.3.0/ioapi-3.2 
make BIN=Linux2_x86_64gfort CPLMODE=nocpl BASEDIR=/home/ramiroespada/libs_gcc_6.3.0/ioapi-3.2 

#Al finalizar el make, en la carpeta $BIN se va a crear libioapi.a
#Tuve problemas con gfortran-11 (gfortran-9 si funciona)
#==============================================================================
# SMOKE
# Dependencias: I/O API, NetCDF (libnetcdf y libnetcdff)
git clone https://github.com/CEMPD/SMOKE
cd SMOKE
SMK_HOME=`pwd`
ln -s ~/libs_gcc_6.3.0/ioapi-3.2 ioapi    #link simbolco a donde instale ioapi
mkdir Linux2_x86_64gfort                 #carpeta donde van a ir los *.o *.mod y exes
cd src
#Editar "Makeinclude" ~~~
#>    BASEDIR=${SMK_HOME}/src
#>    INCDIR=${BASEDIR}/inc
#>    OBJDIR=${SMK_HOME}/${BIN}
#>    
#>    IOBASE=${SMK_HOME}/ioapi
#>    IODIR=${IOBASE}/ioapi
#>    IOBIN=${IOBASE}/${BIN}
#>    IOINC=${IODIR}/fixed_src
#>    
#>    INSTDIR=${OBJDIR}
#>    #(!) No olvidar descomentar flags compilador a usar.
#~~~~~~~~~~~~~~~~~~~~~~~~
#make, make install:
make FC=gfortran BIN=Linux2_x86_64gfort SMK_HOME=${SMK_HOME}
make install FC=gfortran BIN=Linux2_x86_64gfort SMK_HOME=${SMK_HOME}
#==============================================================================
