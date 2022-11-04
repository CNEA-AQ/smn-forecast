#!/bin/bash
# 
# > Instalar librerias para CMAQ en ubuntu.
 
#Dependencias de CMAQ:
#   * C y Fortran Compilers  (GNU> 6.1 | Intel > 17.0)
#   * MPI Library (Minimum versions: IntelMPI 2017.0 | MPICH 3.3.1 | MVAPICH2 2.3.1 | OpenMPI 2.1.0)
#   * NetCDF Library: NetCDF-C y NetCDF-Fortran ( (!) sin HDF4, HDF5,DAP, PnetCDF ni Zlib ) (Minimum versions: NetCDF-C 4.2 | NetCDF-Fortran 4.4.2)
#   * I/O API Library
#==============================================================================
# MPICH (mpich-3.1.4)
sudo apt install mpich
#==============================================================================
# NetCDF-C (netcdf-4.4.0)
sudo apt install libnetcdf-dev netcdf-bin
#------------------------------------------------------------------------------
# NetCDF-Fortran (netcdf-fortran-4.4.3)
sudo apt install libnetcdff
#==============================================================================
#IOAPI (ioapi-3.2.0):
git clone https://github.com/cjcoats/ioapi-3.2
cd ioapi-3.2

cp Makefile.template Makefile
cp ioapi/Makefile.nocpl ioapi/Makefile
cp m3tools/Makefile.nocpl m3tools/Makefile

#Agregar al Makefile:
```
NCFLIBS ="-L$(nc-config --libdir) -lnetcdff -lnetcdf"
```
# Comentar cualquier flag de openMPI en los "Makeinclude.Linux2_x86_64gfort" si es que los hay.
```
OMPFLAGS = # -fopenmp
OMPLIBS = # -fopenmp
```

make configure BIN=Linux2_x86_64gfort CPLMODE=nocpl BASEDIR=/home/usuario/m/libs/ioapi-3.2 
make BIN=Linux2_x86_64gfort CPLMODE=nocpl BASEDIR=/home/usuario/m/libs/ioapi-3.2 

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
