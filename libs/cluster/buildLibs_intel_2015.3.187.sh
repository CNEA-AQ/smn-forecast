################################
# Compilacion de Librerias:
# con compilador 1ntel_215.3.187
################################
# libs's tars:

#(!) para extraer tars:  $> tar -xvf tarfile.tar.gz
# =============================================================================
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
# wget https://download.open-mpi.org/release/open-mpi/v1.10/openmpi-1.10.2.tar.gz
# tar -xzvf openmpi-1.10.2.tar.gz
# cd openmpi-1.10.2
CC=icc CXX=icpc FC=ifort F90=ifort F77=ifort ./configure --prefix=${LIBSDIR}/ompi
#(!) Se puede extender para agregar pmi, mallox, fca, etc.
make -j 16 all
make -j 16 install

export OMPI=$LIBSDIR/ompi
export PATH=$PATH:$OMPI/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMPI/lib
# =============================================================================
#MPICH (mpich-3.1.4)
CC=icc CXX=icpc FC=ifort ./configure --prefix=${LIBSDIR}/mpich
make
make install
# =============================================================================
#CMake (>3.11)  (!) NO FUNCIONA!
tar -xzv TARs/cmake-3.21.4.tar.gz
cd cmake-3.21.4
./bootstrap  # (!) ERROR!  icpc no tiene soporte c++11
make
make install
#-----------------------------------------------------------------------------
# Miscellaneous libs (i will put them into the folder $LIBSDIR/grib2 folder)
misc_path=$LIBSDIR/grib2  #tiene este nombre por razones histÃ³ricas,
                          #la proxima mejor que compile todo mejor llamarle "misc".
mkdir $LIBSDIR/grib2; mkdir $LIBSDIR/grib2/include;mkdir $LIBSDIR/grib2/lib;
# Compilar en este orden:
# ZLIB (zlib-1.2.7), SZIP (szip-2.1.1),  LIBPNG (libpng-1.2.50),  JASPER  (jasper-1.900.1)
# BYACC (byacc-20120115), FLEX (flex-2.5.3), CURL (curl-7.61.1)
export  LDFLAGS="-L${LIBSDIR}/grib2/lib -I${LIBSDIR}/grib2/include"
export CPPFLAGS="-fPIC -I${LIBSDIR}/grib2/include" 
export CFLAGS="-fPIC -I${LIBSDIR}/grib2/include" 
export  FCFLAGS='' 
export  FFLAGSC='-fPIC'
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 ./configure --prefix=$LIBSDIR/grib2 --enable-shared
make;make install

export MISC=$LIBSDIR/grib2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MISC/lib
export PATH=$PATH:$MISC/bin









#-----------------------------------------------------------------------------
#LAPACK & BLAS (10.0.0)
cp INSTALL/make.inc.intel make.inc
module load python/3.6.1/intel64  #entiendo que es solo usa python para el testing
#cambiar en el Makefile: 
TOPSRCDIR=$LIBSDIR/grib2/lib 
make #se generan los archivos: liblapack.a librefblas.a libtmglib.a en $TOPSRCDIR

#-----------------------------------------------------------------------------
# BLITZ

#git clone https://github.com/blitzpp/blitz
cd blitz
mkdir build
cd build

export CPPFLAGS="-I${LIBSDIR}/grib2/include"
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 cmake .. -DCMAKE_CXX_COMPILER=icpc -DCMAKE_INSTALL_PREFIX=/home/ramiroespada/libs_intel_2015.3.187/blitz/build
cmake .. -DCMAKE_CXX_COMPILER=icpc -DCMAKE_INSTALL_PREFIX=/home/ramiroespada/libs_intel_2015.3.187/grib2
make lib
make install

#(puede que salten libs mal instaldas, pejem faltaba libpapi 5.7.0, la compile, linkie y funciono)
#==============================================================================
# HDF5 (hdf5-1.8.16)
# Dependencias: zlib y szip
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 CXXFLAGS='-g -O2 -fPIC' CFLAGS='-g -O2 -fPIC' FFLAGS='-g -fPIC' FCFLAGS='-g -fPIC' LDFLAGS='-fPIC' F90LDFLAGS='-fPIC' FLDFLAGS='-fPIC' ./configure --prefix=$LIBSDIR/hdf5 --enable-parallel -enable-shared --with-szlib='/home/ramiroespada/libs_intel_2015.3.187/grib2' --with-zlib='/home/ramiroespada/libs_intel_2015.3.187/grib2'

make -j 28 install

export HDF5=$LIBSDIR/hdf5
export PATH=$PATH:$HDF5/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib
# -----------------------------------------------------------------------------
# Parallel-NetCDF: (parallel-netcdf-1.10.0)
CC=icc CXX=icpc FC=ifort F90=ifort F77=ifort MPICC=mpicc MPICXX=mpicxx MPIFC=mpif90 MPIF77=mpif77 MPIF90=mpif90 CXXFLAGS='-g -O2 -fPIC' CFLAGS='-g -O2 -fPIC' FFLAGS='-g -fPIC' FCFLAGS='-g -fPIC' LDFLAGS='-fPIC'F90LDFLAGS='-fPIC' FLDFLAGS='-fPIC' ./configure --prefix=$LIBSDIR/pnetcdf --enable-shared --enable-fortran --enable-large-file-test --enable-largefile
make #-j 28 install
make install

export PNET=$LIBSDIR/pnetcdf
export PATH=$PATH:$PNET/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PNET/lib
#==============================================================================
# NetCDF-C (netcdf-4.4.0)
#con pNetCDF:
incs="-I${LIBSDIR}/grib2/include -I${LIBSDIR}/hdf5/include -I${LIBSDIR}/pnetcdf/include"
libs="-L${LIBSDIR}/grib2/lib -L${LIBSDIR}/hdf5/lib -L${LIBSDIR}/pnetcdf/lib"
export CPPFLAGS="${incs} ${libs}" 
export CXXFLAGS="${incs} ${libs}" 
export   CFLAGS="${incs} ${libs}"
export   FFLAGS="${incs} ${libs}" 
export  FCFLAGS="${incs} ${libs}" 
export  LDFLAGS="${incs} ${libs}"
CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort ./configure --prefix=$LIBSDIR/pnetcdf --enable-fortran --enable-shared --with-pic --enable-parallel-tests --enable-large-file-tests --enable-largefile --enable-pnetcdf 
make
make install
#(!) Dice que no soporta pNetCDF.
# Cuando mando el configure dice que necesita una versiÃn de pnetcdf >1.6 pero estoy usando la 1.10. =/
#==============================================================================
# NetCDF-Fortran (netcdf-fortran-4.4.3)
#Usando mismos flags que en NetCDF-C!!
CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort ./configure --prefix=$LIBSDIR/pnetcdf --enable-shared --with-pic --enable-parallel-tests --enable-large-file-tests --enable-largefile
make
make install
##==============================================================================
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

#-----------------------------------------------------------------------------
# grib-api 1.15.0
module load cmake/3.20.5_gcc6.3  

mkdir build
cd build
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 cmake .. -DCMAKE_CXX_COMPILER=icpc -DCMAKE_INSTALL_PREFIX=/home/ramiroespada/libs_intel_2015.3.187/grib2 cmake .. -DCMAKE_INSTALL_PREFIX=/home/ramiroespada/libs_intel_2015.3.187/grib2

#Si no encuentra NetCDF-C hay que poner en el CMakeCache.txt
#NETCDF_netcdf_LIBRARY_RELEASE:FILEPATH=/home/ramiroespada/libs_intel_2015.3.187/pnetcdf/lib/libnetcdf.so
make 
make install

# -----------------------------------------------------------------------------
# EC-CODES
# Dependencias opcionales:
# * NetCDF
# * AEC, eccodes: Support for Adaptive Entropy Coding
# * Jasper
# * OpenJPEG

tar -xvf TARs/eccodes-2.23.0-Source.tar.gz
mkdir eccodes
cd eccodes
module load cmake/3.20.5_gcc6.3   # (!) necesita cmake > 3.11
export LDFLAGS="-L$LIBSDIR/ompi/lib -I$LIBSDIR/ompi/include -L$LIBSDIR/hdf5/lib -I$LIBSDIR/hdf5/include -L$LIBSDIR/pnetcdf/lib -I$LIBSDIR/pnetcdf/include -L$LIBSDIR/grib2/lib -I$LIBSDIR/grib2/include -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz" 
CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort cmake .. -DCMAKE_C_FLAGS="-O2 -Wall -I/home/ramiroespada/libs_intel_2015.3.187/pnetcdf/include -L/home/ramiroespada/libs_intel_2015.3.187/pnetcdf/lib" -DCMAKE_INSTALL_PREFIX=/home/ramiroespada/libs_intel_2015.3.187/eccodes-2.23-Source/build -DENABLE_NETCDF=ON
#luego falla por que no encuentra la lib NetCDF, pero vas a CMake.cache, y editas:
#//netcdf C library
#NetCDF_C_LIBRARY:FILEPATH=/home/ramiroespada/libs_intel_2015.3.187/pnetcdf/lib/libnetcdf.so
#luego repetis el comando anterior
make
ctest
make install
##==============================================================================
