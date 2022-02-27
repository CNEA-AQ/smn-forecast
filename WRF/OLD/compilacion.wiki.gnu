###############################################################
##  C O M PI L A C I O N   DE  L I Bs &  P R O G R A M A S   ##
####======================================================#####

#INICIO:
mkdir /home/ramiroespada/WRF/Build_WRF
#========================================================================================
## LIBRERIAS:
mkdir /home/ramiroespada/libs
DIR= /home/ramiroespada/libs
#----------------------------------------------------------------------------------------
#DESCOMPRIMIR:
tar -xvf nombre_comprimido.tar.gz
#-----------------------------------------------------------------------------------------
#mpich
cd ~/libs/mpich-3.0.4
CC=gcc CXX=g++ FC=gfortran FCFLAGS=-m64 F77=gfortran F90=gfortran FFLAGSC=-m64 ./configure --prefix=$DIR/mpich
make;make install
#
#zlib
cd ~/libs/zlib-1.2.7
CC=gcc CXX=g++ FC=gfortran FCFLAGS=-m64 F77=gfortran FFLAGSC=-m64 LDFLAGS=-LARIES/grib2/lib CPPFLAGS=-I$DIR/grib2/include ./configure --prefix=$DIR/grib2
make;make install
#libpng
cd ~/libs/libpng-1.2.50
CC=gcc CXX=g++ FC=gfortran FCFLAGS=-m64 F77=gfortran ./configure --prefix=$DIR/grib2
make;make install
#jasper
cd jasper-1.900.1
CC=gcc CXX=g++ FC=gfortran FCFLAGS=-m64 F77=gfortran ./configure --prefix=$DIR/grib2
make;make install
#szip 
tar -zxvf szip-2.1.tar.gz
mkdir -p ../szip-2.1
cd szip-2.1
FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ ./configure --prefix=$DIR/szip-2.1
make; make check ; make install
# jpeg-6b
tar -zxvf jpegsrc.v6b.tar.gz
mkdir -p ../jpeg-6b/bin ../jpeg-6b/man/man1 ../jpeg-6b/lib ../jpeg-6b/include
cd jpeg-6b
FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ ./configure --prefix=$DIR/jpeg-6b
make; make check ; make install; make install-lib

#yacc
cd byacc
CC=gcc CXX=g++ FC=gfortran ./configure --prefix=$DIR/byacc
make;make check;make install
#flex
cd flex
CC=gcc CXX=g++ FC=gfortran ./configure --prefix=$DIR/flex
make;make check;make install
#curl
cd curl-7.62.0
CC=gcc CXX=g++ FC=gfortran ./configure --prefix=$DIR/curl-7.62.0
make;make check,;make install
#================================================================
#LIBRERIAS DIFICILES:
#HDF5
tar -zxvf hdf5-1.8.16.tar.gz
cd hdf5-1.8.16/
#FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ ./configure --prefix=$DIR/src/LIBS/hdf5-1.8.16 --with-zlib=$DIR/src/LIBS/zlib-1.2.8 --with-szlib=$DIR/src/LIBS/szip-2.1 --with-jpeg=$DIR/src/LIBS/jpeg-6b --enable-fortran --disable-netcdf --enable-fortran2003
FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ ./configure --prefix=$DIR/hdf5 --with-zlib=$DIR/zlib-1.2.7 --with-szlib=$DIR/szip-2.1 --enable-fortran --enable-hl
make;make check;make install
# NETCDF 4.4 ------------------------#
tar -zxvf netcdf-4.4.1.1.tar.gz       # Este es para C
tar -zxvf netcdf-fortran-4.4.4.tar.gz # Este es para FORTRAN
mkdir -p ../netcdf-4.4 		      # A ambos los pongo en la misma carpeta
# C:
cd netcdf-4.4.1.1
CC=gcc CXX=g++ ./configure --prefix=$DIR/netcdf4.4 --with-zlib=$DIR/zlib-1.2.8 --disable-netcdf-4 --enable-static --disable-shared LIBS="-lm"
make; make check ; make install
#FORTRAN
cd netcdf-fortran-4.4.4
        # El netcdf-fortran tiene como dependencia al netcdf-C (recien instalado), hay que linkiarlo, eso se puede hacer exportando en el .bashrc --> #LD_LIBRARY_PATH=$DIR/LIBRARIES/netcdf-4.4/lib:$LD_LIBRARY_PATH
	#FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ CPPFLAGS=-I$DIR/LIBRARIES/netcdf-4.4/include LDFLAGS=-L$DIR/LIBRARIES/netcdf-4.4/lib ./configure --prefix=$LIB/netcdf-4.4 --with-zlib=$LIB/zlib-1.2.3 --enable-netcdf-4--enable-static --disable-shared LIBS="-lm"
make; make check ; make install
FC=gfortran F77=gfortran F90=gfortran CC=gcc CXX=g++ HDF5_DIR=$DIR/hdf5 CPPFLAGS="-I/home/ramiroespada/libs/netcdf4.4/include" LDFLAGS="-L/home/ramiroespada/libs/netcdf4.4/lib" ./configure --prefix=$DIR/netcdf4.4 LIBS="-lm"
#
#==================================================================
## WRF
tar -xvf WRF...
cd WRF

./clean # (si ya se intentó compilar antes)

DIR=/home/ramiroespada/libs/
export NETCDFPATH=$DIR/netcdf4.4
export NETCDF=$DIR/netcdf4.4
export HDF5=$DIR/hdf5
export MPICH=$DIR/mpich
export ZLIB=$DIR/zlib-1.2.7
export SLIB=$DIR/szip-2.1
export YACC=$DIR/byacc
export FLEX_LIB_DIR=$DIR/flex/lib
export JASPERINC=$DIR/grib2/include
export JASPERLIB=$DIR/grib2/lib

export PATH=$MPICH/bin:$NETCDF/bin:$HDF5/bin:$PATH
export LD_LIBRARY_PATH=$MPICH/lib:$NETCDF/lib:$HDF5/lib$LD_LIBRARY_PATH

export NETCDF4=1 
unset NETCDF_classic
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
#SIN QUIMICA:
export J="-j 4"
export WRF_CHEM=0
export WRF_KPP=0
export WRF_EM_CORE=1
export WRF_NMM_CORE=0

#CON QUIMICA:
 export J="-j 4"
 export WRFIO_NCD_LARGE_FILE_SUPPORT=1
 export WRF_CHEM=1
 export WRF_KPP=0
 export WRF_EM_CORE=1
 export WRF_NMM_CORE=0

./configure # elegir opt=34, y luego opt=1

	#SI NO DETECTA NETCDF:
	# modificar archivo de configure:
	#if [ “`uname`” = “Darwin” ] ; then
	#ans=“`whereis nf-config`”
	#elif [ “`uname`” = “Linux” ] ; then
	#ans=“`which nf-config`”
	#else
./compile em_real >& log.compile

#========================================================================================
## WPS
tar -xvf WPS...
cd WPS
./clean 
JASPERINC=$DIR/grib2/include JASPERLIB=$DIR/grib2/lib ./configure
./compile

	#SI NO COMPILA UNGRIB: (es un problema con namelist y allocatables)
	#vim /ungrib/src/read_namelist.F
	# Cambiar linea 56 por:  real, dimension(250) :: new_plvl
	# Eliminar linea  74
	# Eliminar linea 297
#========================================================================================
# EXPORT LIBRERIAS:
export NETCDF=$DIR/netcdf4.4
export HDF5=$DIR/hdf5
export MPICH=$DIR/mpich
export ZLIB=$DIR/zlib-1.2.7
export SLIB=$DIR/szip-2.1
export YACC=$DIR/byacc

export PATH=$MPICH/bin:$NETCDF/bin:$HDF5/bin:$PATH
export LD_LIBRARY_PATH=$MPICH/lib:$NETCDF/lib:$HDF5/lib$LD_LIBRARY_PATH
