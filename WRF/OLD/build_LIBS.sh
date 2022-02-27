################################
# Compilacion de WRF-4.0
################################
#My build directory tree:   ---------------------------------------------------
# build/
#	libs/  (libraries source)
#	lib/   (compiled libraries)
#	WRF/   (WRF source)
#	WPS/   (WPS source)       
# -----------------------------------------------------------------------------
#Librariess utilized (they are all in /work/mh0735/m300734/tars )
#curl-7.61.1.tar.gz  jasper-1.900.1.tar.gz  netcdf-4.4.0.tar.gz          parallel-netcdf-1.7.0.tar.gz  
#byacc.tar.gz        hdf5-1.8.16.tar  	   openmpi-1.10.2.tar.gz        zlib-1.2.7.tar.gz
#flex.tar.gz         libpng-1.2.50.tar.gz   netcdf-fortran-4.4.3.tar.gz  szip-2.1.1.tar.gz
# -----------------------------------------------------------------------------
#Cargar compiladores (intel-18.0.4)
module purge
module load intel
# -----------------------------------------------------------------------------
# En la carpeta "lib" voy a guardar todas las librerias.
export DIR=/work/mh0735/m300734/build/lib

# OpenMPI   (openmpi-1.10.2)
./configure --prefix=$DIR/ompi
	#(!) Se puede extender para agregar pmi, mallox, fca, etc.
make -j 16 all
make -j 16 install

export MPI_PATH=$DIR/ompi
export PATH=$PATH:$MPI_PATH/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MPI_PATH/lib
# -----------------------------------------------------------------------------
# Miscellaneous libs (i will put them into the folder lib/grib2 )
mkdir $DIR/grib2;mkdir $DIR/grib2/include;mkdir $DIR/grib2/lib;
# Compilar ZLIB (zlib-1.2.7)
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 FCFLAGS='-m64' FFLAGSC='-m64 -fPIC' LDFLAGS='-L/work/mh0735/m300734/build/lib/grib2/lib' CPPFLAGS='-I/work/mh0735/m300734/build/lib/grib2/include' ./configure --prefix=$DIR/grib2
make;make install

# SZIP  LIBPNG  JASPER  ( szip-2.1.1 , libpng-1.2.50,  jasper-1.900.1)
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 FCFLAGS='-m64 -fPIC' FFLAGSC='-m64 -fPIC' ./configure --prefix=$DIR/grib2
make;make install

export GRIB2=$DIR/grib2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRIB2/lib
export PATH=$PATH:$GRIB2/bin
# -----------------------------------------------------------------------------
# HDF5 (hdf5-1.8.16)
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 CXXFLAGS='-g -O2 -fPIC' CFLAGS='-g -O2 -fPIC' FFLAGS='-g -fPIC' FCFLAGS='-g -fPIC' LDFLAGS='-fPIC' F90LDFLAGS='-fPIC' FLDFLAGS='-fPIC' ./configure --prefix=$DIR/hdf5 --enable-parallel -enable-shared --with-szlib='/work/mh0735/m300734/build/lib/grib2' --with-zlib='/work/mh0735/m300734/build/lib/grib2'

make -j 28 install

export HDF5=$DIR/hdf5
export PATH=$PATH:$HDF5/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib
# -----------------------------------------------------------------------------
# Parallel-NetCDF: ( parallel-netcdf-1.7.0)
CC=icc CXX=icpc FC=ifort F90=ifort F77=ifort MPICC=mpicc MPICXX=mpicxx MPIFC=mpif90 MPIF77=mpif77 MPIF90=mpif90 CXXFLAGS='-g -O2 -fPIC' CFLAGS='-g -O2 -fPIC' FFLAGS='-g -fPIC' FCFLAGS='-g -fPIC' LDFLAGS='-fPIC'F90LDFLAGS='-fPIC' FLDFLAGS='-fPIC' ./configure --prefix=$DIR/parallel-netcdf --enable-fortran --enable-large-file-test --enable-largefile

make -j 28 install
	#(!) Creo que me tiro fatal error, por que no podia acceder a una carpeta
	#dentro de --prefix que no estaba creada, y al volver a poner make ya siguio
	# hasta el final.

export PNET=$DIR/parallel-netcdf
export PATH=$PATH:$PNET/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PNET/lib
# -----------------------------------------------------------------------------
# NetCDF-C (netcdf-4.4.0)
export OMPI_MPICC=icc
export OMPI_MPICXX=icpc
export OMPI_MPIFC=ifort
export OMPI_MPIF90=ifort

export WRFIO_NCD_LARGE_FILE_SUPPORT=1

CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort CPPFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' CXXFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' CFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' FFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' FCFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' LDFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/parallel-netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include -L/work/mh0735/m300734/build/lib/hdf5/lib -L/work/mh0735/m300734/build/lib/parallel-netcdf/lib -L/work/mh0735/m300734/build/lib/grib2/lib' ./configure --prefix=$DIR/netcdf --enable-fortran --disable-static --enable-shared --with-pic --enable-parallel-tests -enable-pnetcdf --enable-large-file-tests --enable-largefile
	#(!)Tuve que poner los paths completos, por que anda mal el parser de paths.
make; make install

export NCDIR=$DIR/netcdf
export PATH=$PATH:$NCDIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NCDIR/lib
# -----------------------------------------------------------------------------
# NetCDF-Fortran (netcdf-fortran-4.4.3)

CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 MPICC=icc MPICXX=icpc MPIFC=ifort MPIF77=ifort MPIF90=ifort CPPFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' CXXFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' CFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' FFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' FCFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include' LDFLAGS='-I/work/mh0735/m300734/build/lib/hdf5/include -I/work/mh0735/m300734/build/lib/netcdf/include -I/work/mh0735/m300734/build/lib/grib2/include -L/work/mh0735/m300734/build/lib/hdf5/lib -L/work/mh0735/m300734/build/lib/netcdf/lib -L/work/mh0735/m300734/build/lib/grib2/lib' ./configure --prefix=$DIR/netcdf --disable-static --enable-shared --with-pic --enable-parallel-tests --enable-large-file-tests --enable-largefile
	#(!)Tuve que poner los paths completos, por que anda mal el parser de paths.
make
make install

# -----------------------------------------------------------------------------
# BYACC,   FLEX,   CURL  (byacc-20120115 , flex-2.5.3, curl-7.61.1)
CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77 ./configure --prefix=$DIR/<LIBRERIA>
make; make install

export PATH=$PATH:$DIR/byacc/bin:$DIR/flex/bin:$DIR/curl/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/byacc/lib:$DIR/flex/lib:$DIR/curl/lib
# -----------------------------------------------------------------------------
# Compilar WRF-4.0.3

./clean # (si ya se intentó compilar antes)

#libs:
#export MPICH=$DIR/ompi
export HDF5=$DIR/hdf5
export HDF5_PATH=$DIR/hdf5
export PHDF5=$DIR/hdf5
export NETCDF=$DIR/netcdf
export PNETCDF=$DIR/parallel-netcdf

export ZLIB=$DIR/grib2
export ZLIB_PATH=$DIR/grib2
export LIBPNG=$DIR/grib2
export JASPERINC=$DIR/grib2/include
export JASPERLIB=$DIR/grib2/lib
export SLIB=$DIR/grib2

export YACC=$DIR/byacc/bin/yacc
export CURL_PATH=$DIR/curl
export FLEX=$DIR/flex
export FLEX_LIB_DIR=$DIR/flex/lib

#config:
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export NETCDF4=1
unset NETCDF_classic
export WRF_EM_CORE=1
export WRF_NMM_CORE=0
#chem:
export WRF_CHEM=1
export WRF_KPP=1

./configure 
	# elegir opt=15 (dmpar) INTEL (ifort/icc), y luego opt=1.

        #(!)Si NO detecta NetCDF:
        # vim configure.wrf  #Modificar para que quede:
        #	if [ “`uname`” = “Darwin” ] ; then
        #	ans=“`whereis nf-config`”
        #	elif [ “`uname`” = “Linux” ] ; then
        #	ans=“`which nf-config`”
        #	else
./compile em_real >& log.compile

# -----------------------------------------------------------------------------
# Compilar WPS-4.0.3
./clean -a
./configure
	#elegir 17 (serial) INTEL (ifort/icc)
./compile
        #(!)Si NO compila ungrib por problema con namelist y allocatables.
        # 	vim /ungrib/src/read_namelist.F
        # 	 Cambiar linea 56 por:  real, dimension(250) :: new_plvl
        # 	 Eliminar linea  74
        # 	 Eliminar linea 297

        #(!)Si dice que no encuentra el directorio de WRF.
        # 	vim configure.wps #(despues de hacer configure):
        #        ifneq ($(wildcard $(DEV_TOP)/../WRF), ) # Check for WRF v4.x directory
        #           WRF_DIR         =       ../WRF
        #         else
        #           WRF_DIR         =       ../WRFV3
        #        endif
        #        WRF_DIR=../WRFv4 # <== AGREGAR ESTA LINEA!

	#(!)Si no compila geogrid y metgrid y da errores del tipo: wrf_io.f  undefined reference to ...  <algo con netcdf>.
	#	1) no esta linkeado netcdf de fortran => 
	#		vim configure.wps
	#		WRF_LIB= ... 
	# 			L$(NETCDF)/lib  -lnetcdf -lnetcdff  #<== AGREGAR -lnetcdff
	#
	#	2) problemas con compilacion ompi     => 
	#		vim metgrid|geogrid/src/Makefile:
	#		$(MPI_LIB) -qopenmp -fpp -auto   #(!) si pones -openmp te pide que lo cambies a -qopenmp
# -----------------------------------------------------------------------------
