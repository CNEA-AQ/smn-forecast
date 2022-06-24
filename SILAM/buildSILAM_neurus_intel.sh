#================================================
# SILAM: Build & Run
#================================================
#Dependencias:
#   * Fortran Compilers
module load intel/intel_2015.3.187
#   * MPI Library (OpenMPI)
module load ompi3.1.4_intel_2015.3.187
#   * HDF5 Library
module load hdf5_1.8.17_intel_2015.3.187
#   * NetCDF Library
module load netcdf4.4_intel_2015.3.187_parallel
#   * GRIB_API (us√ ECCODES en su lugar)
module load eccodes-2.23_intel_2015.3.187
#   * OTROS (porlasdudas)
module load miscLibs_intel_2015.3.187

#------------------------------------------------
#(0) Traer Repositorio:
module load http_proxy
git clone https://github.com/fmidev/silam-model

cd SILAM/source
#------------------------------------------------
#(1) Setiar variables en Makefile
SE=intel #gnu
DefaultSetup=intel #DefaultSetup
#------------------------------------------------
#(2) Setiar variables en ../build/options.intel
# Build platform
ARCH = linux_intel
# The compiler
F90C = mpifort
# Flags for optimization, etc.
OPTIMIZATION = -O2 #-vec-report2
OPENMP = -openmp
DEBUG = -traceback -diag-disable 8291 -check bounds #-g -check all 
FIXED = -fixed
PREPROCESS = -cpp
# Include flags
INCLUDE = -I$(OBJDIR) -I$(NCDIR)/include -I$(H5DIR)/include -I$(GRDIR)/include -L$(NCDIR)/lib -L$(H5DIR)/lib -L$(GRDIR)/lib
FFLAGS = -assume byterecl  -static $(OPTIMIZATION) $(OPENMP) $(DEBUG) $(INCLUDE)
#hdf5 install directory
H5DIR=/home/ramiroespada/libs_intel_2015.3.187/hdf5
#netcdf install directory
NCDIR =/home/ramiroespada/libs_intel_2015.3.187/pnetcdf
# grib api install directory
GRDIR =/home/ramiroespada/libs_intel_2015.3.187/eccodes-2.23-Source/build
#GRDIR =/home/ramiroespada/libs_intel_2015.3.187/grib_api-1.15.0-Source/build
#default directory with libraries
LIBDIR =/home/ramiroespada/libs_intel_2015.3.187/grib2/lib
#Math Kernel Library (ifort)
MKLROOT =/share/apps/intel/composer_xe_2015.3.187/mkl
# intel mkl -- oh!
MKL = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
SILAM_LIBS = $(LIBDIR)libfish.a  -L$(H5DIR)/lib -L$(NCDIR)/lib -L$(GRDIR)/lib \
  -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lgrib_api_f90  -lgrib_api -ljasper $(MKL) -lz -openmp #-lblas 
#------------------------------------------------
#(3) MAKE
make

#Error #1: md.silja.mod.linux_intel.f90
./source/md.silja.mod.linux_intel.f90(399): error #6404: This name does not have a type, and must have an explicit type.   [GETPAGESIZE]
my_pagesize = getpagesize()
#Soluci√n: Falta en agregar antes del CONTAINS de md.silja.mod.linux_intel.f90:
# interface
#       integer(c_int) function getpagesize () bind (c, NAME='getpagesize')
#       use, intrinsic :: ISO_C_BINDING
#       end function
# end interface
#---------------------
#Error #2: grib_api_io.silam.mod.f90
../source/grib_api_io.silam.mod.f90(414): error #6364: The upper bound shall not be omitted in the last dimension of a reference to an assumed size array.   [PRECISION4EDITION]
    real, dimension(*), parameter :: precision4edition = (/1e-3, 1e-6, 1e-6/) !!lon and lat precision
#Solucion: Tengo que cambiarlo por:
    real, dimension(3), parameter :: precision4edition = (/1e-3, 1e-6, 1e-6/) !!lon and lat precision
#---------------------
#Error #3: grib_api_io.silam.mod.f90
../source/grib_api_io.silam.mod.f90(592): error #6446: This name has already been used as a construct name.   [GRIB_GRID_CACHE]
              grGridPtr => grib_grid_cache(iCacheInd)%grGrid
#Solucion: Cambiar el nombre de $OMP CRITICAL(grib_grid_cache) y !$OMP END CRITICAL (grib_grid_cache_1)  por otro, por ejemplo: $OMP CRITICAL(grib_grid_cache_1)
#---------------------
#Error #4: advection_eulerian_v5.silam.mod.f90
#(a)
../source/advection_eulerian_v5.silam.mod.f90(2680): error #8527: If a bound remapping list is specified, data target must be simply contiguous or of rank one.   [PASSENGERS]
                  passDiffBak(0:maxpass,0:maxz) => mystuff%passengers(:,:,iSpecies,iSrc)
                  -----------------------------------------------------------^

#(b)
../source/advection_eulerian_v5.silam.mod.f90(2684): error #8527: If a bound remapping list is specified, data target must be simply contiguous or of rank one.   [PASSENGERS_OUT]
                                    passDiffBak(0:maxpass,0:maxz) => mystuff%passengers_out(:,:)
#Solucion:
#El error de (b) se resuelve simplemente comentando mystuff%passengers_out!(:,:)
#El error de (a) entiendo que otros compialdores saben como manejar este problema. Pero con el intel 2015 todav√a no encontr√© soluci√n.

#---------------------
#Error #5: silam_main.f90

# ISO_FORTRAN_ENV no encuentra funciones. Ya fue, no es miportante lo comento
#Solucion

#---------------------
#Error #6: bzip.f90
#bzip.f90:(.text+0x7c4): undefined reference to `gerror_'

#Solucion
# Comentar una serie de lineas en bzip.f90 que son muy villeras.

#------------------------------------------------
#(4) RUN
