#================================================
# SILAM: Build & Run
#================================================
#Dependencias:
#   * Fortran Compilers (gfortran)
#   * MPI Library (OpenMPI)
#   * HDF5 Library
#   * NetCDF Library
#   * GRIB_API (usa ECCODES)
#   * Otros: (LAPACK/BLAS, PROJ)
#------------------------------------------------
#(0) Instalar dependencias:
sudo apt install make python gfortran libeccodes-dev libnetcdf-dev libnetcdff-dev liblapack-dev libblas-dev libbz2-dev libproj-dev
#(1) Traer Repositorio:
git clone https://github.com/fmidev/silam-model
#------------------------------------------------
#(2) Compilar
cd silam-model/source/
make gnu
make
#Deberia crearse un binario en:
#$ ../bin/silam_v5_7pub.gnu
#Si hay problemas con alguna libreria cambiar flags en ../build/options.gnu y volver a intentar.
#En ubuntu hay un error: "grib_api.mod not found" que en el readme.md del repo explican como resolver.
