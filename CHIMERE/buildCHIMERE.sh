#================================================
# CHIMERE: Build & Run
#================================================
#Dependencias:
#   * C y Fortran Compilers
module load intel/intel_2015.3.187
#   * MPI Library (OpenMPI)
module load ompi3.1.4_intel_2015.3.187
#   * HDF5 with parallelism enabled
module load hdf5_1.8.17_intel_2015.3.187
#   * NetCDF Library
module load netcdf4.4_intel_2015.3.187_parallel
#   * grib2, jasper blitz
module load miscLibs_intel_2015.3.187
#------------------------------------------------
#(0) Traer Repositorio:
tar -xzv stage/chimere_v2020r3.tar.gz

cd chimere_v2020r3
#------------------------------------------------
#(1) Setiar variables en script mychimere/mychimere.sh
cp mychimere/mychimere-ciclad.ifort mychimere/mychimere.sh  #agarré el m�s parecido a lo que necesito
#lo edité: para neurus lo llame mychimer-neurusS.sh (deberia haber una copia en la misma carpeta de este script)

./build-chimere.sh --avail            #ver compilaciones disponibles
./build-chimere.sh --mychimere neurus #compilar

#------------------------------------------------
#(2) RUN



