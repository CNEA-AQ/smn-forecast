#================================================
# CMAQ: Build & Run
#================================================
#Dependencias:
#   * C y Fortran Compilers
module load intel/intel_2015.3.187
#   * MPI Library (MPICH)
#module load mpich3.1.4_intel_2015.3.187  #no la cargo pero si le tengo que pasar los paths a CMAQ para que lo use
module load ompi3.1.4_intel_2015.3.187    #la nacesito por que NetCDF y IOAPI fueron compiladas con ompi
#   * NetCDF Library
module load netcdf4.4_intel_2015.3.187_parallel
#   * I/O API Library
module load ioapi-3.2_intel_2015.3.197

#------------------------------------------------
#(0) Traer Repositorio:
module load http_proxy
git clone https://github.com/USEPA/CMAQ

cd CMAQ
#------------------------------------------------
#(1) Setiar variables en script bldit_project.csh
``` set CMAQ_HOME = /home/ramiroespada/CMAQ ```
#Ejecutar:
csh bldit_project.csh

#------------------------------------------------
#(2) Setiar variables en scripts config_cmaq.csh
```
setenv NCDIR  /home/ramiroespada/libs_intel_2015.3.187/pnetcdf
setenv NFDIR  /home/ramiroespada/libs_intel_2015.3.187/pnetcdf
setenv NETCDF netcdf_combined_directory_path 
setenv IOAPI  /home/ramiroespada/libs_intel_2015.3.187/ioapi-3.2/Linux2_x86_64ifortmpi

setenv MPI_INCL_DIR /home/ramiroespada/libs_intel_2015.3.187/mpich/include #> MPI Include directory path
setenv MPI_LIB_DIR  /home/ramiroespada/libs_intel_2015.3.187/mpich/lib     #> MPI Lib directory patha

#setenv MPI_INCL_DIR /home/ramiroespada/libs_intel_2015.3.187/ompi/include #> MPI Include directory path
#setenv MPI_LIB_DIR  /home/ramiroespada/libs_intel_2015.3.187/ompi/lib     #> MPI Lib directory patha
```
#(!) revisar que todo estébien setiado
# Ejecutar script:
csh config_cmaq.csh intel

#------------------------------------------------
#(3) Ir a CCTM/scripts
cd CCTM/scripts

csh bldit_cctm.csh

cd BLD_CCTM_v533_intel/
make

#(!)ERROR, no puede abrir STATE3.EXT
#Resulta que STATE3.EXT y otros archivos estaba en /home/ramiroespada/libs_intel_2015.3.187/ioapi-3.2/ioapi/fixed_src. Los copié en la carpeta que le indique como include al CMAQ. y anduvo!




#(4) RUN


### Requerimientos:
#- Output de modelo meteorol�gico regional, ej: WRF
#- Meteorology-Chemistry Interface Processor (MCIP)
#- Emisiones, procesados con el soft: Sparse Matrix Operator Kernel for Emissions (SMOKE)
#- Init cond (ICON)
#- Boundary cond (BCON)


#El CMAQ (CCTM). MCIP, ICON and BCON estan incluidos en el repo de CMAQ.
#Mientras que SMOKE y FEST-C y "Spatial allocator tools" son softs externos.

### Inputs:
