#================================================
# CMAQ: Build & Run
#================================================
#Dependencias de CMAQ:
#   * C y Fortran Compilers  (GNU> 6.1 | Intel > 17.0)
module purge
module load gcc/6.3.0
#   * MPI Library (IntelMPI>2017.0 | MPICH>3.3.1 | MVAPICH2>2.3.1 | OpenMPI>2.1.0)
module load ompi3.1.4_gcc_6.3.0    
#   * NetCDF-C y NetCDF-Fortran (NetCDF-C>4.2 | NetCDF-Fortran>4.4.2) (!) sin HDF4, HDF5,DAP, PnetCDF ni Zlib 
module load netcdf4.4_gcc_6.3.0
#   * I/O API Library
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
setenv NCDIR  /home/ramiroespada/libs_gcc_6.3.0/netcdf
setenv NFDIR  /home/ramiroespada/libs_gcc_6.3.0/netcdf
setenv NETCDF netcdf_combined_directory_path 
setenv IOAPI  /home/ramiroespada/libs_gcc_6.3.0/ioapi-3.2/Linux2_x86_64gfort
setenv WRF_ARCH 34                              # [1-75] Optional, ONLY for WRF-CMAQ  
#> I/O API, netCDF, and MPI library locations
setenv IOAPI_INCL_DIR ${IOAPI}  #> I/O API include header files
setenv IOAPI_LIB_DIR  ${IOAPI}  #> I/O API libraries

setenv MPI_INCL_DIR /home/ramiroespada/libs_gcc_6.3.0/mpich/include #> MPI Include directory path
setenv MPI_LIB_DIR  /home/ramiroespada/libs_gcc_6.3.0/mpich/lib     #> MPI Lib directory patha
```
#(!) revisar que todo est√©bien setiado
# Ejecutar script:
csh config_cmaq.csh gcc

#------------------------------------------------
#(3) Ir a CCTM/scripts
cd CCTM/scripts
csh bldit_cctm.csh gcc

# Si todo anda bien: se crea una carpeta "BLD_CCTM_v533_gcc" adentro est· el ejecutable
cd BLD_CCTM_v533_gcc

#(!)ERROR, no puede abrir STATE3.EXT
#Resulta que STATE3.EXT (y otros archivos) estan en ../ioapi-3.2/ioapi/fixed_src  
#Hay que copiarlos dentro de ioapi-3.2/Linux2_x86_64gfort  para que funcione
#------------------------------------------------
#(4) RUN

### Requerimientos:
#- Output de modelo meteorol√gico regional, ej: WRF
#- Meteorology-Chemistry Interface Processor (MCIP)
#- Emisiones, procesados con el soft: Sparse Matrix Operator Kernel for Emissions (SMOKE)
#- Init cond (ICON)
#- Boundary cond (BCON)


#El CMAQ (CCTM). MCIP, ICON and BCON estan incluidos en el repo de CMAQ.
#Mientras que SMOKE y FEST-C y "Spatial allocator tools" son softs externos.

### Inputs:

