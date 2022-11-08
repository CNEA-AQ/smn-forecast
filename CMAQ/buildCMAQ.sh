#================================================
# CMAQ v5.4: Build & Run
#================================================
#Dependencias de CMAQ:
#   * C y Fortran Compilers  (GNU> 6.1 | Intel > 17.0)
#   * MPI Library (IntelMPI>2017.0 | MPICH>3.3.1 | MVAPICH2>2.3.1 | OpenMPI>2.1.0)
#   * NetCDF-C y NetCDF-Fortran (NetCDF-C>4.2 | NetCDF-Fortran>4.4.2) (!) sin HDF4, HDF5,DAP, PnetCDF ni Zlib 
#   * I/O API Library
#------------------------------------------------
#(0) Traer Repositorio:
git clone https://github.com/USEPA/CMAQ
cd CMAQ
#------------------------------------------------
#(1) Setiar variables en script bldit_project.csh
``` set CMAQ_HOME = /home/usuario/m/CMAQ ```
#Ejecutar:
tcsh bldit_project.csh
#------------------------------------------------
#(2) Setiar variables en scripts config_cmaq.csh
```
        setenv NETCDF     /usr #netcdf_root_gcc # Note please combine netCDF-C & Fortran Libraries
        setenv IOAPI      /home/usuario/m/libs/ioapi-3.2/Linux2_x86_64gfort #ioapi_root_gcc
        setenv WRF_ARCH   34 # [1-75]

        #> I/O API, netCDF, and MPI Library Locations -- used in CMAQ
        setenv IOAPI_INCL_DIR   ${IOAPI}             #> I/O API include header files
        setenv IOAPI_LIB_DIR    ${IOAPI}             #> I/O API libraries
        setenv NETCDF_LIB_DIR   /usr/lib/x86_64-linux-gnu                 #> netCDF C directory path
        setenv NETCDF_INCL_DIR  /usr/include         #> netCDF C directory path
        setenv NETCDFF_LIB_DIR  /usr/lib/x86_64-linux-gnu          #> netCDF Fortran directory path
        setenv NETCDFF_INCL_DIR /usr/include         #> netCDF Fortran directory path
        setenv MPI_LIB_DIR      /lib/x86_64-linux-gnu    			#> MPI Lib directory path
        setenv MPI_INCL_DIR     /lib/x86_64-linux-gnu/mpich/include             #> MPI Include directory path
        #setenv MPI_INCL_DIR     /lib/x86_64-linux-gnu/openmpi/include          #> MPI Include directory path
```
#(!) revisar que todo esta bien setiado
#(!) si ioapi está compilado con mpich, usar ese.
# Ejecutar script con csh ó tcsh:
tcsh config_cmaq.csh gcc
#(!) Asegurarse de setiar el compilador correspondiente al que se quiere usar (por ejemplo: setenv myFC mpifort.mpich)

#------------------------------------------------
#(3) Ir a CCTM/scripts y ejecutar bldit_cctm.csh gcc con csh ó tcsh
cd CCTM/scripts
tcsh bldit_cctm.csh gcc

# Si todo anda bien: se crea una carpeta "BLD_CCTM_v533_gcc" que contiene el ejecutable

#(!)ERROR, no puede abrir STATE3.EXT
#Resulta que STATE3.EXT (y otros archivos) estan en ../ioapi-3.2/ioapi/fixed_src  
#Hay que copiarlos dentro de ioapi-3.2/Linux2_x86_64gfort  para que funcione
#------------------------------------------------
#(4) Compilar preprocesadores

# (a) Compilar ICON (initial conditions)
#Ir a PREP/mcip/src
cd ~/CMAQ/PREP/icon/scripts
tcsh bldit_icon.csh gcc

# (b) Compilar BCON (boundary conditions)
#Ir a PREP/mcip/src
cd ~/CMAQ/PREP/bcon/scripts
tcsh bldit_icon.csh gcc

# (c) Compilar MCIP (preprocesador de met)
#Ir a PREP/mcip/src
cd ~/CMAQ/PREP/mcip/src

#Editar Makefile:
```
#...gfortran
FC  = gfortran
NETCDF=/usr
IOAPI_ROOT=/home/usuario/m/libs/ioapi-3.2
FFLAGS  = -O3 -I$(NETCDF)/include -I$(IOAPI_ROOT)/Linux2_x86_64gfort
LIBS    = -L$(IOAPI_ROOT)/Linux2_x86_64gfort -lioapi  \
          -L$(NETCDF)/x86_64-linux-gnu -lnetcdff -lnetcdf
```
tcsh #activo tcsh
source ~/CMAQ/config_cmaq.csh
make

#------------------------------------------------
### Corrida de Benchmark:
#Descargar datos de prueba: CMAQv5.3.2_Benchmark_2Day_Input.tar.gz
#Ir a ~/CMAQ/CCTM/scripts
#copiar sbatch.cmaq editar (sobre todo la variable CMAQ_DATA que es donde se van a guardar los outs y se van a leer los inps)
 
module purge
module load cmaq5.3.3_gcc_6.3.0
sbatch sbatch.cmaq

#================================================
#(4) RUN

### Datos requeridos:
#- Output de modelo meteorologico regional, ej: WRF  (si es wrf4.x ó wrfchem hay que editarle el TITLE en el header de los wrfout)
#- Emisiones, procesados con el soft: Sparse Matrix Operator Kernel for Emissions (SMOKE)

# Preprocesadores:
#- Meteorology-Chemistry Interface Processor (MCIP)
#- Initial conditions (ICON)
#- Boundary conditions (BCON)

#El CMAQ (CCTM). MCIP, ICON and BCON estan incluidos en el repo de CMAQ.
#Mientras que SMOKE y FEST-C y "Spatial allocator tools" son softs externos.

## Pasos:

#(1) correr WRF con la meteo.


#(2) procesar los wrfout y geo_em.d01.nc con el MCIP. La salida sirve para el CMAQ y el SMOKE.
cd PREP/mcip/scripts
tcsh run_mcip.csh gcc   #Edito archivo run_mcip.csh para incluir paths, y fecha de corrida.

#(!) al correr, quiza no le gusta el TITLE del wrfout. Se puede cambiar usando: ncatted -O -h -a TITLE,global,m,c," OUTPUT FROM WRF V4    MODEL" wrfout wrfout_modif

#(3) crear archivos de emisiones con SMOKE.
# Hice un script con python: prepEmis.py para este fin (Hay que trabajarlo).

#(4) crear archivo de condicion de borde con BCON. Vamos a necesitar un modelo global.
cd PREP/bcon/scripts
#Editar archivo run_icon.csh para incluir paths, fecha de corrida, compilador, nombre de grilla y ubicacion de GRIDDESC. Por el momento IBTYPE=profile.
tcsh run_bcon.csh gcc

#(5) crear archivo de condiciones iniciales ICON. Vamos a necesitar un modelo global.
cd PREP/icon/scripts
#Editar archivo run_icon.csh para incluir paths, fecha de corrida, compilador, nombre de grilla y ubicacion de GRIDDESC. Por el momento ICTYPE=profile.
tcsh run_icon.csh gcc

#(6) Correr CMAQ:
cd CCTM/scripts

tcsh run_cctm.sh  #Hay varios scripts que preparan la corrida. editar las variables y correr.


