===== Compilar librerias de WRFChem o cargar modulo =====
[[WRFCHEM: Compilar WRFCHEM4.4|Compilar WRFCHEM4.4 con intel]] \\
module load wrfchem/wrfchem4

====== Compilar otras librerias ======
<code> --- //[[mdiaz@cnea.gov.ar]] 2023/09/27 14:38// --- //[[mdiaz@cnea.gov.ar]] 2023/09/27 14:38// --- //[[mdiaz@cnea.gov.ar]] 2023/09/27 14:38//
CC=mpiicc CXX=mpiicpc FC=mpiifort F90=mpiifort F77=mpiifort cmake .. -DCMAKE_C_FLAGS="-O2 -Wall -I/home/mdiaz/pquimica_mdiaz/libs_wrf/netcdf/include -L/home/mdiaz/pquimica_mdiaz/libs_wrf/netcdf/lib" -DCMAKE_INSTALL_PREFIX=/home/mdiaz/pquimica_mdiaz/libs_wrf/eccodes-2.23-Source/build -DENABLE_NETCDF=ON
</code>

====== Compilar chimere =====
<code>
./build-chimere.sh --arch neurus.ifort

./build-wrf.sh --arch neurus.ifort
</code>

===== mychimere.neurus.ifort =====
<code>
#!/bin/bash

#---------------------------------------------------------------------------------
#       Architecture file for compiling and running CHIMERE     
#       Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------
## Modules
module purge
module load wrfchem/wrfchem4netcdf4.4_intel
#module load eccodes-2.23_intel_2015.3.187

export MiDIR=/home/mdiaz/pquimica_mdiaz
export MiLIBS=$MiDIR/libs_wrf
#---------------------------------------------------------------------------------
#       Compilers
#---------------------------------------------------------------------------------
export my_compilerF90=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/ifort #`which ifort`   # Path to Fortran 90 compiler
export my_compilerC=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/icc     #`which icc    `                 # Path to C compiler
export my_compilerCpp=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/icpc  #`which icpc   `               # Path to C++ compiler

#---------------------------------------------------------------------------------
#       MPI - parallel execution of chimere
#---------------------------------------------------------------------------------
export  my_mpiframe=openmpi                                             # implementaion of MPI norm [ ompi / ccrt ] TO REMOVE
export  my_mpibin=/share/apps/intel/impi/5.0.3.048/intel64/bin          # Path to MPI binary directory
export  my_mpirun=/share/apps/intel/impi/5.0.3.048/intel64/bin/mpirun   # Path to mpirun to execute parallel job in MPI
export  my_mpif90=/share/apps/intel/impi/5.0.3.048/intel64/bin/mpiifort # Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=/share/apps/intel/impi/5.0.3.048/intel64/bin/mpiicc    # Wrapper to my_compilerC to link with MPI library
export  my_mpilib=/share/apps/intel/impi/5.0.3.048/intel64/lib          # Path to MPI libraries directory
export  my_mpiinc=/share/apps/intel/impi/5.0.3.048/intel64/include      # Path to MPI include files directory

#export  my_mpibin=$MiLIBS/ompi/bin                     # Path to MPI binary directory
#export  my_mpirun=$MiLIBS/ompi/bin/mpirun   # Path to mpirun to execute parallel job in MPI
#export  my_mpif90=$MiLIBS/ompi/bin/mpifort  # Wrapper to my_compilerF90 to link with MPI library
#export  my_mpicc=$MiLIBS/ompi/bin/mpicc    # Wrapper to my_compilerC to link with MPI library
#export  my_mpilib=$MiLIBS/ompi/lib                     # Path to MPI libraries directory
#export  my_mpiinc=$MiLIBS/ompi/include         # Path to MPI include files directory

#---------------------------------------------------------------------------------
#       HDF5  - parallel version        
#---------------------------------------------------------------------------------
export my_hdflib=$MiLIBS/hdf5/lib               # Path to HDF5 parallel library directory
export my_hdfinc=$MiLIBS/hdf5/include   # Path to HDF5 parallel include files directory

#---------------------------------------------------------------------------------
#       NETCDF-C  - link with HDF5 parallel 
#---------------------------------------------------------------------------------
export my_netcdfCbin=$MiLIBS/netcdf/bin                 # Path to NETCDF-C (linked with HDF5 parallel) binaries directory 
export my_netcdfClib=$MiLIBS/netcdf/lib         # Path to NETCDF-C (linked with HDF5 parallel) library directory

#---------------------------------------------------------------------------------
#       NETCDF-Fortran  - link with HDF5 parallel and NETCDF-C
#---------------------------------------------------------------------------------
export my_netcdfF90bin=$MiLIBS/netcdf/bin      # PATH to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) binaries  directory
export my_netcdfF90lib=$MiLIBS/netcdf/lib               # Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) library  directory
export my_netcdfF90inc=$MiLIBS/netcdf/include   # Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) include files  directory

#---------------------------------------------------------------------------------
#       GRIB  - link with jasper 
#---------------------------------------------------------------------------------
#eccodes
export my_griblib=$MiLIBS/eccodes-2.23-Source/build/lib64       # Path to GRIB library directory
export my_gribinc=$MiLIBS/eccodes-2.23-Source/build/include     # Path to GRIB include files directory

export my_jasperlib=$MiLIBS/grib2/lib                   # Path to JASPER library directory
export my_jasperinc=$MiLIB/grib2/include                # Path to JASPER include files directory

#---------------------------------------------------------------------------------
# agregados por mi probando si toma jasper
#---------------------------------------------------------------------------------
export FLEX_LIB_DIR=$LIBSDIR/grib2/lib
export YACC="${LIBSDIR}/grib2/bin/yacc -y -d"
export JASPERINC=$LIBSDIR/grib2/include
export JASPERLIB=$LIBSDIR/grib2/lib
#---------------------------------------------------------------------------------
#       BLITZ
#---------------------------------------------------------------------------------
export my_blitzinc=$MiLIBS/grib2/include #/blitz                 # Path to BLITZ include files 
#---------------------------------------------------------------------------------
#       Utilities       
#---------------------------------------------------------------------------------
export my_make=make                             # Path to make 
export my_awk=awk                                         # Path to awk
export my_ncdump=$MiLIBS/netcdf/bin/ncdump              # Path to ncdump

#---------------------------------------------------------------------------------
#       Makefile header needed to compile CHIMERE and WRF 
#            - with this architecture configuration -   
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.ifort-64-ompi                        # Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name=configure.wrf.ifort              # Makefile header to compile WRF in config_wrf directory
export configure_wps_file_name=configure_ifort.wps              # Makefile header to compile WPS in config_wps directory

#---------------------------------------------------------------------------------
#       Export of Shared Library to be available at run time    
#---------------------------------------------------------------------------------
export LD_LIBRARY_PATH=${my_hdflib}:${my_netcdfF90lib}:${my_netcdfClib}:${my_griblib}:${my_mpilib}:${my_mpilib}/ompi:$LD_LIBRARY_PATH
                                                                                                                                            97,1          Bot
</code>

===== sbatch Testcase online =====
Completar los  \_\_FILL\_\_
<code>
#!/bin/bash
#SBATCH -J CHIMERE_DCAP_ONLINE          # Job NAME (human)
#SBATCH -p __FILL__           # Job QUEUE (human)
#SBATCH -A __FILL__
#SBATCH -n 30                    # number of cores (human)
#SBATCH --wait-all-nodes=1       # Do not begin execution until all nodes are ready for use
#SBATCH -o %j.slurm.out          # STDOUT
#SBATCH -e %j.slurm.err          # STDERR
#SBATCH --mem 0                  # 0 = all node memory allocated
#SBATCH --cpu-freq high
#SBATCH --verbose
#SBATCH --threads-per-core=1
############################## MAIN ################################


#Limpia enviroment
module purge

#Modulos propios del usuario
module load wrfchem/wrfchem4netcdf4.4_intel

#Modulo MKL y Compiladores intel
module load intel/intel_2015.3.187

#Modulo Intel MPI 5
module load intel/impi_5.0.3.048

#Modulo para usar solo infiniband
module load intel/impi_fabric_ib

#Modulo para cuando se usa hypertherading
module load intel/impi_pinning_HT

#Modulo para cuando se quiere verbose y stats de intel mpi
module load intel/impi_verbose


#Show user limits
echo
echo "User limits:"
echo
ulimit -a
echo
#Show slurm var:
echo
echo "Slurm enviroment:"
echo
if [ -z "$SLURM_NPROCS" ] ; then
    if [ -z "$SLURM_NTASKS_PER_NODE" ] ; then
      SLURM_NTASKS_PER_NODE=1
    fi
    SLURM_NPROCS=$(( $SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE ))
fi
echo "SLURM_JOBID             = " $SLURM_JOBID
echo "SLURM_JOB_NODELIST      = " $SLURM_JOB_NODELIST
echo "SLURM_JOB_PARTITION     = " $SLURM_JOB_PARTITION
echo "SLURM_MEM_PER_CPU       = " $SLURM_MEM_PER_CPU
echo "SLURM_MEM_PER_NODE      = " $SLURM_MEM_PER_NODE
echo "SLURM_ARRAY_JOB_ID      = " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID     = " $SLURM_ARRAY_TASK_ID
echo "SLURMTMPDIR             = " $SLURMTMPDIR
echo "SLURM_TASKS_PER_NODE    = " $SLURM_TASKS_PER_NODE
echo "SLURM_NTASKS            = " $SLURM_NTASKS
echo "SLURM_NTASKS_PER_CORE   = " $SLURM_NTASKS_PER_CORE
echo "SLURM_NTASKS_PER_NODE   = " $SLURM_NTASKS_PER_NODE
echo "SLURM_NTASKS_PER_SOCKET = " $SLURM_NTASKS_PER_SOCKET
echo "SLURM_JOB_NUM_NODES     = " $SLURM_JOB_NUM_NODES 
echo "SLURM_NNODES            = " $SLURM_NNODES 
echo "SLURM_CPUS_PER_TASK     = " $SLURM_CPUS_PER_TASK
echo "SLURM_NPROCS            = " $SLURM_NPROCS
echo "SLURM_SUBMIT_DIR        = " $SLURM_SUBMIT_DIR
echo "SLURM_SUBMIT_HOST       = " $SLURM_SUBMIT_HOST
echo "SLURM_CPU_FREQ_REQ      = " $SLURM_CPU_FREQ_REQ

echo cd $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

srun --nodes=${SLURM_NNODES} bash -c 'hostname'> $SLURM_JOBID.machines
srun --nodes=${SLURM_NNODES} bash -c 'hostname' | sort -r | uniq > $SLURM_JOBID.mpd.hosts

</code>

===== Modificar en scripts chimere-step2.sh y chimere-meteo.sh =====
**chimere-meteo.sh**
<code>
time ${my_mpirun} -n ${nproc_rea} ${real_exe} || exit 1
</code>
por
<code>
time srun --verbose -A ${accountneurus} -p ${partitionneurus} ${mpiparams} -n ${nproc_rea} ${real_exe} || exit 1      
</code>
y

<code>
time ${my_mpirun} -n ${nproc_rea} ${real_exe} || exit 1
</code>
por
<code>
time srun --verbose -A ${accountneurus} -p ${<code>
time ${my_mpirun} -n ${nproc_rea} ${ndown_exe} || exit 1
</code>
por
<code>
time srun --verbose -A ${accountneurus} -p ${partitionneurus} ${mpiparams} -n ${nproc_rea} ${ndown_exe} || exit 1         
</code>

**chimere-step2.sh**
En offline mode \\
<code>
time srun --verbose --threads-per-core=1 -A ${accountneurus} -p ${partitionneurus} -n ${nproc_chimere} ${mpiparams} ./chimere.e
</code>

En online mode \\
<code>
         echo "Running on Neurus-CNEA with hardcoded SRUN :"
         echo "PWD: "`pwd`
         touch chimwrf.conf
         rm chimwrf.conf
         echo "0-$((nproc_wrf-1))" ${wrf_exe} > chimwrf.conf
         echo "${nproc_wrf}-$((nproc_wrf+nproc_chimere-1))" ./chimere.e >> chimwrf.conf
         chmod +x chimwrf.conf
         echo "SRUN MPDM config file for NEURUS-CNEA chimwrf.conf in PWD:" `pwd`
         echo "----------------------------------------------------------"
         cat chimwrf.conf
         echo "----------------------------------------------------------"
         #SLURM MPMD (multiple program multiple data) MPI program
         #cmd_debug="-vv --slurmd-debug=debug"
         cmd="srun -A ${accountneurus} -p ${partitionneurus} -n $((nproc_wrf+nproc_chimere)) ${cmd_debug} -l --multi-prog chimwrf.conf"
         echo "CMD: " $cmd
         #take off
         unset I_MPI_PMI_LIBRARY.
         export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
         #cmd="echo hello word"
         ${cmd[@]}
         [ $? -eq 0 ] || { echo "Abnormal termination of chimwrf.conf"; exit 1; }
</code>

En runwrfonly \\
<code>
         time srun --verbose -A ${accountneurus} -p ${partitionneurus} ${mpiparams} -n ${nproc_wrf} ${wrf_exe}
        # time ${my_mpirun} ${mpiparams} -np ${nproc_wrf} ${wrf_exe} 
</code>
