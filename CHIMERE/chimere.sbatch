#!/bin/bash
#SBATCH -J CHIMERE               # Job NAME (human)
#SBATCH -p PQUIMICA              # Job QUEUE (human)
#SBATCH -A pquimica
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

cd $SLURM_SUBMIT_DIR

#echo $SLURM_JOB_NODELIST > $SLURM_JOBID.machines
#echo $SLURM_JOB_NODELIST | uniq > $SLURM_JOBID.mpd.hosts
#    NUM_PROCS=`cat $SLURM_JOBID.machines|wc -l`
#    NUM_NODES=`cat $SLURM_JOBID.mpd.hosts|wc -l`

srun --nodes=${SLURM_NNODES} bash -c 'hostname'> $SLURM_JOBID.machines
srun --nodes=${SLURM_NNODES} bash -c 'hostname' | sort -r | uniq > $SLURM_JOBID.mpd.hosts

#NUM_PROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
#NUM_NODES=$SLURM_NNODES

NUM_PROCS=`cat $SLURM_JOBID.machines|wc -l`
NUM_NODES=`cat $SLURM_JOBID.mpd.hosts|wc -l`

echo NUM_PROCS = $NUM_PROCS
echo NUM_NODES = $NUM_NODES

#mpdbootoption=" --file=$SLURM_SUBMIT_DIR/$SLURM_JOBID.mpd.hosts -n $NUM_NODES --rsh=ssh --chkup "
#mpdbootoption=" -n $NUM_NODES --rsh=ssh --chkup "
#mpioption=" -machinefile $SLURM_SUBMIT_DIR/$SLURM_JOBID.machines -np $NUM_PROCS"
#mpioption=" -bootstrap slurm  -np $NUM_PROCS "

#source ./mychimere/statcodes_paths.sh

firstdate=2017010100
lastdate=2017010100
incday=1
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

cd $SLURM_SUBMIT_DIR

#echo $SLURM_JOB_NODELIST > $SLURM_JOBID.machines
#echo $SLURM_JOB_NODELIST | uniq > $SLURM_JOBID.mpd.hosts
#    NUM_PROCS=`cat $SLURM_JOBID.machines|wc -l`
#    NUM_NODES=`cat $SLURM_JOBID.mpd.hosts|wc -l`

srun --nodes=${SLURM_NNODES} bash -c 'hostname'> $SLURM_JOBID.machines
srun --nodes=${SLURM_NNODES} bash -c 'hostname' | sort -r | uniq > $SLURM_JOBID.mpd.hosts

#NUM_PROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
#NUM_NODES=$SLURM_NNODES

NUM_PROCS=`cat $SLURM_JOBID.machines|wc -l`
NUM_NODES=`cat $SLURM_JOBID.mpd.hosts|wc -l`

echo NUM_PROCS = $NUM_PROCS
echo NUM_NODES = $NUM_NODES

#mpdbootoption=" --file=$SLURM_SUBMIT_DIR/$SLURM_JOBID.mpd.hosts -n $NUM_NODES --rsh=ssh --chkup "
#mpdbootoption=" -n $NUM_NODES --rsh=ssh --chkup "
#mpioption=" -machinefile $SLURM_SUBMIT_DIR/$SLURM_JOBID.machines -np $NUM_PROCS"
#mpioption=" -bootstrap slurm  -np $NUM_PROCS "

#source ./mychimere/statcodes_paths.sh
firstdate=2017010100
lastdate=2017010100
incday=1

parfile=chimere.par.test.offline
#----------------------------------------------------------------------
curdate=${firstdate}
nh=`expr ${incday} \* 24`
#----------------------------------------------------------------------
# prod or devel
typmod=prod
#----------------------------------------------------------------------
# simulation
# First block, no restart
chimrestart=no
#while [ $curdate -le $lastdate ] ; do
#    datestr=`date +"%Y%m%d%H%M%S"`
#    curlog=${logdir}/${logfile}.${datestr}.log
#    srun ./chimere.sh ${parfile} f ${curdate} ${nh} --${typmod} --restart ${chimrestart} 2>&1 | tee ${curlog} || exit 1
#    echo "Log file for this simulation: "${curlog}
#    chimrestart=yes
#    curdate=$(date -u -d "${curdate:0:8} ${curdate:8:2} ${nh} hour" +%Y%m%d%H)
#done
bin="chimere.sh"
#param="chimere.par.test.offline s 2017010100 24"
echo
echo
echo -n "pwd: "
pwd
echo -n "which mpiexec: "
which mpiexec
echo -n "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH
#echo mpdboot  $mpdbootoption 
#mpdboot $mpdbootoption  
#echo mpiexec $mpioption  $bin
echo -n "BEGIN MPI RUN date:"
date
echo
#mpiexec $mpioption  $bin
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

echo srun --verbose -n $SLURM_NTASKS $bin
#srun --verbose -n $SLURM_NTASKS $bin ${parfile} f ${curdate} ${nh} --${typmod} --restart ${chimrestart}
bash $bin ${parfile} f ${curdate} ${nh} --${typmod} --restart ${chimrestart}
echo
echo -n "END MPI RUN date:"
date
