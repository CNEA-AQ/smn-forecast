#!/bin/bash
#
# Runs silam model 
#
exe=silam_v5_8pub.mahti_mpi 
template=arg-cb5-control #

export date=$(date +"%Y%m%d")
export ntasks=8 	     #this should be == x_divisions *  y_divisions
partition=medium

nohup srun --nodes=1 --tasks-per-node=$ntasks --account=project_2004363 --cpus-per-task=32 --threads-per-core=2 -J silam_arg --partition=${partition} --time=23:00:00 ${exe} ${template} 2>&1> out${date}.log &
