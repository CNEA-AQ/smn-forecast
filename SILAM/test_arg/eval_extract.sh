#!/bin/bash
#
# Extract model values for evaluation.
#

set -e
set -u
#set -x

script=`readlink -f "$0"`
scriptdir=`dirname $script`

variablesMET=""
variables="cnc_O3_gas cnc_NO_gas cnc_NO2_gas cnc_SO2_gas cnc_PM2_5 cnc_PM10 cnc_CO_gas"; dssuff=""

cd $scriptdir
for dataset in $@; do
  path_to_ds=../output
  ctl_file=$path_to_ds/${dataset}_nc/0PM.nc.ctl
  #extract to joint datasets
  outd=${dataset}${dssuff}
  mkdir -p ../TimeVars/${outd}
  for v in $variables $variablesMET; do
      #echo true extracting $v db $db dataset $dataset 
      dbvar=`echo $v |sed -e s/PM2_5/PM25/ -e s/cnc_// -e s/_gas//  -e 's/\*.*$//'`
      ncvar=`echo $v |sed -e 's/cnc_NOX_gas/\(cnc_NO_gas*1.53+cnc_NO2_gas\)/' -e 's/cnc_OX_gas/\(cnc_NO2_gas*1.04+cnc_O3_gas\)/'`

      outfile=../TimeVars/${outd}/${outd}_${dbvar}.nc
      echo "python3 extract_ts2nc.py $ctl_file $outfile \\\"$ncvar\\\""
  done 
done |  srun --partition=test --nodes=1 --account=project_2004363 --ntasks-per-node=1 --cpus-per-task=32 --time 30:00 xargs -r -t -P 32 -I XXX sh -c "XXX"

