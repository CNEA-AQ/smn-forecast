#!/bin/bash

date=`date +'%Y%m%d'`
#for ff in `seq -f %03.0f 24 3 48  `; do
  inf=https://thredds.silam.fmi.fi/thredds/dodsC/dailysilam_glob06_v5_8/runs/daily_silam_glob06_v5_8_RUN_2024-06-06T00:00:00Z
  #inf=https://thredds.silam.fmi.fi/thredds/dodsC/silam_glob06_v5_8/files/SILAM-AQ-glob06_v5_8_2023111200_${ff}.nc4
  ouf=`basename $inf`
  echo ncks -v 'a,b,a_half,b_half,vmr_dust.*' $inf $ouf
  ncks -v 'a,b,a_half,b_half,vmr_dust.*' $inf $ouf
  #ncatted -O -a units,time,o,c,"hours since 2009-03-10 00:00:00.000 UTC" bc.nc bc_modif_date.nc
  #ncatted -O -a SIMULATION_START_DATE,global,o,c,"2009-03-15T00:00:00"  bc_modif_date.nc bc_modif_date.nc

#done | xargs -I{} -t -P8 sh -c '{}'
