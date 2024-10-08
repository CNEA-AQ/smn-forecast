#!/bin/sh

for ff in `seq -f %03.0f 24 3 48  `; do
  inf=SILAM-AQ-glob06_v5_8_2023062900_${ff}.nc4
  a=$(echo "( ${ff}-24+18)" | bc )
  aa=$(printf %03.0f $a)
  ouf=out/SILAM-AQ-glob06_v5_8_2009031500_${aa}.nc4

  echo "$inf => $ouf"
  ncatted -O -a units,time,o,c,"seconds since 2009-03-14 18:00:00.000 UTC" $inf  $ouf
  ncatted -O -a SIMULATION_START_DATE,global,o,c,"2009-03-15T18:00:00"   $ouf  $ouf
  ncatted -O -a _CoordinateModelRunDate,global,o,c,"2009-03-15T18:00:00Z" $ouf  $ouf


done

#Lastly, merge all the files into one.
cdo mergetime out/SILAM-AQ* out/bcon.nc
