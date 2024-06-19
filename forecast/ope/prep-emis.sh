#!/bin/bash

today="2024-06-11" #`date +'%Y-%m-%d'`
yesterday=`date -d "${today} - 1 day" +'%Y-%m-%d'`
tomorrow=`date -d "${today}  + 1 day" +'%Y-%m-%d'`

#Prep. Emissions:
echo "* * * * * * * * * * * * "
echo "* Prep emission files   "

#============================================================
# Silam:
basedir=${HOME}/github/CNEA-AQ/smn-forecast/forecast
#wrf_dir=${basedir}/ope/wrf
wrf_dir=${HOME}/forecast/ope/wrf
emis_dir=${HOME}/data/silam/emis

if [ ! -d slm ]; then mkdir slm; mkdir slm/emis; fi
cd slm/emis

#----------------------
# "Static" inventories 
echo -e "Processing \e[35m static emission files\e[0m (CAMs v5.3, dust, soil-NO, sea salt)"
#Static emision files:
 ln -s ${emis_dir}/CAMS_GLOB_v53_ARG_CB5v2_as_2022 .   # Anthro (CAMS  v5.3)
 ln -s ${emis_dir}/dust-simple                     .   # dust emis             
 ln -s ${emis_dir}/soil_no                         .   # soil NO            
 ln -s ${emis_dir}/sslt                            .   # sea salt and DMS   
#ln -s ${emis_dia}/lightning                       .   # thunders           
#ln -s ${emis_dia}/aircraft                        .   # airplanes

#----------------------
## BVOCs  (MEGAN v3.3)
echo -e "Processing \e[32m BVOCs (w/MEGAN) v3.3\e[0m"
 cp -r ${emis_dir}/megan3.3    .   # sea salt and DMS   
 cd megan3.3
 #edito namelist:
 sed -i "s/start_date.*/start_date=\"${today} 00:00:00\"/" megan_namelist
 sed -i "s/end_date.*/end_date=\"${tomorrow} 00:00:00\"/" megan_namelist
 sed -i "s@met_files.*@met_files=\"${wrf_dir}/wrfout_d01_<date>_00:00:00\"@" megan_namelist
 #ejecuto:
 ./megan_v3.3.exe < megan_namelist                                        # prep-megan ya deberia estár corrido en emis_dir
 
cd ..
#----------------------
## Fires  (FINN  v2.5)
echo -e "Processing \e[31m Fires (w/FINN v2.5)\e[0m "
 cp -r ${emis_dir}/finn2.5  .   # sea salt and DMS   
 cd finn2.5
 #edito namelist
 sed -i "s/start_date=.*/start_date=\"${yesterday}\"/" fires_namelist     #Acá estoy haciendo "trampa" por que uso los fuegos de ayer!
 sed -i "s/end_date=.*/end_date=\"${yesterday}\"/"     fires_namelist     #todavia no hay disponible un pronóstico de fuegos..
 #ejecuto:
 ./finn2silam.exe < fires_namelist

 #newStart=`date +'seconds since %Y-%m-%d 00:00:00 UTC'`
 newStart=`date +'%Y-%m-%d'`
 ncatted -h -O -a units,time,m,c,"seconds since ${newStart} 00:00:00 UTC" emis_fires*
 ncatted -h -O -a _CoordinateModelRunDate,global,m,c,"${newStart} 00:00:00Z" emis_fires_*
 mv emis_fires_* emis_fires_${newStart}.nc

cd ..
#============================================================
# CMAQ
#
# newTitle=" OUTPUT FROM WRF V4    MODEL"
# ncatted -h -O -a title,global,m,c,"${newTitle}" ${wrf_dir}/wrfout_d01_*

#============================================================

echo "* * * * * * * * * * * * "
