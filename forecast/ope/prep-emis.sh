#!/bin/bash

today=`date +'%Y-%m-%d'`
yesterday=`date -d "yesterday" +'%Y-%m-%d'`
tomorrow=`date -d "tomorrow" +'%Y-%m-%d'`
#Prep. Emissions:
echo "* * * * * * * * * * * * "
echo "* Prep emission files   "

#============================================================
# Silam:
basedir=${HOME}/github/CNEA-AQ/smn-forecast/forecast
wrf_dir=${basedir}/ope/wrf
emis_dir=${HOME}/data/silam/emis

if [ ! -d slm ]; then mkdir slm; mkdir slm/emis; fi
cd slm/emis

#----------------------
# "Static" inventories 
echo -e "Processing \e[35m static emission files\e[0m (CAMs v5.3, dust, soil-NO, sea salt)"
#Static emision files:
 ln -s ${emis_dir}/CAMS_GLOB_v53_ARG_CB5v2_as_2022 .   # Anthro (CAMS  v5.3)
 ln -s ${emis_dir}/dust-simple                     .   # dust emis             
 ln -s ${emis_dir}/soil-NO                         .   # soil NO            
 ln -s ${emis_dir}/sslt                            .   # sea salt and DMS   
#ln -s ${emis_dia}/lightning                       .   # thunders           
#ln -s ${emis_dia}/aircraft                        .   # airplanes

#----------------------
## BVOCs  (MEGAN v3.3)
echo -e "Processing \e[32m BVOCs (w/MEGAN) v3.3\e[0m"
 cp -r ${emis_dir}/megan3.3    .   # sea salt and DMS   
 cd megan3.3

 sed -i "s/start_date.*/start_date=\"${today} 00:00:00\"/" megan_namelist
 sed -i "s/end_date.*/end_date=\"${tomorrow} 00:00:00\"/" megan_namelist
 sed -i "s@met_files.*@met_files=\"${wrf_dir}/wrfout_d01_<date>_00:00:00\"@" megan_namelist

 ./megan_v3.3.exe < megan_namelist
 if [ ! -f emis_bio_d01* ] ; then ./megan_v3.3.exe < megan_namelist; fi
 
cd ..
#----------------------
## Fires  (FINN  v2.5)
echo -e "Processing \e[31m Fires (w/FINN v2.5)\e[0m "
 cp -r ${emis_dir}/finn2.5  .   # sea salt and DMS   
 cd finn2.5
 
 sed -i "s/start_date=.*/start_date=\"${yesterday}\"/" fires_namelist
 sed -i "s/end_date=.*/end_date=\"${yesterday}\"/" fires_namelist
 ./finn2silam.exe < fires_namelist

 newStart=`date +'seconds since %Y-%m-%d 00:00:00 UTC'`
 ncatted -h -O -a units,time,m,c,"${newStart}" emis_fires*
cd ..
#============================================================
# CMAQ
#
# newTitle=" OUTPUT FROM WRF V4    MODEL"
# ncatted -h -O -a title,global,m,c,"${newTitle}" ${wrf_dir}/wrfout_d01_*


#============================================================

echo "* * * * * * * * * * * * "
