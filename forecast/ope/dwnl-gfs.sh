#!/usr/bin/env csh
#
# c-shell script to download selected files from rda.ucar.edu using WGET
set dia=`date +'%Y-%m-%d'`
set anio=`date -d$dia +%Y`
set mes=`date -d$dia +%m`
set dia=`date -d$dia +%d`

set CC = (00)
#set ForecastHours = (000 0001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024)
set ForecastHours = (000 003 006 009 012 015 018 021 024)
set outDir = 'GFS/'

#wget options:
set opts = "-N -P "$outDir #(!) Do NOT use the -b (--background) option - simultaneous file downloads can cause your data access to be blocked

# Check wget version.  Set the --no-check-certificate option 
# if wget version is 1.10 or higher
set v = `wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
set a = `echo $v | cut -d '.' -f 1`
set b = `echo $v | cut -d '.' -f 2`
if(100 * $a + $b > 109) then
  set cert_opt = "--no-check-certificate"
else
  set cert_opt = ""
endif

foreach FFF ($ForecastHours)
  set file='gfs.t'$CC'z.pgrb2.0p25.f'$FFF
  echo $outDir$file
  if ( ! -f $outDir$file ) then
     set link='https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.'$anio$mes$dia'/'$CC'/atmos/'$file
     wget $cert_opt $opts $link
  endif
end
   
