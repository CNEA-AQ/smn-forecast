#!/usr/bin/env csh
#
# c-shell script to download selected files from rda.ucar.edu using WGET

set fechas = ('2024-06-03')
set horas = (00 06 12 18)
set dias_buffer = 1           #dias antes y despues

set outDir = 'GFS/'

#wget options:
set opts = "-N -P "$outDir
#              (!) Do NOT use the -b (--background) option - simultaneous file downloads can cause your data access to be blocked

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

foreach fecha ($fechas)
   foreach i (`seq -$dias_buffer $dias_buffer`)
      set dia=`date -d$fecha"+"$i" days" "+%Y-%m-%d"`
      set anio=`date -d$dia +%Y`
      set mes=`date -d$dia +%m`
      set dia=`date -d$dia +%d`
      foreach hora ($horas)
        set file='fnl_'$anio$mes$dia'_'$hora'_00.grib2'
        echo $outDir$file
        if ( ! -f $outDir$file ) then
           set link='https://stratus.rda.ucar.edu/ds083.2/grib2/'$anio'/'$anio'.'$mes'/'$file
           wget $cert_opt $opts $link
        endif
      end
   end
end
