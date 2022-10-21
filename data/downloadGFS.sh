#! /bin/csh -f
#
# c-shell script to download selected files from rda.ucar.edu using Wget
# NOTE: if you want to run under a different shell, make sure you change
#       the 'set' commands according to your shell's syntax
# after you save the file, don't forget to make it executable
#   i.e. - "chmod 755 <name_of_script>"
#
# Experienced Wget Users: add additional command-line flags here
#   Use the -r (--recursive) option with care
#   Do NOT use the -b (--background) option - simultaneous file downloads
#       can cause your data access to be blocked
set opts = "-N"
#
# Replace xxxxxx with your rda.ucar.edu password on the next uncommented line
# IMPORTANT NOTE:  If your password uses a special character that has special
#                  meaning to csh, you should escape it with a backslash
#                  Example:  set passwd = "my\!password"

set usr_email = 'ramiroespada@cnea.gov.ar' #'mdiaz@cnea.gov.ar'
set passwd = 'kx405MUG' #'aguara'

set num_chars = `echo "$passwd" |awk '{print length($0)}'`
if ($num_chars == 0) then
  echo "You need to set your password before you can continue"
  echo "  see the documentation in the script"
  exit
endif
@ num = 1
set newpass = ""
while ($num <= $num_chars)
  set c = `echo "$passwd" |cut -b{$num}-{$num}`
  if ("$c" == "&") then
    set c = "%26";
  else
    if ("$c" == "?") then
      set c = "%3F"
    else
      if ("$c" == "=") then
        set c = "%3D"
      endif
    endif
  endif
  set newpass = "$newpass$c"
  @ num ++
end
set passwd = "$newpass"
#
set cert_opt = ""
# If you get a certificate verification error (version 1.10 or higher),
# uncomment the following line:
#set cert_opt = "--no-check-certificate"
#
if ("$passwd" == "xxxxxx") then
  echo "You need to set your password before you can continue - see the documentation in the script"
  exit
endif
 
# authenticate - NOTE: You should only execute this command ONE TIME.
# Executing this command for every data file you download may cause
# your download privileges to be suspended.

wget $cert_opt -O auth_status.rda.ucar.edu --save-cookies auth.rda.ucar.edu.$$ --post-data="email=$usr_email&passwd=$passwd&action=login" https://rda.ucar.edu/cgi-bin/login
#
# download the file(s)
# NOTE:  if you get 403 Forbidden errors when downloading the data files, check
#        the contents of the file 'auth_status.rda.ucar.edu'
		
#set fechas = ('2009-03-16' '2009-03-28')
#set horas = (00 06 12 18)
set fechas = ('2022-10-21')
set horas = (00)
set dias_buffer=1	#dias antes y despues
foreach fecha ($fechas) 
	foreach i (`seq -$dias_buffer $dias_buffer`)
		set dia=`date -d$fecha"+"$i" days" "+%Y-%m-%d"`

		set anio=`date -d$dia +%Y`
		set mes=`date -d$dia +%m`
		set dia=`date -d$dia +%d`
		foreach hora ($horas)
			set link='https://rda.ucar.edu/data/ds083.2/grib2/'$anio'/'$anio'.'$mes'/fnl_'$anio$mes$dia'_'$hora'_00.grib2'
			#echo $link
			wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ $link
		end
	end
end

#wget $cert_opt $opts --load-cookies auth.rda.ucar.edu.$$ https://rda.ucar.edu/data/ds083.2/grib2/2007/2007.12/fnl_20071206_18_00.grib2

#Clean:
rm auth.rda.ucar.edu.$$ auth_status.rda.ucar.edu

#Si se los quiere transformar a netcdf: cdo -f nc copy fnl_*.grib2 filetry.nc

