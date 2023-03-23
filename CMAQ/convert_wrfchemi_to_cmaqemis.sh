#!/bin/bash

## Set the input and output directories
#input_dir="wrfout_d01_01-0102"
#output_dir="wrfout_d01_01-0102_modifTitle"

files=($(ls wrfchemi*))
for inFile in "${files[@]}"; do

	file=emis_mole_$inFile.nc

cp $inFile $file
#file=wrfchemi_test 

echo $file

polluts=($(ncdump -h $file | grep "float E_" | sed 's/(/ /g' | awk '{printf("%s ",$2)}END{print ""}'))
dx=$(ncdump -h wrfchemi_test | grep "DX" | sed 's/^.* = //;s/;//;s/f//')
dy=$(ncdump -h wrfchemi_test | grep "DY" | sed 's/^.* = //;s/;//;s/f//')

for pollut in "${polluts[@]}";
do
	pollutNameCmaq=${pollut/E_/}	#es el nombre más probable, (luego vemos si hay que hacer modificaciones)

	echo $pollut
	#Cambiar nombre de pollut:

	#Cambiar unidad:	
	
	unit=$(ncdump -h $file | grep "${pollut}:" | grep "units" | sed 's/^.* = //;s/;//;s/\"//')
	

  	#if [[ $unit == "mol km^-2 hr^-1" ]]; then
  	if [[ $unit == *mol\ km* ]]; then
	echo $unit
		ncap2 -O -s "$pollut=$pollut*$dx*$dy*1e6" $file $file			#mol/km2.h -> moles/s
		ncatted -O -a units,$pollut,o,c,"moles/s         " $file

	elif [[ $unit == *ug\ m* ]]; then
	echo $unit
		ncap2 -O -s "$pollut=$pollut*$dx*$dy*1e6" $file $file			#ug/m2.s   -> g/s
		ncatted -O -a units,$pollut,o,c,"g/s             " $file
	fi;
done

#nombres de dimensiones:                    #netcdf wrfchemi_test2022-10-19_18\:00\:00 {        netcdf emis_mole_test.nc20221019 {
                                            #dimensions:                                        dimensions:
ncrename -d Time,TSTEP $file                #	Time = 24 ;                                     TSTEP = UNLIMITED ; // (24 currently)
ncrename -d DateStrLen,DATE-TIME,t $file    #	DateStrLen = UNLIMITED ; // (19 currently)      DATE-TIME = 2 ;
ncrename -d south_north,ROW $file 	    #	west_east = 199 ;                               COL = 197 ;
ncrename -d west_east,COL $file             #	south_north = 239 ;                             ROW = 237 ;
ncrename -d emissions_zdim_stag,LAY $file        # 	emissions_zdim = 1 ;                            LAY = 1 ;
                                            #                                                   VAR = 62 ;
#crear dimension VAR:                        
ncap2 -s "defdim(\"VAR\",${#polluts[@]})" -A $file

#modificar dimension TSTEP (remove DateStrLen y variable

ncks -O -x -v Times $file $file	 #saco la variable times para que no me traiga conflicto
#ncpdq -O -a DATE-TIME,2 $file	 #cambio dimension de DATE-TIME

#Creo variable TFLAG


#Metadata:

#        ncatted -O -h -a TITLE,global,m,c," OUTPUT FROM WRF V4    MODEL" "${input_dir}/${file}" "${output_dir}/${file}"

done


