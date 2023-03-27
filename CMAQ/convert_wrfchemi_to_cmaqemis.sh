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
dx=$(ncdump -h $file | grep "DX" | sed 's/^.* = //;s/;//;s/f//')
dy=$(ncdump -h $file | grep "DY" | sed 's/^.* = //;s/;//;s/f//')


#Ver como usar array/dicionario para mapear bien las especies:
#declare -A emis2cmaq=( ["BCI"]="" ["BCJ"]="" ["E_BENZENE"]="BENZ" ["E_BIGALK"]="" ["E_BIGENE"]="" ["E_C10H16"]="" ["E_C2H2"]="" ["E_C2H4"]="" ["E_C2H5OH"]="" ["E_C2H6"]="" ["E_C3H6"]="" ["E_C3H8"]="" ["E_CH2O"]="" ["E_CH3CHO"]="" ["E_CH3COCH3"]="" ["E_CH3COOH"]="" ["E_CH3OH"]="" ["E_CO"]="CO" ["E_HCOOH"]="" ["E_ISOP"]="ISOP" ["E_MEK"]="" ["E_NH3"]="NH3" ["E_NO"]="NO" ["E_NO2"]="NO2" ["E_OCI"]="" ["E_OCJ"]="" ["E_SO2"]="SO2" ["E_SULF"]="SULF" ["E_TOLUENE"]="TOL" ["E_XYLENES"]="") 

polluts=(${emis2cmaq[@]})

for pollut in "${polluts[@]}";
do
	pollutEmis=emis2cmaq[$pollut]
	#pollutNameCmaq=${pollut/E_/}	#es el nombre más probable, (luego vemos si hay que hacer modificaciones)

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
		ncatted -O -a var_desc,$pollut,o,c,"$pollutNameCmaq[1]   " $file	#agregar descripcion de variable
		#long_name=$(printf "%16s" $pollutNameCmaq
		long_name=$(printf "%-16s " $pollutNameCmaq )
		ncatted -O -a long_name,$pollut,o,c,"$long_name" $file	#agregar descripcion de variable

	ncrename -v $pollut,$pollutNameCmaq $file
done

#nombres de dimensiones:                    #netcdf wrfchemi_test2022-10-19_18\:00\:00 {        netcdf emis_mole_test.nc20221019 {
                                            #dimensions:                                        dimensions:
ncrename -d Time,TSTEP $file                #	Time = 24 ;                                     TSTEP = UNLIMITED ; // (24 currently)
#ncrename -d DateStrLen,DATE-TIME $file      #	DateStrLen = UNLIMITED ; // (19 currently)      DATE-TIME = 2 ;
ncrename -d south_north,ROW $file 	    #	west_east = 199 ;                               COL = 197 ;
ncrename -d west_east,COL $file             #	south_north = 239 ;                             ROW = 237 ;
ncrename -d emissions_zdim_stag,LAY $file   # 	emissions_zdim = 1 ;                            LAY = 1 ;
                                            #                                                   VAR = 62 ;
#crear dimension VAR:                        
ncap2 -s "defdim(\"VAR\",${#polluts[@]})" -A $file
ncap2 -s "defdim(\"DATETIME\",2)" -A $file

#modificar dimension TSTEP (remove DateStrLen y variable

dates=($(ncdump -v Times $file | sed -e '1,/data:/d;s/["\,\;]//g;' -e '$d' | awk 'NR>2{printf("%s",$0)}'))

#ncks -O -x -v Times $file $file       	 #saco la variable times para que no me traiga conflicto
#ncpdq -O -a DateStrLen,2 $file	 	#cambio dimension de DATE-TIME

#Creo variable TFLAG
ncap2 -A -s 'TFLAG = array(0,0,/$TSTEP,$VAR,$DATETIME/);' $file

for date in ${dates[@]}
do
	current_date=$(date -d "${date//_/ }" +"%Y-%m-%d %H")
	echo "Hora: "$current_date

        read YYYY MM DD HH <<< ${current_date//[-:\/_ ]/ }
        DDD=$(date -d "$YYYY/$MM/$DD" +"%j")             #dia juliano
        hh=$(echo $HH | awk '{print $0*1}');    #hora como integer
        i=$(( $hh % 12 ));                      #indice de variable time

	ncap2 -O -s "TFLAG($i,:,0) = ${YYYY}${DDD}; TFLAG($i,:,1) =  ${HH}0000;" $file -o $file
done

#Cambio nombre de dimension a DATE-TIME:
ncrename -d DATETIME,DATE-TIME $file    
 

#Metadata:
 #ncatted -O -h -a TITLE,global,m,c," OUTPUT FROM WRF V4    MODEL" ${file} ${file}

 ncatted -O -h -a IOAPI_VERSION,global,c,c,"ioapi-3.2: $Id: init3" $file $file
 ncatted -O -h -a EXEC_ID,global,c,c,"????????????????                                                                " $file $file
 ncatted -O -h -a FILEDESC,global,c,c,"Merged emissions output file from Mrggrid" $file $file  
 
 var_list=$(printf "%-16s " ${polluts[@]//E_/})
 ncatted -O -h -a VAR-LIST,global,c,c,"$var_list" $file $file
 
 #del metcro:
 ncatted -O -h -a FTYPE,global,c,f,1 $file $file
 ncatted -O -h -a CDATE,global,c,f,2023079 $file $file
 ncatted -O -h -a CTIME,global,c,f,141036 $file $file
 ncatted -O -h -a WDATE,global,c,f,2023079 $file $file
 ncatted -O -h -a WTIME,global,c,f,141036 $file $file
 ncatted -O -h -a SDATE,global,c,f,2019001 $file $file
 ncatted -O -h -a STIME,global,c,f,20000 $file $file
 ncatted -O -h -a TSTEP,global,c,f,10000 $file $file
 ncatted -O -h -a NTHIK,global,c,f,1 $file $file
 ncatted -O -h -a NCOLS,global,c,f,260 $file $file
 ncatted -O -h -a NROWS,global,c,f,337 $file $file
 ncatted -O -h -a NLAYS,global,c,f,32 $file $file
 ncatted -O -h -a NVARS,global,c,f,16 $file $file
 ncatted -O -h -a GDTYP,global,c,f,7 $file $file
 ncatted -O -h -a P_ALP,global,c,f,0. $file $file
 ncatted -O -h -a P_BET,global,c,f,0. $file $file
 ncatted -O -h -a P_GAM,global,c,f,-61. $file $file
 ncatted -O -h -a XCENT,global,c,f,-61. $file $file
 ncatted -O -h -a YCENT,global,c,f,0. $file $file
 ncatted -O -h -a XORIG,global,c,f,-3510000. $file $file
 ncatted -O -h -a YORIG,global,c,f,-7631603. $file $file
 ncatted -O -h -a XCELL,global,c,f,27000. $file $file
 ncatted -O -h -a YCELL,global,c,f,27000. $file $file
 ncatted -O -h -a VGTYP,global,c,f,-9999 $file $file
 ncatted -O -h -a VGTOP,global,c,f,5000. $file $file
 ncatted -O -h -a VGLVLS,global,c,f,1.,0.9938147  $file $file
 ncatted -O -h -a GDNAM,global,c,c,"PAPILAGRID" $file $file
 ncatted -O -h -a UPNAM,global,c,c,"OUTGM3IO      " $file $file


done


#polluts emis:
declare -A emis2cmaq=( ["BCI"]="" ["BCJ"]="" ["E_BENZENE"]="BENZ" ["E_BIGALK"]="" ["E_BIGENE"]="" ["E_C10H16"]="" ["E_C2H2"]="" ["E_C2H4"]="" ["E_C2H5OH"]="" ["E_C2H6"]="" ["E_C3H6"]="" ["E_C3H8"]="" ["E_CH2O"]="" ["E_CH3CHO"]="" ["E_CH3COCH3"]="" ["E_CH3COOH"]="" ["E_CH3OH"]="" ["E_CO"]="CO" ["E_HCOOH"]="" ["E_ISOP"]="ISOP" ["E_MEK"]="" ["E_NH3"]="NH3" ["E_NO"]="NO" ["E_NO2"]="NO2" ["E_OCI"]="" ["E_OCJ"]="" ["E_SO2"]="SO2" ["E_SULF"]="SULF" ["E_TOLUENE"]="TOL" ["E_XYLENES"]="") 
#polluts cmaq:
#declare -A cmaqPolluts=(
#["ACET"]            
#["ACROLEIN" ]       
#["ALD2" ]           
#["ALD2_PRIMARY"]   
#["ALDX"]          
#["BENZ"]           	E_BENZENE
#["BUTADIENE13"]     
#["CH4"]            
#["CH4_INV"]        
#["CL2"]             
#["CO"] 			E_CO             
#["CO2_INV"]         
#["ETH"]             
#["ETHA"]            
#["ETHY"]            
#["ETOH"]            
#["FORM"]            
#["FORM_PRIMARY"]    
#["HCL"]             
#["HONO"]            
#["IOLE"]            
#["ISOP"]            	E_ISOP
#["KET"]             
#["MEOH"]            
#["N2O_INV"]         
#["NAPH"]            
#["NH3"]             	E_NH3
#["NH3_FERT"]        
#["NO"]              	E_NO
#["NO2"]             	E_NO2
#["NVOL"]            
#["OLE"]             
#["PAL"]             
#["PAR"]             
#["PCA"]             
#["PCL"]             
#["PEC"]             
#["PFE"]             
#["PH2O"]            
#["PK"]              
#["PMC"]             
#["PMG"]             
#["PMN"]             
#["PMOTHR"]         
#["PNA"]             
#["PNCOM"]           
#["PNH4"]            
#["PNO3"]            
#["POC"]             
#["PRPA"]            
#["PSI"]             
#["PSO4"]            
#["PTI"]             
#["SO2"]             	E_SO2
#["SOAALK"]          
#["SULF"]            	E_SULF
#["TERP"]            
#["TOL"]             	E_TOLUENE
#["UNK"]             
#["UNR"]             
#["VOC_INV"]         
#["XYLMN"]
#)
