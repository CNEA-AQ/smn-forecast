#!/bin/bash

files=($(ls wrfchemi*))

for inFile in "${files[@]}"; do

	file=emis_mole_$inFile.nc
	echo $file

	cp $inFile $file
	
	polluts=($(ncdump -h $file | grep "float E_" | sed 's/(/ /g' | awk '{printf("%s ",$2)}END{print ""}'))
	dx=$(ncdump -h $file | grep "DX" | sed 's/^.* = //;s/;//;s/f//')
	dy=$(ncdump -h $file | grep "DY" | sed 's/^.* = //;s/;//;s/f//')
	
	#saco variables innecesarias
	ncks -O -x -v XLAT $file $file       	#saco la variable para que no me traiga conflicto
	ncks -O -x -v XLONG $file $file       	#saco la variable para que no me traiga conflicto
	
	#Lista de variables y asignación:
	declare -A emis2cmaq=(["ACET"]="" ["ACROLEIN"]="" ["ALD2"]="" ["ALD2_PRIMARY"]="" ["ALDX"]="" ["BENZ"]="E_BENZENE" ["BUTADIENE13"]="" ["CH4"]="" ["CH4_INV"]="" ["CL2"]="" ["CO"]="E_CO" ["CO2_INV"]="" ["ETH"]="" ["ETHA"]="" ["ETHY"]="" ["ETOH"]="" ["FORM"]="" ["FORM_PRIMARY"]="" ["HCL"]="" ["HONO"]="" ["IOLE"]="" ["ISOP"]="E_ISOP" ["KET"]="" ["MEOH"]="" ["N2O_INV"]="" ["NAPH"]="" ["NH3"]="E_NH3" ["NH3_FERT"]="" ["NO"]="E_NO" ["NO2"]="E_NO2" ["NVOL"]="" ["OLE"]="" ["PAL"]="" ["PAR"]="" ["PCA"]="" ["PCL"]="" ["PEC"]="" ["PFE"]="" ["PH2O"]="" ["PK"]="" ["PMC"]="" ["PMG"]="" ["PMN"]="" ["PMOTHR"]="" ["PNA"]="" ["PNCOM"]="" ["PNH4"]="" ["PNO3"]="" ["POC"]="" ["PRPA"]="" ["PSI"]="" ["PSO4"]="" ["PTI"]="" ["SO2"]="E_SO2" ["SOAALK"]="" ["SULF"]="E_SULF" ["TERP"]="" ["TOL"]="E_TOLUENE" ["UNK"]="" ["UNR"]="" ["VOC_INV"]="" ["XYLMN"]="")

	polluts=(${!emis2cmaq[@]})
	npolluts=${#emis2cmaq[@]}
	
	for pollut in "${polluts[@]}";
	do
		echo $pollut
		if [[ ${emis2cmaq[$pollut]} == "" ]]
	       	then 
			#ncks -O -x -v $pollut $file $file       	#saco la variable times para que no me traiga conflicto
			#ncap2 -A -s 'swap = array(0,0,/$Time,$emissions_zdim_stag,$south_north,$west_east/);' $file
			#ncrename -v swap,$pollut $file			#Cambiar nombre de pollut:
			#falta unidades, long_name y var desc
			continue; 
		else
			ncrename -v ${emis2cmaq[$pollut]},$pollut $file		#Cambiar nombre de pollut:
			#Cambiar unidad:	
			unit=$(ncdump -h $file | grep "${pollut}:" | grep "units" | sed 's/^.* = //;s/;//;s/\"//')
	
	  		if [[ $unit == *mol\ km* ]]; then
				ncap2 -O -s "$pollut=float($pollut*$dx*$dy*1e6)" $file $file	#mol/km2.h -> moles/s
				unit="moles/s"
	
			elif [[ $unit == *ug\ m* ]]; then
				ncap2 -O -s "$pollut=float($pollut*$dx*$dy*1e6)" $file $file	#ug/m2.s   -> g/s
				unit="g/s"
			fi;

			ncatted -O -h -a ".,$pollut,d,,;" $file $file

			units=$(printf "%-16s" $unit )
			var_desc=$(printf "%-80s" $pollut[1] )
			long_name=$(printf "%-16s" $pollut )
			ncatted -O -a long_name,$pollut,o,c,"$long_name" $file	#agregar atributo de variable
			ncatted -O -a var_desc,$pollut,o,c,"$var_desc" $file	#agregar atributo de variable
			ncatted -O -a units,$pollut,o,c,"$units" $file		#agregar atributo de variable
			
		fi;
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
	ncap2 -s "defdim(\"VAR\",${npolluts})" -A $file
	ncap2 -s "defdim(\"DATETIME\",2)" -A $file
	
	#modificar dimension TSTEP (remove DateStrLen y variable
	
	dates=($(ncdump -v Times $file | sed -e '1,/data:/d;s/["\,\;]//g;' -e '$d' | awk 'NR>2{printf("%s",$0)}'))
	
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
	
		#Aca tuve que forzar por que el emis de entrada es de otra fecha:
		YYYY=2019;DDD=001
	
		ncap2 -O -s "TFLAG($i,:,0) = ${YYYY}${DDD}; TFLAG($i,:,1) = ${HH}0000;" $file -o $file
	done
	#metadata de TFLAG:
	ncatted -O -a units,TFLAG,o,c,"<YYYYDDD,HHMMSS>" $file	
        ncatted -O -a long_name,TFLAG,o,c,"TFLAG           " $file
        ncatted -O -a var_desc,TFLAG,o,c,"Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " $file

	ncrename -d DATETIME,DATE-TIME $file    #Cambio nombre de dimension a DATE-TIME:
	
	ncks -O -x -v Times $file $file       	#saco la variable times para que no me traiga conflicto

	#METADATA:
	#Primero borro toda la metadata
	 ncatted -O -h -a ".,global,d,,;" $file $file
	 
	 ncatted -O -h -a TITLE,global,m,c,"EMISSION FILE               " ${file} ${file}
        
         ncatted -O -h -a IOAPI_VERSION,global,c,c,"\$Id: @(#) ioapi library version 3.1 $                                           " $file $file   #ioapi-3.2: $Id: init3"
         ncatted -O -h -a EXEC_ID,global,c,c,"????????????????                                                                " $file $file

 	#Agrego a metadata info sacada del METCRO*:
	 ncatted -O -h -a FTYPE,global,c,i,1 $file $file
	 ncatted -O -h -a CDATE,global,c,i,2023079 $file $file
	 ncatted -O -h -a CTIME,global,c,i,141036 $file $file
	 ncatted -O -h -a WDATE,global,c,i,2023079 $file $file
	 ncatted -O -h -a WTIME,global,c,i,141036 $file $file
	 ncatted -O -h -a SDATE,global,c,i,2019001 $file $file
	 ncatted -O -h -a STIME,global,c,i,010000 $file $file
	 ncatted -O -h -a TSTEP,global,c,i,10000 $file $file
	 ncatted -O -h -a NTHIK,global,c,i,1 $file $file
	 ncatted -O -h -a NCOLS,global,c,i,260 $file $file
	 ncatted -O -h -a NROWS,global,c,i,337 $file $file
	 ncatted -O -h -a NLAYS,global,c,i,1 $file $file
	 ncatted -O -h -a NVARS,global,c,i,$npolluts $file $file
	 ncatted -O -h -a GDTYP,global,c,i,7 $file $file
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
        
	 IFS=$'\n' sorted_vars=($(sort <<<"${polluts[*]}")); unset IFS;
	 var_list=$(printf "%-16s" ${sorted_vars[@]})
	 ncatted -O -h -a VAR-LIST,global,c,c,"$var_list" $file $file
	 
	 ncatted -O -h -a FILEDESC,global,c,c,"Merged emissions output file from Mrggrid                                       " $file $file  
	 ncatted -O -h -a HISTORY,global,c,c,"" $file $file
	 
done

#Merge todos en un solo archivo.
rm out_emis_mole_d01.nc

#Merge todos los archivos en uno
ncrcat -4 emis_mole* out_emis_mole_d01.nc

#Elimino archivos individuales
rm emis_mole*

#Filtro fechas
ncks -h -O -d TSTEP,1,23 out_emis_mole_d01.nc out_emis_mole_d01.nc
ncks -h -O -d COL,1,-2 out_emis_mole_d01.nc out_emis_mole_d01.nc
ncks -h -O -d ROW,1,-2 out_emis_mole_d01.nc out_emis_mole_d01.nc





#declare -A emis2cmaq=(["BCI"]="" ["BCJ"]="" ["E_BENZENE"]="" ["E_BIGALK"]="" ["E_BIGENE"]="" ["E_C10H16"]="" ["E_C2H2"]="" ["E_C2H4"]="" ["E_C2H5OH"]="" ["E_C2H6"]="" ["E_C3H6"]="" ["E_C3H8"]="" ["E_CH2O"]="" ["E_CH3CHO"]="" ["E_CH3COCH3"]="" ["E_CH3COOH"]="" ["E_CH3OH"]="" ["E_CO"]="CO" ["E_HCOOH"]="" ["E_ISOP"]="ISOP" ["E_MEK"]="" ["E_NH3"]="NH3" ["E_NO"]="NO" ["E_NO2"]="NO2" ["E_OCI"]="" ["E_OCJ"]="" ["E_SO2"]="SO2" ["E_SULF"]="SULF" ["E_TOLUENE"]="TOL" ["E_XYLENES"]="") 


