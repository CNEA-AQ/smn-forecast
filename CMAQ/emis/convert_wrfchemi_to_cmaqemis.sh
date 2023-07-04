#!/bin/bash

files=($(ls wrfchemi*))

start_date="2019-01-01"
  end_date="2019-02-01"

outfile=emis_mole_d01.nc
tmpfile=tmp_

for inFile in "${files[@]}"; do

	file=${tmpfile}emis_mole_$inFile.nc
	echo $file

	cp $inFile $file
	
	polluts=($(ncdump -h $file | grep "float E_" | sed 's/(/ /g' | awk '{printf("%s ",$2)}END{print ""}'))
	dx=$(ncdump -h $file | grep "DX" | sed 's/^.* = //;s/;//;s/f//')
	dy=$(ncdump -h $file | grep "DY" | sed 's/^.* = //;s/;//;s/f//')
	
	#saco variables innecesarias
	ncks -O -x -v XLAT $file $file       	#saco la variable para que no me traiga conflicto
	ncks -O -x -v XLONG $file $file       	#saco la variable para que no me traiga conflicto
	
	#Lista de variables CMAQ y asignación: (para desplegar en vim: q+q+f+(space)+a+(enter)+(Esc)+q  -> 60@q )
	declare -A var_all=(["ACET"]="" ["ACROLEIN"]="" ["ALD2"]="" ["ALD2_PRIMARY"]="" ["ALDX"]="" ["BENZ"]="E_BENZENE" ["BUTADIENE13"]="" ["CH4"]="" ["CH4_INV"]="" ["CL2"]="" ["CO"]="E_CO" ["CO2_INV"]="" ["ETH"]="" ["ETHA"]="" ["ETHY"]="" ["ETOH"]="" ["FORM"]="" ["FORM_PRIMARY"]="" ["HCL"]="" ["HONO"]="" ["IOLE"]="" ["ISOP"]="E_ISOP" ["KET"]="" ["MEOH"]="" ["N2O_INV"]="" ["NAPH"]="" ["NH3"]="E_NH3" ["NH3_FERT"]="" ["NO"]="E_NO" ["NO2"]="E_NO2" ["NVOL"]="" ["OLE"]="" ["PAL"]="" ["PAR"]="" ["PCA"]="" ["PCL"]="" ["PEC"]="" ["PFE"]="" ["PH2O"]="" ["PK"]="" ["PMC"]="" ["PMG"]="" ["PMN"]="" ["PMOTHR"]="" ["PNA"]="" ["PNCOM"]="" ["PNH4"]="" ["PNO3"]="" ["POC"]="" ["PRPA"]="" ["PSI"]="" ["PSO4"]="" ["PTI"]="" ["SO2"]="E_SO2" ["SOAALK"]="" ["SULF"]="E_SULF" ["TERP"]="" ["TOL"]="E_TOLUENE" ["UNK"]="" ["UNR"]="" ["VOC_INV"]="" ["XYLMN"]="")

	polluts=(${!var_all[@]})
	npolluts=0 #${#var_all[@]}
	
	for pollut in "${polluts[@]}";
	do
		if [[ ${var_all[$pollut]} == "" ]]
	       	then 
			#ncks -O -h -x -v $pollut $file $file       	#saco la variable times para que no me traiga conflicto
			#ncap2 -A -h -s 'swap = array(0.0,0,/$Time,$emissions_zdim_stag,$south_north,$west_east/);' $file
			#ncrename -h -v swap,$pollut $file			#Cambiar nombre de pollut:
			#falta unidades, long_name y var desc
			#unit="moles/s"
			continue; 
		else
			echo $pollut
			npolluts=$(($npolluts + 1 ))
			ncrename -h -v ${var_all[$pollut]},$pollut $file		#Cambiar nombre de pollut:
			#Cambiar unidad:	
			unit=$(ncdump -h $file | grep "${pollut}:" | grep "units" | sed 's/^.* = //;s/;//;s/\"//')
	
	  		if [[ $unit == *mol\ km* ]]; then
				ncap2 -O -h -s "$pollut=float($pollut*$dx*$dy*1e-6)" $file $file	#mol/km2.h -> moles/s
				unit="moles/s"
	
			elif [[ $unit == *ug\ m* ]]; then
				ncap2 -O -h -s "$pollut=float($pollut*$dx*$dy*1e-6)" $file $file	#ug/m2.s   -> g/s
				unit="g/s"
			fi;

			ncatted -O -h -a ".,$pollut,d,,;" $file $file	#borro atributos preexistentes

			units=$(printf "%-16s" $unit )
			var_desc=$(printf "%-80s" $pollut[1] )
			long_name=$(printf "%-16s" $pollut )
			ncatted -O -h -a long_name,$pollut,o,c,"$long_name" $file	#agregar atributo de variable
			ncatted -O -h -a units,$pollut,o,c,"$units" $file		#agregar atributo de variable
			ncatted -O -h -a var_desc,$pollut,o,c,"$var_desc" $file	#agregar atributo de variable
		fi;
	done
	#Elimino variables del wrfchemi que no tienen asignación                                                                   	
	ncks -O -x -v E_. $file $file 	#saco la variable times para que no me traiga conflicto	#nombres de dimensiones:                   

	#dimensiones:
	ncrename -h -d Time,TSTEP $file                #	Time           -> TSTEP
	#ncrename -d -h  DateStrLen,DATE-TIME $file     #	DateStrLen     -> DATE-TIME
	ncrename -h -d south_north,ROW $file 	    #	west_east      -> COL
	ncrename -h -d west_east,COL $file             #	south_north    -> ROW
	ncrename -h -d emissions_zdim_stag,LAY $file   # 	emissions_zdim => LAY
	
	# crear dimension VAR:                        
	ncap2 -s "defdim(\"VAR\",${npolluts})" -A $file
	ncap2 -s "defdim(\"DATETIME\",2)" -A $file
	
	#Creo variable TFLAG
	ncap2 -A -s 'TFLAG = array(0,0,/$TSTEP,$VAR,$DATETIME/);' $file
	
	dates=($(ncdump -v Times $file | sed -e '1,/data:/d;s/["\,\;]//g;' -e '$d' | awk 'NR>2{printf("%s",$0)}'))
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
	
		ncap2 -O -h -s "TFLAG($i,:,0) = ${YYYY}${DDD}; TFLAG($i,:,1) = ${HH}0000;" $file -o $file
	done
	#metadata de TFLAG:
	ncatted -O -h -a units,TFLAG,o,c,"<YYYYDDD,HHMMSS>" $file	
        ncatted -O -h -a long_name,TFLAG,o,c,"TFLAG           " $file
        ncatted -O -h -a var_desc,TFLAG,o,c,"Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " $file

	ncrename -d DATETIME,DATE-TIME $file    #Cambio nombre de dimension a DATE-TIME:

	ncks -O -h -x -v Times $file $file       	#saco la variable times para que no me traiga conflicto
done

#Merge todos en un solo archivo.
rm $outfile		#por las dudas de que ya exista lo borro.
ncrcat -4 ${tmpfile}* $outfile

#Filtro fechas
ncks -h -O -d TSTEP,1,23 $outfile $outfile
ncks -h -O -d COL,1,-2 $outfile $outfile
ncks -h -O -d ROW,1,-2 $outfile $outfile

#METADATA:
#Primero borro toda la metadata
ncatted -O -h -a ".,global,d,,;" $outfile

ncatted -O -h -a TITLE,global,m,c,"EMISSION FILE                " $outfile
ncatted -O -h -a IOAPI_VERSION,global,c,c,"ioapi-3.2: $Id: init3" $outfile
ncatted -O -h -a EXEC_ID,global,c,c,"????????????????           " $outfile

#Agrego a metadata info sacada del METCRO*:
ncatted -O -h -a FTYPE,global,c,i,1 $outfile
ncatted -O -h -a CDATE,global,c,i,2023001 $outfile
ncatted -O -h -a CTIME,global,c,i,000000 $outfile
ncatted -O -h -a WDATE,global,c,i,2023001 $outfile
ncatted -O -h -a WTIME,global,c,i,000000 $outfile
ncatted -O -h -a SDATE,global,c,i,2019001 $outfile
ncatted -O -h -a STIME,global,c,i,010000 $outfile
ncatted -O -h -a TSTEP,global,c,i,10000 $outfile
ncatted -O -h -a NTHIK,global,c,i,1 $outfile
ncatted -O -h -a NCOLS,global,c,i,260 $outfile
ncatted -O -h -a NROWS,global,c,i,337 $outfile
ncatted -O -h -a NLAYS,global,c,i,1 $outfile
ncatted -O -h -a NVARS,global,c,i,$npolluts $outfile
ncatted -O -h -a GDTYP,global,c,i,7 $outfile
ncatted -O -h -a P_ALP,global,c,f,0. $outfile
ncatted -O -h -a P_BET,global,c,f,0. $outfile
ncatted -O -h -a P_GAM,global,c,f,-61. $outfile
ncatted -O -h -a XCENT,global,c,f,-61. $outfile
ncatted -O -h -a YCENT,global,c,f,0. $outfile
ncatted -O -h -a XORIG,global,c,f,-3510000. $outfile
ncatted -O -h -a YORIG,global,c,f,-7631603. $outfile
ncatted -O -h -a XCELL,global,c,f,27000. $outfile
ncatted -O -h -a YCELL,global,c,f,27000. $outfile
ncatted -O -h -a VGTYP,global,c,f,-9999 $outfile
ncatted -O -h -a VGTOP,global,c,f,5000. $outfile
ncatted -O -h -a VGLVLS,global,c,f,1.,0.9938147  $outfile
ncatted -O -h -a GDNAM,global,c,c,"PAPILAGRID_CROSS" $outfile
ncatted -O -h -a UPNAM,global,c,c,"OUTGM3IO      " $outfile

variables=($(ncdump -h $file | grep "float" | sed 's/^.*float \(.*\)(.*/\1 /' | tr '\n' ' '))
IFS=$'\n' sorted_vars=($(sort <<<"${variables[*]}")); unset IFS;
var_list=$(printf "%-16s" ${sorted_vars[@]})
ncatted -O -h -a VAR-LIST,global,c,c,"$var_list" $outfile

ncatted -O -h -a FILEDESC,global,c,c,"Merged emissions output file from Mrggrid                                       " $outfile  
ncatted -O -h -a HISTORY,global,c,c,"" $outfile


#Elimino archivos temporales
rm ${tmpfile}*
