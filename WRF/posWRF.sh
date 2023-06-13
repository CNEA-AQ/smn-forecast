#!/bin/bash

name=test
date="hoy"
clat=$1 #-39.704213;
clon=$2 #-61.833573;
DX=$3	       #Ancho dominio [km]
DY=$4	       #Alto  dominio [km]

#SMN:------------------
ref_lat=-35.00
ref_lon=-65.00
truelat1=-35.0
truelat2=-35.0
stand_lon=-65.0
dx=4000
dy=4000
#----------------------

clat=$ref_lat
clon=$ref_lon
#dx=25000.0; #11000.0; #   #Resolucion-X [m]
#dy=25000.0; #11000.0; #   #Resolución-Y [m]

srs_latlon="+proj=longlat +ellps=WGS84 +no_defs" 
srs_local="+proj=lcc lat_1=${truelat1} lat_2=${truelat2} lat_0=${ref_lat} lon_0=${ref_lon}"

echo $srs_local.
epsg_latlon=4326        #latlon wgs84
epsg_local=32721        #WGS-84 UTM Southern Hemisphere


#Mesh:

e_ew=$(bc -l <<< "scale=5;o=$DX*1000/$dx;scale=0;o/1.0") # == nx
e_sn=$(bc -l <<< "scale=5;o=$DY*1000/$dy;scale=0;o/1.0") # == ny

echo "e_ex: $e_ew e_sn; $e_sn; n_elementos= $( bc -l <<< "${e_ew}*${e_sn}" )"

#read xc yc zc <<<$(gdaltransform -s_srs epsg:${epsg_latlon} -t_srs epsg:${epsg_local} <<< $( echo "${clon} ${clat}" ) )
read xc yc zc <<<$(gdaltransform -s_srs "${srs_latlon}" -t_srs "${srs_local}" <<< $( echo "${clon} ${clat}" ) )

xini=$(bc -l <<<"$xc-$DX*1000*0.5")
yini=$(bc -l <<<"$yc-$DY*1000*0.5")
xfin=$(bc -l <<<"$xc+$DX*1000*0.5")
yfin=$(bc -l <<<"$yc+$DY*1000*0.5")


printf "id;wkt\n" > ${name}_${date}_${DX}x${DY}.csv
for ((j=0;j<${e_sn};j++))
	do
        for ((i=0;i<${e_ew};i++)) 
        do
		X=$(bc -l <<<"$xini+$i*$dx") 
		Y=$(bc -l <<<"$yini+$j*$dy") 

		echo "$X $Y" >> testXY.txt

	done;
done;

#cat testXY.txt | gdaltransform -s_srs epsg:${epsg_local} -t_srs epsg:${epsg_latlon}
cat testXY.txt | gdaltransform -s_srs "${srs_local}" -t_srs "${srs_latlon}" > testLL.txt

##Plot:
##gnuplot meshPlot.gnu
#gnuplot << EOF
#set title "GRILLA DE MODELADO.\n Experimento: ${name}, Fecha:${date}\n CLAT:${clat}   CLON:${clon}. \n Extensión: ${DX}x${DY}km, Resolución: ${dx}x${dy}m"
#unset key
##set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 800 
#set terminal svg enhanced background rgb 'white' size 480, 720
#set output '${name}_${date}_${DX}x${DY}km.svg'
#set xrange [-80:-50]
#set yrange [-57:-20]
#set format x "%D %E" geographic
#set format y "%D %N" geographic
##set cbrange [ * : * ] noreverse writeback
##set rrange [ * : * ] noreverse writeback
#set pointsize 0.3
#set style increment default
#set style data lines
#set yzeroaxis
##set grid 
#set grid xtics ytics mxtics mytics lc rgb "#555555", lc rgb "#bbbbbb"
#set mytics 10
#set mxtics 10
#plot 'data/world_50m.txt' with lines lc rgb "gray" , \\
#     'data/provincias.txt' using 2:3 w filledcu fc rgb "#f0f0f0" fs solid 0.5 border lc rgb "#303030",\\
#     '${name}_${date}_${DX}x${DY}.mesh' with points lt 1 pt 7 lc rgb "#ff0000", 
#pause 1 "la"
#EOF
#
#;;
#
#esac

#Para obtener provincias.txt, bajé shapefile de internet. con Qgis lo exporte a .csv (usando plugin MMQGIS) y luego lo procese:
#sed 's/"//g' temp-nodes.csv | awk 'NR>1 {if(prev != $1){printf("\n%s\n",$0); prev=$1} else{printf("%s\n",$0);prev=$1}}' >provincias.txt

