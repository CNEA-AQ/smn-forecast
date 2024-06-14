#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"

#GRIDDESC FILE:
#  COORD-NAME
#    COORDTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT
#  GRID-NAME
#    COORD-NAME, XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK

    wrfout_file=$1

    if [[ ! -f "$wrfout_file" ]]; then
        echo "Error: WRFOUT file '$wrfout_file' does not exist."
        return 1
    fi

    # Extract grid information from the WRFOUT file
    map_proj=$(ncdump -h "$wrfout_file" | grep 'MAP_PROJ ='                    | awk '{printf("%s"  , $3)}')

    nx=$(      ncdump -h "$wrfout_file" | grep 'WEST-EAST_GRID_DIMENSION ='    | sed 's/f//g' | awk '{print $3}')
    ny=$(      ncdump -h "$wrfout_file" | grep 'SOUTH-NORTH_GRID_DIMENSION ='  | sed 's/f//g' | awk '{print $3}')


    dx=$(      ncdump -h "$wrfout_file" | grep 'DX ='           | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    dy=$(      ncdump -h "$wrfout_file" | grep 'DY ='           | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    ref_lat=$( ncdump -h "$wrfout_file" | grep 'MOAD_CEN_LAT =' | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    ref_lon=$( ncdump -h "$wrfout_file" | grep 'STAND_LON ='    | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    truelat1=$(ncdump -h "$wrfout_file" | grep 'TRUELAT1 ='     | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    truelat2=$(ncdump -h "$wrfout_file" | grep 'TRUELAT2 ='     | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    clat=$(    ncdump -h "$wrfout_file" | grep 'CEN_LAT ='      | sed 's/f//g' | awk '{printf("%.3f", $3)}')
    clon=$(    ncdump -h "$wrfout_file" | grep 'CEN_LON ='      | sed 's/f//g' | awk '{printf("%.3f", $3)}')
  
#echo $dx,$dy, $ref_lat, $ref_lon, $truelat1 $truelat2 $clat $clon
  #map_proj =
  #1: Lambert Conformal
  #2: Polar Stereographic
  #3: Mercator
  #6: latitude and longitude (including global)

  #if (( $map_proj == 1 )) #'lambert' ) 
  #then
    coordtype=2
    alp=$truelat1; bet=$truelat2; gam=$ref_lon
    proj_str="+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=$truelat1 +lat_2=$truelat2 +lon_0=$ref_lon"
#echo $proj_str
  #elif (( $map_proj == 2 )) 
  #   coordtype=4
  #   proj_str=""
  #elif (( $map_proj == 3 ))  
  #   coordtype=3
  #   proj_str="+proj=merc +a=6370000.0 +b=6370000.0 +lon_0=${ref_lon} +lat_ts="
  #fi

  read xc yc ellipsoidh <<<$( gdaltransform -s_srs epsg:4326 -t_srs "${proj_str}" <<< $( echo "${clon} ${clat}" ) )
  xmin=$(bc -l <<<"$xc-$dx*$nx*1.0")
  ymin=$(bc -l <<<"$yc-$dy*$ny*0.5")
#echo "xc yc "$xc, $yc
#echo "xmin, ymin" $xmin, $ymin
  
  grid_name=TEST_GRID
  coord_name=TEST_COORD
  # Format and write the GRIDDESC file
 echo "'${coord_name}'                                  "
 echo "${coordtype} $alp $bet $gam $clon $clat        "
 echo "'${grid_name}'                                   "
 echo "'${coord_name}'  $xmin $ymin $dx $dy $nx $ny $nz "

printf "%-12s \n" ${coord_name}
printf "%-6d %.4f %9.6f %9.6f %9.6f %9.5f\n" ${coordtype} $alp $bet $gam $clon $clat        
printf "%-12s \n" ${grid_name}
printf "%.2f %.2f %.1f %.1f %5d %5d %5d\n" $xmin $ymin $dx $dy $nx $ny 1


