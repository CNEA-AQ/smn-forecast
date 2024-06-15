#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"

#GRIDDESC FILE:-----------------------------------------------------
#  COORD-NAME                                                      |   
#    COORDTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT                  |  
#  GRID-NAME                                                       |  
#    COORD-NAME, XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK   |
#-------------------------------------------------------------------
    wrfout_file=$1

    if [[ ! -f "$wrfout_file" ]]; then
        echo "Error: WRFOUT file '$wrfout_file' does not exist."
        return 1
    fi

    # Extract grid information from the WRFOUT file
    ncdump -h "$wrfout_file"  > tmp_ncdump_-h
    map_proj=$(tail -n 200 tmp_ncdump_-h | grep 'MAP_PROJ ='                    | awk '{printf("%s"  , $3)}')
    nx=$(      tail -n 200 tmp_ncdump_-h | grep 'WEST-EAST_GRID_DIMENSION ='    | awk '{print $3}')
    ny=$(      tail -n 200 tmp_ncdump_-h | grep 'SOUTH-NORTH_GRID_DIMENSION ='  | awk '{print $3}')
    dx=$(      tail -n 200 tmp_ncdump_-h | grep 'DX ='           | sed 's/f//g' | awk '{printf("%8.3f",$3)}')
    dy=$(      tail -n 200 tmp_ncdump_-h | grep 'DY ='           | sed 's/f//g' | awk '{printf("%8.3f",$3)}')
    clat=$(    tail -n 200 tmp_ncdump_-h | grep 'CEN_LAT ='      | sed 's/f//g' | awk '{printf("%8.5f",$3)}')
    clon=$(    tail -n 200 tmp_ncdump_-h | grep 'CEN_LON ='      | sed 's/f//g' | awk '{printf("%8.5f",$3)}')
     ref_lat=$(tail -n 200 tmp_ncdump_-h | grep 'MOAD_CEN_LAT =' | sed 's/f//g' | awk '{printf("%8.3f",$3)}')
     ref_lon=$(tail -n 200 tmp_ncdump_-h | grep 'STAND_LON ='    | sed 's/f//g' | awk '{printf("%8.3f",$3)}')
    truelat1=$(tail -n 200 tmp_ncdump_-h | grep 'TRUELAT1 ='     | sed 's/f//g' | awk '{printf("%8.5f",$3)}')
    truelat2=$(tail -n 200 tmp_ncdump_-h | grep 'TRUELAT2 ='     | sed 's/f//g' | awk '{printf("%8.5f",$3)}')
    rm tmp_ncdump_-h

  #map_proj =
  #1: Lambert Conformal
  #2: Polar Stereographic
  #3: Mercator
  #6: latitude and longitude (including global)

  #if (( $map_proj == 1 )) #'lambert' ) 
  #then
    coordtype=2
    alp=${truelat1}
    bet=${truelat2}
    gam=${ref_lon}
    proj_str="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lat_0=${ref_lat} +lon_0=${ref_lon} +unit=m" #+a=6370000.0 +b=6370000.0 
    proj_str=${proj_str//= /=}
  #elif (( $map_proj == 2 )) 
  #   coordtype=4
  #   proj_str=""
  #elif (( $map_proj == 3 ))  
  #   coordtype=3
  #   proj_str="+proj=merc +a=6370000.0 +b=6370000.0 +lon_0=${ref_lon} +lat_ts=${ref_lat}"
  #fi
echo "proj:" $proj_str
  read xc yc ellipsoidh <<<$( gdaltransform -s_srs epsg:4326 -t_srs "${proj_str}" <<< $( echo "${clon} ${clat}" ) )
  xmin=$(bc -l <<<"$xc-$dx*$nx*0.5")
  ymin=$(bc -l <<<"$yc-$dy*$ny*0.5")
  echo "xc yc "$xc, $yc
  echo "xmin, ymin" $xmin, $ymin
  
  grid_name=TEST_GRID
  coord_name=TEST_COORD
  # Format and write the GRIDDESC file
#@ echo "'${coord_name}'                                  "
#@ echo "${coordtype} $alp $bet $gam $clon $clat        "
#@ echo "'${grid_name}'                                   "
#@ echo "'${coord_name}' $xmin $ymin $dx $dy $nx $ny $nz "

printf "%-12s \n" ${coord_name}
printf "%-8d %9.5f %9.5f %9.5f %9.5f %10.7f \n" ${coordtype} ${alp} ${bet} ${gam} ${clon} ${clat};
printf "%-12s \n" ${grid_name}
printf "%-12s %.1f %.1f %.2f %.2f %5d %5d %5d\n" ${coord_name} $xmin $ymin $dx $dy $nx $ny 1


