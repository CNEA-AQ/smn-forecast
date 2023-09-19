#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"
   srsInp="epsg:4326"	#Los archivos vienen en latlon.
#Input files:
 GRIDDESC="GRIDDESC"    #"GRIDDESC_USA" #
 GRIDNAME="papila_grid" #"12US1"        #

coastline=coastline.gpkg	        #polygons to make a landmask
 landpoly=input/shp/land_10m_global.gpkg	#polygons to make a landmask
chlo_nc=input/CHLO_monthly_global.nc     #chlorophyl conc. files (from MODIS (MY1DMM)
 dms_nc=input/DMS_monthly_global.nc      #dms conc. files (DMS climatology of Lana et al. (2011))
#-------------------------------------------------
#(0) Get grid & proj parameters from GRIDDESC:
read projName xorig yorig dx dy nx ny nz <<< $( sed -n "/${GRIDNAME}/{n;p;q}" "$GRIDDESC" )
read COORDTYPE P_ALP P_BET P_GAM XCENT YCENT <<< $( sed -n "/${projName}/{n;p;q}" "$GRIDDESC" )
#read COORDTYPE truelat1 truelat2 stand_lon ref_lon ref_lat     #(en el namelist de wrf)
#truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
  if [ $COORDTYPE == 1 ]; then           #Geographic:
   srsOut="+proj=latlong +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 2 ]; then     #Lambert Conformal Conic:
   srsOut="+proj=lcc +lat_1=$P_ALP +lat_2=$P_BET +lon_0=$P_GAM +lat_0=$YCENT +a=6370000.0 +b=6370000.0 +units=m"
elif [ $COORDTYPE == 3 ]; then     #General Mercator
   srsOut="+proj=merc +lat_ts=$P_ALP +lon_0=$P_GAM +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 4 ]; then     #General tangent Stereografic
   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
elif [ $COORDTYPE == 5 ]; then     #UTM
   echo  "proyección: 5 (Universal Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 6 ]; then     #Polar Secant Stereographic
   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
elif [ $COORDTYPE == 7 ]; then     #Equatorial Mercator
   srsOut="+proj=merc +lat_ts=$P_ALP +lon_0=$P_GAM +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 8 ]; then     #Transverse Mercator
   echo  "proyección: 8 (Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 9 ]; then     #Lambert Azimuthal Equal-Area
   echo  "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada en esta aplicación."; stop
else
   echo  "codigo de proyección invalido. COORDTYPE"; stop
fi;

echo "SRS of Input  Grids: $srsInp"
echo "SRS of Output Grids: $srsOut"
#------------------------------------------------
#Parametros intermedios:
xmin=$(bc -l <<<"$xorig")
xmax=$(bc -l <<<"$xmin+$dx*$nx ")
ymin=$(bc -l <<<"$yorig")
ymax=$(bc -l <<<"$ymin+$dy*$ny ")

#bbox latlon: para hacer clippear el archivo de fuegos al dominio:
read lonmin latmin ellipsoidh <<<$( gdaltransform  -s_srs "${srsOut}" -t_srs "${srsInp}" <<< $( echo "${xmin} ${ymin}" ) )
read lonmax latmax ellipsoidh <<<$( gdaltransform  -s_srs "${srsOut}" -t_srs "${srsInp}" <<< $( echo "${xmax} ${ymax}" ) )
#echo "lon lat min,  $lonmin $latmin lon lat max: $lonmax $latmax"
#------------------------------------------------
# SURF:
#(1) Rasterizo las lineas de costas "muy fino" (100m o 200m), asigno 1/0 si toca o no costa. Luego regrid a grilla de modelado calculando "average".
dx_small=300; dy_small=300
echo "Calculating SURF.."
if [ ! -f surf_${GRIDNAME}.tif ] 
then
    
    if [ ! -f ${coastline} ]
    then
    	ogr2ogr -f "GPKG" ${coastline} $landpoly -nln coastline -nlt MULTILINESTRING -t_srs "${srsOut}" -clipsrc "${lonmin}" "${latmin}" "${lonmax}" "${latmax}"
    fi
    gdal_rasterize -q -l coastline -at -burn 1 -a_srs "${srsOut}" -tr ${dx_small} ${dy_small} -te $xmin $ymin $xmax $ymax -of GTiff ${coastline} tmp.tif
    gdalwarp -r average -t_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax tmp.tif surf_${GRIDNAME}.tif
    gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" surf_${GRIDNAME}.tif surf_${GRIDNAME}.nc 
    rm tmp.tif
else
    echo "surf.tif already exists!"
fi
# LAND: (same procedure than with SURF)
echo "Calculating LAND.."
if [ ! -f land_${GRIDNAME}.tif ] 
then
    ogr2ogr -f "GPKG" land.gpkg $landpoly -nln land -nlt MULTIPOLYGON -t_srs "${srsOut}" -clipsrc "${lonmin}" "${latmin}" "${lonmax}" "${latmax}"
    gdal_rasterize -q -l land -at -burn 1 -a_srs "${srsOut}" -tr ${dx_small} ${dy_small} -te $xmin $ymin $xmax $ymax -of GTiff land.gpkg tmp.tif
    gdalwarp -r average -t_srs "${srsOut}" -tr ${dx} ${dy} -te $xmin $ymin $xmax $ymax tmp.tif land_${GRIDNAME}.tif
    rm tmp.tif
else
   echo "land.tif already exists!"
fi
# OPEN = 1 - LAND - SURF
echo "Calculating OPEN.."
if [ ! -f open_${GRIDNAME}.tif ] 
then
    gdal_calc.py -S surf_${GRIDNAME}.tif --S_band=1 -L land_${GRIDNAME}.tif --L_band=1 --calc="1-L-S" --outfile=open_${GRIDNAME}.tif
    gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" open_${GRIDNAME}.tif open_${GRIDNAME}.nc 
else
    echo "open.tif already exists!"
fi
#-------------------------------------------------
#(1) Re-gridding: Agarrar los inputs y regrillarlos (e interpolar) segun GRIDDESC.

#Merge them all in a NetCDF file:
# Crear NetCDF base:
variables=("OPEN" "SURF" "CHLO" "DMS")
var_list=$(printf "%-16s" ${variables[@]})

#Crear estructura del NetCDF de salida
cat  > netcdf_emission_template.cdl <<EOF
netcdf ocean_1 {
dimensions:
    TSTEP = 1; DATE-TIME = 2; COL = $nx; ROW = $ny; LAY = 1; VAR= ${#variables[@]};
variables:
  float OPEN(TSTEP, LAY, ROW, COL);
  	OPEN:long_name = "OPEN            ";
  	OPEN:units = "UNKNOWN         ";
  	OPEN:var_desc = "OPEN                                                                            " ;
  float SURF(TSTEP, LAY, ROW, COL);
  	SURF:long_name = "SURF            " ;
  	SURF:units = "UNKNOWN         " ;
  	SURF:var_desc = "SURF                                                                            " ;
  float DMS(TSTEP, LAY, ROW, COL);
  	DMS:long_name = "DMS             " ;
  	DMS:units = "nM              " ;
  	DMS:var_desc = "DMS                                                                             " ;
  float CHLO(TSTEP, LAY, ROW, COL);
  	CHLO:long_name = "CHLO            " ;
  	CHLO:units = "mg m^-3         " ;
  	CHLO:var_desc = "Chlorophyll Concentration, OCI Algorithm                                        " ;
  int TFLAG(TSTEP, VAR, DATE-TIME) ;
  	TFLAG:units = "<YYYYDDD,HHMMSS>" ;
  	TFLAG:long_name = "TFLAG           " ;
  	TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
  float Band1(ROW, COL);
   //int TFLAG(TSTEP, VAR, DATE_TIME) ;
   //float OPEN(TSTEP, LAY, ROW, COL);
   //float SURF(TSTEP, LAY, ROW, COL);
   //float CHLO(TSTEP, LAY, ROW, COL);
   //float DMS(TSTEP, LAY, ROW, COL);
// global attributes:
 :IOAPI_VERSION = "ioapi-3.2: \$Id: init3" ; :EXEC_ID = "???????????????? " ; 
 :FTYPE = 1 ; :SDATE = -635; :STIME = 000000 ; :WDATE = 2023001 ; :WTIME = 000000 ; :CDATE = 2023001 ; :CTIME = 000000;
 :TSTEP = 0; :NTHIK = 1;    
 :NCOLS = ${nx}; :NROWS = ${ny}; :NLAYS = ${nz}; :NVARS = ${#variables[@]}; :GDTYP = ${COORDTYPE};
 :P_ALP = ${P_ALP};:P_BET = ${P_BET} ;:P_GAM = ${P_GAM}; 
 :XCENT = ${XCENT}; :YCENT = ${YCENT}; :XORIG = ${xorig}; :YORIG = ${yorig}; :XCELL = ${dx};  :YCELL = ${dy};  
 :VGTYP = -9999 ;:VGTOP = 5000.f ;:VGLVLS = 1.f, 0.9938147f;    
 :GDNAM = "${GRIDNAME}";:UPNAM = "CONVERT" ;   :VAR-LIST = "${var_list}";      
 :FILEDESC = "OCEAN file" ; :HISTORY = "" ;
}
EOF

for mm in $(seq --format="%02.0f" 1 12) 
do
   #Creo NetCDF destino:
   file_out=OCEAN_${GRIDNAME}_${mm}.nc 
   echo "Creating: ${file_out} .."
   ncgen -o $file_out netcdf_emission_template.cdl
   #OPEN
   ncks -4 -A -h -v Band1 open_${GRIDNAME}.nc ${file_out}
   ncap2 -A -h -s "OPEN(0,0,:,:)=float(Band1(:,:))" $file_out $file_out
   #SURF
   ncks -4 -A -h -v Band1 surf_${GRIDNAME}.nc ${file_out}
   ncap2 -A -h -s "SURF(0,0,:,:)=float(Band1(:,:))" $file_out $file_out
   ##CHLO
   gdalwarp -dstnodata 99999 -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r average -f "NetCDF" NETCDF:${chlo_nc}:chlo${mm} tmp0.nc
   ncks -4 -A -h -v Band1 tmp0.nc ${file_out}
   ncap2 -A -h -s "CHLO(0,0,:,:)=float(Band1); where (CHLO == 99999 ) CHLO = 0" $file_out $file_out
   ncatted -h -O -a missing_value,CHLO,o,f,99999 ${file_out}
   rm tmp0.nc
   ##DMS
   gdalwarp -dstnodata -999 -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r average -f "NetCDF" NETCDF:${dms_nc}:dms${mm} tmp0.nc
   ncks -4 -A -h -v Band1 tmp0.nc ${file_out}
   ncap2 -A -h -s "DMS(0,0,:,:)=float(Band1); where ( DMS < 0 ) DMS = 0" $file_out $file_out
   #TFLAG
   ncap2 -O -h -s "TFLAG(:,:,:) = 0"  $file_out $file_out

   ncks -O -x -h -v Band1,mercator ${file_out} ${file_out}
   ncatted -h -O -a missing_value,global,o,s,"" ${file_out}
done
rm tmp0.nc 
