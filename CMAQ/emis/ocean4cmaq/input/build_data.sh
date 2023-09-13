#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
   landpoly=shp/land_10m_global.gpkg	                 #polygons to make a landmask
  chlor_dir=/home/rama/data/ocean/MY1DMM_CHLORA_2019     #chlorophyl conc. files (from MODIS, MY1DMM)
   dms_file=dms_lana_2011.nc                             #dms conc. files (DMS climatology of Lana et al. (2011))
#------------------------------------------------
# CHLO
#
# I first download the data as geotiff (floating point) from:
# https://neo.gsfc.nasa.gov/view.php?datasetId=MY1DMM_CHLORA&year=2019
#
if [ ! -f CHLO_monthly_global.nc ]
then
   for mm in $(seq --format="%02.0f" 1 12)
   do
       m=${mm#0} #"$(printf "%d" ${mm})";
       file=${chlor_dir}/MY1DMM_CHLORA_2019-${mm}-01_rgb_3600x1800.FLOAT.TIFF

       gdal_translate -q -a_srs "$srsOut" -ot Float32 -of netCDF -co "FORMAT=NC4" -co "COMPRESS=DEFLATE" -co "ZLEVEL=9" ${file} tmp0.nc 
       ncks -4 -A -v Band1 tmp0.nc tmp.nc

       ncrename -v Band1,chlo${mm} tmp.nc
       ncatted -h -O -a units,chlo${mm},o,c,"mg m^-3" tmp.nc
       ncatted -h -O -a long_name,chlo${mm},o,c,"CHLO" tmp.nc
       ncatted -h -O -a standard_name,chlo${mm},o,c,"mass_concentration_of_chlorophyll_in_sea_water" tmp.nc
       ncatted -h -O -a var_desc,chlo${mm},o,c,"Chlorophyll Concentration, from MY1DMM MODIS product." tmp.nc
       ncatted -h -O -a missing_value,chlo${mm},o,s,"99999" tmp.nc
   done
   ncatted -h -O -a Source_Software,global,d,c,"" tmp.nc
   ncatted -h -O -a Conventions,global,d,c,"" tmp.nc
   ncatted -h -O -a history,global,o,c,"" tmp.nc
   ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
   
   mv tmp.nc CHLO_monthly_global.nc
   rm tmp0.nc
else
   echo "CHLO_monthly_global.nc already exists!"
fi
 
#------------------------------------------------
# DMS
#  I already have this data from SILAM dataset (dms_lana_2011.nc).
#  DMS is in M and i should convert it to nM
if [ ! -f DMS_monthly_global.nc ]
then

cat << EOF > scriptfile.tmp
   dms01=cnc_DMS_gas(0 ,:,:)*1e9;dms01@units="nM";
   dms02=cnc_DMS_gas(1 ,:,:)*1e9;dms02@units="nM";
   dms03=cnc_DMS_gas(2 ,:,:)*1e9;dms03@units="nM";
   dms04=cnc_DMS_gas(3 ,:,:)*1e9;dms04@units="nM";
   dms05=cnc_DMS_gas(4 ,:,:)*1e9;dms05@units="nM";
   dms06=cnc_DMS_gas(5 ,:,:)*1e9;dms06@units="nM";
   dms07=cnc_DMS_gas(6 ,:,:)*1e9;dms07@units="nM";
   dms08=cnc_DMS_gas(7 ,:,:)*1e9;dms08@units="nM";
   dms09=cnc_DMS_gas(8 ,:,:)*1e9;dms09@units="nM";
   dms10=cnc_DMS_gas(9 ,:,:)*1e9;dms10@units="nM";
   dms11=cnc_DMS_gas(10,:,:)*1e9;dms11@units="nM";
   dms12=cnc_DMS_gas(11,:,:)*1e9;dms12@units="nM";
EOF
   sed -i 's/;/;\n/g' scriptfile.tmp
   ncap2 -4 -L5 -O -v -S scriptfile.tmp ${dms_file} DMS_monthly_global.nc
else
   echo "DMS_monthly_global.nc already exists!"
fi

