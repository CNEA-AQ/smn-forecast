#!/bin/sh
#
# Extracts and calculates PM and other measured species from silam output.
#
set -e
set -u
export OMP_NUM_THREADS=20

inp_files="output/arg_*\.nc"       # silam's output netcdf files directory.
outdir="output_2"
files=($(ls -d ${inp_files}))  
year=2022; hours=$((${#files[@]}*24))

#ncdump: (asegurarse que las variables a extraer esten incluidas acá!)
#cnc_AVB0_gas    cnc_XYL_gas   cnc_IOLE_gas    cnc_XO2N_gas     cnc_PM_m6_0 ;	cnc_BVB0_gas    cnc_HCO3_gas  cnc_HCHO_gas     cnc_sslt_m20 ;		
#cnc_AVB1e0_gas  cnc_PAN_gas   cnc_OLE5_gas    cnc_NO_gas       cnc_TOLUENE_gas cnc_BVB1e0_gas  cnc_ROOH_gas  cnc_O3_gas       cnc_NO3_c_m3_0 ;
#cnc_AVB1e1_gas  cnc_HONO_gas  cnc_C5H8_2_gas  cnc_OH_gas       cnc_PM_m2_5     cnc_BVB1e1_gas  cnc_MGLY_gas  cnc_MEO2_gas     cnc_EC_m_50 ;
#cnc_AVB1e2_gas  cnc_H2O2_gas  cnc_CRES_gas    cnc_NH3_gas      cnc_dust_m_30   cnc_BVB1e2_gas  cnc_CRO_gas   cnc_CXO3_gas     cnc_NOX_gas ;
#cnc_AVB1e3_gas  cnc_FACD_gas  cnc_C5H8_gas    cnc_NH4NO3_m_70  cnc_dust_m1_5   cnc_BVB1e3_gas  cnc_PANX_gas  cnc_XO2_gas      cnc_AVB0_m_50 ;
#cnc_BENZENE_gas cnc_AACD_gas  cnc_ISPD_gas    cnc_sslt_m_05    cnc_dust_m6_0   cnc_SESQ_gas    cnc_ROR5_gas  cnc_NO3_gas      cnc_AVB1e0_m_50 ;
#cnc_AVB1e4_gas  cnc_PACD_gas  cnc_NTR_gas     cnc_sslt_m_50    cnc_dust_m20    cnc_O1D_gas     cnc_MEPX_gas  cnc_C2O3_gas     cnc_AVB1e1_m_50 ;
#cnc_AVB1e5_gas  cnc_PNA_gas   cnc_ALDX_gas    cnc_sslt_m3_0    cnc_N2O5_gas    cnc_TOL_gas     cnc_CO_gas    cnc_O_gas        cnc_AVB1e2_m_50 ;
#cnc_AVB1e6_gas  cnc_TO2_gas   cnc_ALD2_gas    cnc_sslt_m9_0    cnc_ETHA_gas    cnc_OPEN_gas    cnc_NO2_gas   cnc_AVB1e3_m_50  
#cnc_HO2_gas     cnc_SO2_gas   cnc_PAR5_gas    cnc_CH3Cl_gas    cnc_MEOH_gas    cnc_HNO3_gas    cnc_ETOH_gas  cnc_ETH_gas    

#Grouping:
  sslt_Fine="cnc_sslt_m_05*1e9f cnc_sslt_m_50*1e9f cnc_sslt_m3_0*1e9f"
sslt_Coarse="cnc_sslt_m9_0*1e9f cnc_sslt_m20*1e9f"
  dust_Fine="cnc_dust_m_30*1e9f cnc_dust_m1_5*1e9f"
dust_Coarse="cnc_dust_m6_0*1e9f"
  fire_Fine="cnc_PM_m2_5*1e9f"
fire_Coarse="cnc_PM_m6_0*1e9f"
avb_species="cnc_AVB0_m_50*1e9f cnc_AVB1e0_m_50*1e9f cnc_AVB1e1_m_50*1e9f cnc_AVB1e2_m_50*1e9f cnc_AVB1e3_m_50*1e9f"
bvb_species="cnc_BVB0_m_50*1e9f cnc_BVB1e0_m_50*1e9f cnc_BVB1e1_m_50*1e9f cnc_BVB1e2_m_50*1e9f cnc_BVB1e3_m_50*1e9f"
nh4_species="cnc_NH4NO3_m_70*80e6f" #cnc_NH415SO4_m_20*117e6f cnc_NH415SO4_m_70*117e6f 
#PM:
PM25_SPECIES="${sslt_Fine} ${dust_Fine} ${fire_Fine} ${avb_species} ${nh4_species} cnc_EC_m_50*1e9f" # ${bvb_species} cnc_SO4_m_20*96e6f cnc_SO4_m_70*96e6fmineral_m_50*1e9f"
PM10_SPECIES="${sslt_Coarse} ${dust_Coarse} ${fire_Coarse} cnc_NO3_c_m3_0*62e6f"

#air_dens=air_dens;
cat << EOF > scriptfile.tmp
cnc_CO_gas=cnc_CO_gas*28e6f;cnc_CO_gas@units ="ug/m3";
cnc_NO_gas=cnc_NO_gas*30e6f;cnc_NO_gas@units ="ug/m3";
cnc_NO2_gas=cnc_NO2_gas*46e6f;cnc_NO2_gas@units="ug/m3";
cnc_SO2_gas=cnc_SO2_gas*64e6f;cnc_SO2_gas@units="ug/m3";
cnc_O3_gas=cnc_O3_gas*48e6f;cnc_O3_gas@units ="ug/m3";
cnc_PM10=${PM10_SPECIES// /+};cnc_PM10@long_name="Concentration in air PM10";cnc_PM10@components="${PM10_SPECIES//cnc_/}";cnc_PM10@substance_name="PM10"; cnc_PM10@silam_amount_unit="kg";cnc_PM10@mode_distribution_type="NO_MODE";
cnc_PM2_5=${PM25_SPECIES// /+};cnc_PM2_5@long_name="Concentration in air PM2_5";cnc_PM2_5@components="${PM25_SPECIES//cnc_/}";cnc_PM2_5@substance_name="PM2_5";cnc_PM2_5@silam_amount_unit="kg";cnc_PM2_5@mode_distribution_type="NO_MODE";
EOF

sed -i 's/;/;\n/g' scriptfile.tmp

pairs_inout=""
for f in ${files[@]}
do
	of="${outdir}/$(basename $f)4"
	pairs_inout="$pairs_inout $f $of"
done
echo "${pairs_inout[@]}"

if [ ! -d $outdir ]; then mkdir $outdir; fi;

if [ -n "$pairs_inout" ]; then
    echo -n $pairs_inout | xargs -r -n 2 echo  ncap2 -4 -L5 -O -v -S scriptfile.tmp | srun --input=0 -J post --ntasks=1 --account=project_2004363 --cpus-per-task=128 --time 30:00 xargs -IXXX -t -P128 sh -c "XXX"
    cat scriptfile.tmp; echo Done!
fi

cat > $outdir/0PM.nc.ctl <<EOF
DSET ^arg_%y4%m2%d2.nc4
DTYPE NETCDF
UNDEF -999998980358144.
OPTIONS TEMPLATE
TDEF time ${hours} LINEAR      00:00Z01jan${year}  1hr
EOF

#Alternativa?
#for f in ${files[@]}
#do
#	of="${outdir}/$(basename $f)4"
#	echo "ncap2 -4 -L5 -O -v -S scriptfile.tmp $f $outfile"
#
#done | xargs -I{} -t -P8 sh -c '{}' | srun  --input=0 -J wrfmet --ntasks=1 --account=project_2004363 --cpus-per-task=128 --time 30:00 xargs -IXXX -t -P128 sh -c "XXX"
