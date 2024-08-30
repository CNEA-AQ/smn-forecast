#!/bin/csh -f
# ======================= BCONv5.3.X Run Script ======================== 
# Usage: run.bcon.csh >&! bcon.log &                                
# ==================================================================== 
#> Runtime Environment Options
# ==================================================================
#Path a donde este compilado CMAQ.
setenv CMAQ_HOME /home/ramiroespada/shared_divqa_data/m/CMAQ

set APPL       = prueba
set CoordName  = LCC_TAN_ARG   #> como quieran nombrar a la proyeccion. 16-character maximum
set GridName   = ARG_GRD       #> como quieran nombrar a la grilla    . 16-character maximum

set DataPath   = /home/ramiroespada/shared_divqa_data/tests
set MCIP_DATA  = $DataPath/cmaq/mcip
set OUTDIR     = $DataPath/cmaq/bcon

#Ubicacion de ejecutable y nombre de ejecutable:
set VRSN       = v532          #> Code Version 
set compiler   = gcc           #> compilador
set BLD        = ${CMAQ_HOME}/PREP/bcon/scripts/BLD_BCON_${VRSN}_${compiler}
set EXEC       = BCON_${VRSN}.exe  
#Fecha
set DATE       = "2022-10-19"  #fecha de cond iniciales

set BCTYPE     = profile       #> Boundary condition type [profile|regrid]
#
# ==================================================================
#> Set the build directory:
cat $BLD/BCON_${VRSN}.cfg; echo " "; set echo

#> Horizontal grid definition 
setenv GRID_NAME $GridName               #> check GRIDDESC file for GRID_NAME options
setenv GRIDDESC ${MCIP_DATA}/GRIDDESC    #> path to GRIDDESC file
setenv IOAPI_ISPH 20                     #> GCTP spheroid, use 20 for WRF-based modeling

#> I/O Controls
setenv IOAPI_LOG_WRITE F     #> turn on excess WRITE3 logging [ options: T | F ]
setenv IOAPI_OFFSET_64 YES   #> support large timestep records (>2GB/timestep record) [ options: YES | NO ]
setenv EXECUTION_ID $EXEC    #> define the model execution id

# =====================================================================
#> BCON Configuration Options
#
# BCON can be run in one of two modes:                                     
#     1) regrids CMAQ CTM concentration files (BC type = regrid)     
#     2) use default profile inputs (BC type = profile)
# =====================================================================

setenv BCON_TYPE ` echo $BCTYPE | tr "[A-Z]" "[a-z]" `

# =====================================================================
#> Input/Output Directories
# =====================================================================
# =====================================================================
#> Input Files
#  
#  Regrid mode (BC = regrid) (includes nested domains, windowed domains,
#                             or general regridded domains)
#     CTM_CONC_1 = the CTM concentration file for the coarse domain          
#     MET_CRO_3D_CRS = the MET_CRO_3D met file for the coarse domain
#     MET_BDY_3D_FIN = the MET_BDY_3D met file for the target nested domain
#                                                                            
#  Profile mode (BC type = profile)
#     BC_PROFILE = static/default BC profiles 
#     MET_BDY_3D_FIN = the MET_BDY_3D met file for the target domain 
#
# NOTE: SDATE (yyyyddd), STIME (hhmmss) and RUNLEN (hhmmss) are only 
#       relevant to the regrid mode and if they are not set,  
#       these variables will be set from the input MET_BDY_3D_FIN file
# =====================================================================
#> Output File
#     BNDY_CONC_1 = gridded BC file for target domain
# =====================================================================
set YYYYJJJ  = `date -ud "${DATE}" +%Y%j`   #> Convert YYYY-MM-DD to YYYYJJJ
set YYMMDD   = `date -ud "${DATE}" +%y%m%d` #> Convert YYYY-MM-DD to YYMMDD
set YYYYMMDD = `date -ud "${DATE}" +%Y%m%d` #> Convert YYYY-MM-DD to YYYYMMDD

if ( $BCON_TYPE == regrid ) then 
   setenv CTM_CONC_1     /path/to/gridded/bcon/CCTM_CONC_v53_intel18.0_2016_CONUS_test_${YYYYMMDD}.nc
   setenv MET_CRO_3D_CRS $MCIP_DATA/METCRO3D.12US1.35L.${YYMMDD}
   setenv MET_BDY_3D_FIN $MCIP_DATA/METBDY3D_${YYMMDD}.nc
   setenv BNDY_CONC_1    "$OUTDIR/BCON_${VRSN}_${APPL}_${BCON_TYPE}_${YYYYMMDD} -v"
endif

if ( $BCON_TYPE == profile ) then
   setenv BC_PROFILE $BLD/avprofile_cb6r3m_ae7_kmtbr_hemi2016_v53beta2_m3dry_col051_row068.csv
   setenv MET_BDY_3D_FIN $MCIP_DATA/METBDY3D_${APPL}.nc
   setenv BNDY_CONC_1    "$OUTDIR/BCON_${VRSN}_${APPL}_${BCON_TYPE}_${YYYYMMDD} -v"
endif

# =====================================================================
#> Output File
# =====================================================================
if ( ! -d "$OUTDIR" ) mkdir -p $OUTDIR

ls -l $BLD/$EXEC; size $BLD/$EXEC
unlimit
limit

#> Executable call:
time $BLD/$EXEC

exit() 
