#!/bin/csh -f

# ======================= ICONv5.3.X Run Script ========================
# Usage: run.icon.csh >&! icon.log &                                   
# ==================================================================
#> Runtime Environment Options
# ==================================================================
#Path a donde estácompilado CMAQ.
setenv CMAQ_HOME /home/ramiroespada/shared_divqa_data/m/CMAQ

set APPL       = prueba
set CoordName  = LCC_TAN_ARG #> como quieran nombrar a la proyeccion. 16-character maximum
set GridName   = ARG_GRD     #> como quieran nombrar a la grilla    . 16-character maximum

set DataPath   = /home/ramiroespada/shared_divqa_data/tests
set MCIP_DATA  = $DataPath/cmaq/mcip
set OUTDIR     = $DataPath/cmaq/icon

#Ubicacion de ejecutable y nombre de ejecutable:
set VRSN       = v532        #> Code Version 
set compiler   = gcc         #> compilador
set BLD        = ${CMAQ_HOME}/PREP/icon/scripts/BLD_ICON_${VRSN}_${compiler}
set EXEC       = ICON_${VRSN}.exe  

#Fecha
set DATE       = "2022-10-19"  #fecha de cond iniciales

#> Set General Parameters for Configuring the Simulation
set ICTYPE   = profile #regrid         #> Initial conditions type [profile|regrid]
#
# ==================================================================
#> Set the working directory:
cat $BLD/ICON_${VRSN}.cfg; echo " "; set echo

#> Horizontal grid definition 
setenv GRID_NAME $GridName   #> como quieran nombrar a la grilla    . 16-character maximum
setenv GRIDDESC ${MCIP_DATA}/GRIDDESC
setenv IOAPI_ISPH 20                     #> GCTP spheroid, use 20 for WRF-based modeling

#> I/O Controls
setenv IOAPI_LOG_WRITE F     #> turn on excess WRITE3 logging [ options: T | F ]
setenv IOAPI_OFFSET_64 YES   #> support large timestep records (>2GB/timestep record) [ options: YES | NO ]
setenv EXECUTION_ID $EXEC    #> define the model execution id

# =====================================================================
#> ICON Configuration Options
#
# ICON can be run in one of two modes:                                     
#     1) regrids CMAQ CTM concentration files (IC type = regrid)     
#     2) use default profile inputs (IC type = profile)
# =====================================================================
setenv ICON_TYPE ` echo $ICTYPE | tr "[A-Z]" "[a-z]" ` 

# =====================================================================
#> Input/Output Directories
# =====================================================================
#> Input Files
#  
#  Regrid mode (IC = regrid) (includes nested domains, windowed domains,
#                             or general regridded domains)
#     CTM_CONC_1 = the CTM concentration file for the coarse domain          
#     MET_CRO_3D_CRS = the MET_CRO_3D met file for the coarse domain
#     MET_CRO_3D_FIN = the MET_CRO_3D met file for the target nested domain 
#                                                                            
#  Profile Mode (IC = profile)
#     IC_PROFILE = static/default IC profiles 
#     MET_CRO_3D_FIN = the MET_CRO_3D met file for the target domain 
#
# NOTE: SDATE (yyyyddd) and STIME (hhmmss) are only relevant to the
#       regrid mode and if they are not set, these variables will 
#       be set from the input MET_CRO_3D_FIN file
# =====================================================================
#> Output File
#     INIT_CONC_1 = gridded IC file for target domain
# =====================================================================

 set YYYYJJJ  = `date -ud "${DATE}" +%Y%j`   #> Convert YYYY-MM-DD to YYYYJJJ
 set YYMMDD   = `date -ud "${DATE}" +%y%m%d` #> Convert YYYY-MM-DD to YYMMDD
 set YYYYMMDD = `date -ud "${DATE}" +%Y%m%d` #> Convert YYYY-MM-DD to YYYYMMDD

 if ( $ICON_TYPE == regrid ) then
    setenv CTM_CONC_1 0 #/work/MOD3EVAL/sjr/CCTM_CONC_v53_intel18.0_2016_CONUS_test_${YYYYMMDD}.nc
    setenv MET_CRO_3D_CRS ${MCIP_DATA}/METCRO3D_${APPL}.nc 
    setenv MET_CRO_3D_FIN ${MCIP_DATA}/METCRO3D_${APPL}.nc 
 endif

 if ( $ICON_TYPE == profile ) then
    setenv IC_PROFILE $BLD/avprofile_cb6r3m_ae7_kmtbr_hemi2016_v53beta2_m3dry_col051_row068.csv
    setenv MET_CRO_3D_FIN ${MCIP_DATA}/METCRO3D_${APPL}.nc 
 endif
 
 setenv INIT_CONC_1  "$OUTDIR/ICON_${VRSN}_${APPL}_${ICON_TYPE}_${YYYYMMDD} -v"
#>- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if ( ! -d "$OUTDIR" ) mkdir -p $OUTDIR

 ls -l $BLD/$EXEC; size $BLD/$EXEC
 unlimit
 limit

#> Executable call:
 time $BLD/$EXEC

 exit() 
