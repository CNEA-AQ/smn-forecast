#!/bin/csh -f 

#-----------------------------------------------------------------------
# Set identification for input and output files.
#
#   APPL       = Application Name (tag for MCIP output file names)
#   CoordName  = Coordinate system name for GRIDDESC
#   GridName   = Grid Name descriptor for GRIDDESC
#   InMetDir   = Directory that contains input meteorology files
#   InGeoDir   = Directory that contains input WRF "GEOGRID" file to
#                provide fractional land-use categories if "LANDUSEF"
#                was not included in the WRFOUT files.
#   OutDir     = Directory to write MCIP output files
#   ProgDir    = Directory that contains the MCIP executable
#   WorkDir    = Working Directory for Fortran links and namelist
#-----------------------------------------------------------------------

#Path a donde estácompilado CMAQ.
setenv CMAQ_HOME /home/ramiroespada/shared_divqa_data/m/CMAQ


set APPL       = prueba
set CoordName  = LCC_TAN_ARG # como quieran nombrar a la proyeccion. 16-character maximum
set GridName   = ARG_GRD     # como quieran nombrar a la grilla    . 16-character maximum

set DataPath   = /home/ramiroespada/shared_divqa_data/tests          
set InMetDir   = $DataPath/met
set InGeoDir   = $DataPath/met
set OutDir     = $DataPath/cmaq/mcip
set ProgDir    = $CMAQ_HOME/PREP/mcip/src
set WorkDir    = $OutDir

#-----------------------------------------------------------------------
# Set name(s) of input meteorology file(s)
#
#   File name(s) must be set inside parentheses since "InMetFiles" is
#   a C-shell script array.  Multiple file names should be space-
#   delimited.  Additional lines can be used when separated by a
#   back-slash (\) continuation marker.  The file names can be as
#   they appear on your system; MCIP will link the files in by a
#   Fortran unit number and the explicit name via a namelist.  The
#   files must be listed in chronological order.  The maximum number
#   of input meteorology files must be less than or equal to the number
#   in MAX_MM in file_mod.F (default is 367).
#
#   Example:
#     set InMetFiles = ( $InMetDir/wrfout_d01_date1 \
#                        $InMetDir/wrfout_d01_date2 )
#
#-----------------------------------------------------------------------
set InMetFiles = ( $InMetDir/wrfout_d01_2022-10-19_18:00:00 )

set IfGeo      = "T"
set InGeoFile  = $InGeoDir/geo_em.d01.nc

#-----------------------------------------------------------------------
# Set user control options.
#
#   LPV:     0 = Do not compute and output potential vorticity
#            1 = Compute and output potential vorticity
#
#   LWOUT:   0 = Do not output vertical velocity
#            1 = Output vertical velocity
#
#   LUVBOUT: 0 = Do not output u- and v-component winds on B-grid
#            1 = Output u- and v-component winds on B-grid (cell corner)
#                in addition to the C-grid (cell face) output
#-----------------------------------------------------------------------

set LPV     = 0
set LWOUT   = 0
set LUVBOUT = 1

#-----------------------------------------------------------------------
# Set run start and end date.  (YYYY-MO-DD-HH:MI:SS.SSSS)
#   MCIP_START:  First date and time to be output [UTC]
#   MCIP_END:    Last date and time to be output  [UTC]
#   INTVL:       Frequency of output [minutes]
#-----------------------------------------------------------------------

set MCIP_START = 2022-10-19-19:00:00.0000  # [UTC]
set MCIP_END   = 2022-10-20-18:00:00.0000  # [UTC]

set INTVL      = 60 # [min]

#-----------------------------------------------------------------------
# Choose output format.
#   1 = Models-3 I/O API
#   2 = netCDF
#-----------------------------------------------------------------------

set IOFORM = 1

#-----------------------------------------------------------------------
# Set number of meteorology "boundary" points to remove on each of four
# horizontal sides of MCIP domain.  This affects the output MCIP domain
# dimensions by reducing meteorology domain by 2*BTRIM + 2*NTHIK + 1,
# where NTHIK is the lateral boundary thickness (in BDY files), and the
# extra point reflects conversion from grid points (dot points) to grid
# cells (cross points).  Setting BTRIM = 0 will use maximum of input
# meteorology.  To remove MM5 lateral boundaries, set BTRIM = 5.
#
# *** If windowing a specific subset domain of input meteorology, set
#     BTRIM = -1, and BTRIM will be ignored in favor of specific window
#     information in X0, Y0, NCOLS, and NROWS.
#-----------------------------------------------------------------------

set BTRIM = 0

#-----------------------------------------------------------------------
# Define MCIP subset domain.  (Only used if BTRIM = -1.  Otherwise,
# the following variables will be set automatically from BTRIM and
# size of input meteorology fields.)
#   X0:     X-coordinate of lower-left corner of full MCIP "X" domain
#           (including MCIP lateral boundary) based on input MM5 domain.
#           X0 refers to the east-west dimension.  Minimum value is 1.
#   Y0:     Y-coordinate of lower-left corner of full MCIP "X" domain
#           (including MCIP lateral boundary) based on input MM5 domain.
#           Y0 refers to the north-south dimension.  Minimum value is 1.
#   NCOLS:  Number of columns in output MCIP domain (excluding MCIP
#           lateral boundaries).
#   NROWS:  Number of rows in output MCIP domain (excluding MCIP
#           lateral boundaries).
#-----------------------------------------------------------------------

set X0    =  13
set Y0    =  94
set NCOLS =  89
set NROWS = 104

#-----------------------------------------------------------------------
# Set coordinates for cell for diagnostic prints on output domain.
# If coordinate is set to 0, domain center cell will be used.
#-----------------------------------------------------------------------

set LPRT_COL = 0
set LPRT_ROW = 0

#-----------------------------------------------------------------------
# Optional:  Set WRF Lambert conformal reference latitude.
#            (Handy for matching WRF grids to existing MM5 grids.)
#            If not set, MCIP will use average of two true latitudes.
# To "unset" this variable, set the script variable to "-999.0".
# Alternatively, if the script variable is removed here, remove it
# from the setting of the namelist (toward the end of the script).
#-----------------------------------------------------------------------

 set WRF_LC_REF_LAT = "-999.0"

#=======================================================================
#=======================================================================
# Set up and run MCIP.
#   Should not need to change anything below here.
#=======================================================================
#=======================================================================

set PROG = mcip

date

#-----------------------------------------------------------------------
# Make sure directories exist.
#-----------------------------------------------------------------------

if ( ! -d $InMetDir ) then
  echo "No such input directory $InMetDir"
  exit 1
endif

if ( ! -d $OutDir ) then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if ( $status != 0 ) then
    echo "Failed to make output directory, $OutDir"
    exit 1
  endif
endif

if ( ! -d $ProgDir ) then
  echo "No such program directory $ProgDir"
  exit 1
endif

#-----------------------------------------------------------------------
# Make sure the input files exist.
#-----------------------------------------------------------------------

if ( $IfGeo == "T" ) then
  if ( ! -f $InGeoFile ) then
    echo "No such input file $InGeoFile"
    exit 1
  endif
endif

foreach fil ( $InMetFiles )
  if ( ! -f $fil ) then
    echo "No such input file $fil"
    exit 1
  endif
end

#-----------------------------------------------------------------------
# Make sure the executable exists.
#-----------------------------------------------------------------------

if ( ! -f $ProgDir/${PROG}.exe ) then
  echo "Could not find ${PROG}.exe"
  exit 1
endif

#-----------------------------------------------------------------------
# Create a work directory for this job.
#-----------------------------------------------------------------------

if ( ! -d $WorkDir ) then
  mkdir -p $WorkDir
  if ( $status != 0 ) then
    echo "Failed to make work directory, $WorkDir"
    exit 1
  endif
endif

cd $WorkDir

#-----------------------------------------------------------------------
# Set up script variables for input files.
#-----------------------------------------------------------------------

if ( $IfGeo == "T" ) then
  if ( -f $InGeoFile ) then
    set InGeo = $InGeoFile
  else
    set InGeo = "no_file"
  endif
else
  set InGeo = "no_file"
endif

set FILE_GD  = $OutDir/GRIDDESC

#-----------------------------------------------------------------------
# Create namelist with user definitions.
#-----------------------------------------------------------------------

set MACHTYPE = `uname`
if ( ( $MACHTYPE == "AIX" ) || ( $MACHTYPE == "Darwin" ) ) then
  set Marker = "/"
else
  set Marker = "&END"
endif

cat > $WorkDir/namelist.${PROG} << !

 &FILENAMES
  file_gd    = "$FILE_GD"
  file_mm    = "$InMetFiles[1]",
!

if ( $#InMetFiles > 1 ) then
  @ nn = 2
  while ( $nn <= $#InMetFiles )
    cat >> $WorkDir/namelist.${PROG} << !
               "$InMetFiles[$nn]",
!
    @ nn ++
  end
endif

if ( $IfGeo == "T" ) then
cat >> $WorkDir/namelist.${PROG} << !
  file_geo   = "$InGeo"
!
endif

cat >> $WorkDir/namelist.${PROG} << !
  ioform     =  $IOFORM
 $Marker

 &USERDEFS
  lpv        =  $LPV
  lwout      =  $LWOUT
  luvbout    =  $LUVBOUT
  mcip_start = "$MCIP_START"
  mcip_end   = "$MCIP_END"
  intvl      =  $INTVL
  coordnam   = "$CoordName"
  grdnam     = "$GridName"
  btrim      =  $BTRIM
  lprt_col   =  $LPRT_COL
  lprt_row   =  $LPRT_ROW
  wrf_lc_ref_lat = $WRF_LC_REF_LAT
 $Marker

 &WINDOWDEFS
  x0         =  $X0
  y0         =  $Y0
  ncolsin    =  $NCOLS
  nrowsin    =  $NROWS
 $Marker

!

#-----------------------------------------------------------------------
# Set links to FORTRAN units.
#-----------------------------------------------------------------------

rm fort.*
if ( -f $FILE_GD ) rm -f $FILE_GD

ln -s $FILE_GD                   fort.4
ln -s $WorkDir/namelist.${PROG}  fort.8

set NUMFIL = 0
foreach fil ( $InMetFiles )
  @ NN = $NUMFIL + 10
  ln -s $fil fort.$NN
  @ NUMFIL ++
end

#-----------------------------------------------------------------------
# Set output file names and other miscellaneous environment variables.
#-----------------------------------------------------------------------

setenv IOAPI_CHECK_HEADERS  T
setenv EXECUTION_ID         $PROG

setenv GRID_BDY_2D          $OutDir/GRIDBDY2D_${APPL}.nc
setenv GRID_CRO_2D          $OutDir/GRIDCRO2D_${APPL}.nc
setenv GRID_DOT_2D          $OutDir/GRIDDOT2D_${APPL}.nc
setenv MET_BDY_3D           $OutDir/METBDY3D_${APPL}.nc
setenv MET_CRO_2D           $OutDir/METCRO2D_${APPL}.nc
setenv MET_CRO_3D           $OutDir/METCRO3D_${APPL}.nc
setenv MET_DOT_3D           $OutDir/METDOT3D_${APPL}.nc
setenv LUFRAC_CRO           $OutDir/LUFRAC_CRO_${APPL}.nc
setenv SOI_CRO              $OutDir/SOI_CRO_${APPL}.nc
setenv MOSAIC_CRO           $OutDir/MOSAIC_CRO_${APPL}.nc

if ( -f $GRID_BDY_2D ) rm -f $GRID_BDY_2D
if ( -f $GRID_CRO_2D ) rm -f $GRID_CRO_2D
if ( -f $GRID_DOT_2D ) rm -f $GRID_DOT_2D
if ( -f $MET_BDY_3D  ) rm -f $MET_BDY_3D
if ( -f $MET_CRO_2D  ) rm -f $MET_CRO_2D
if ( -f $MET_CRO_3D  ) rm -f $MET_CRO_3D
if ( -f $MET_DOT_3D  ) rm -f $MET_DOT_3D
if ( -f $LUFRAC_CRO  ) rm -f $LUFRAC_CRO
if ( -f $SOI_CRO     ) rm -f $SOI_CRO
if ( -f $MOSAIC_CRO  ) rm -f $MOSAIC_CRO

if ( -f $OutDir/mcip.nc      ) rm -f $OutDir/mcip.nc
if ( -f $OutDir/mcip_bdy.nc  ) rm -f $OutDir/mcip_bdy.nc

#-----------------------------------------------------------------------
# Execute MCIP.
#-----------------------------------------------------------------------

$ProgDir/${PROG}.exe

if ( $status == 0 ) then
  rm fort.*
  exit 0
else
  echo "Error running $PROG"
  exit 1
endif
