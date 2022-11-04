#!/bin/csh -f

set exitstat = 0
set here = `pwd`
set tmpfile = .install

/bin/rm -rf $tmpfile

# Check that downloaded files are available
set file = smoke_v48.nctox.data.tar.gz
echo "Checking for $file file..."

ls $file > $tmpfile
if ( $status > 0 ) then
    echo "ERROR: Could not find $file"
    echo "       Please download the file $file and try again"
    set exitstat = 1
endif

set file = smoke_v48.Linux2_x86_64ifort.tar.gz
echo "Checking for $file file..."

ls $file > $tmpfile
if ( $status > 0 ) then
    echo "ERROR: Could not find $file"
    echo "       Please download the file $file and try again"
    set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Check that SMK_HOME is set
if ( ! $?SMK_HOME ) then
    echo "ERROR: You must define the SMK_HOME environment variable"
    echo "       before running this script.  Use the following command"
    echo "       to set your SMK_HOME variable:"
    echo "          setenv SMK_HOME <your chosen SMOKE installation dir>"
    exit( 1 )
endif

echo "SMOKE v3.1 will be installed in the following directory:"
echo "      $SMK_HOME"
echo " "

# Install files
cd $SMK_HOME

set file = smoke_v48.nctox.data.tar.gz
echo "Unpacking file $file..."
tar xzvf $here/$file
if ( $status > 0 ) then
   echo "ERROR: Could not unpack the file $file"
   echo "       Confirm that the 'tar' command is available and"
   echo "       that you have the correct permissions"
   set exitstat = 1
endif

set file = smoke_v48.Linux2_x86_64ifort.tar.gz
echo "Unpacking file $file..."
tar xzvf $here/$file
if ( $status > 0 ) then
   echo "ERROR: Could not unpack the file $file"
   echo "       Confirm that the 'tar' command is available and"
   echo "       that you have the correct permissions"
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Source the assigns file
cd subsys/smoke/assigns
set file = ASSIGNS.nctox.cmaq.cb05_soa.us12-nc 
source $file 
if ( $status > 0 ) then
   echo "ERROR: Could not source the Assigns file $file"
   echo "       Please contact the CMAS Help Desk at http://www.cmascenter.org"
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit( $exitstat )
endif

# Create the inventory list files
echo "Creating inventory list files..."

cd $INVDIR/area 
echo "#LIST" > arinv.area.lst
ls $INVDIR/area/arinv_lm_no_c3_cap2002v3_orl_nc.txt >> arinv.area.lst
ls $INVDIR/area/arinv_nonpt_pf4_cap_nopfc_orl_nc.txt >> arinv.area.lst

cd $INVDIR/edgar
echo "#LIST GRID" > arinv.edgar.lst
echo "#SCC,InvPollName,EDGAR_VarName,Month,FileName" >> arinv.edgar.lst
echo "#DESC Sample test case for TRANSPORT EDGAR sector based on HTAP version 2.0" >> arinv.edgar.lst
echo "TRANSPORT,CO,emi_co,0,"$INVDIR/edgar/edgar_HTAP_CO_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,NOX,emi_nox,0,"$INVDIR/edgar/edgar_HTAP_NOx_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,SO2,emi_so2,0,"$INVDIR/edgar/edgar_HTAP_SO2_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,NH3,emi_nh3,0,"$INVDIR/edgar/edgar_HTAP_NH3_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,PM2_5,emi_pm2.5,0,"$INVDIR/edgar/edgar_HTAP_PM2.5_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,PM10,emi_pm10,0,"$INVDIR/edgar/edgar_HTAP_PM10_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst
echo "TRANSPORT,VOC,emi_nmvoc,0,"$INVDIR/edgar/edgar_HTAP_NMVOC_emi_TRANSPORT_2010.0.1x0.1.nc >> arinv.edgar.lst

cd $INVDIR/nonroad
echo "#LIST" > arinv.nonroad.lst
ls $INVDIR/nonroad/arinv_nonroad_caps_2005v2_jul_orl_nc.txt >> arinv.nonroad.lst

cd $INVDIR/point
echo "#LIST" > ptinv.point.lst
ls $INVDIR/point/ptinv_ptipm_cap2005v2_orl_nc.txt >> ptinv.point.lst
ls $INVDIR/point/ptinv_ptnonipm_xportfrac_cap2005v2_orl_nc.txt >> ptinv.point.lst

cd $INVDIR/rateperdistance
echo "#LIST" > mbinv.rateperdistance.lst
ls $INVDIR/rateperdistance/mbinv.VMT.nc.txt >> mbinv.rateperdistance.lst 
ls $INVDIR/rateperdistance/mbinv.SPEED.nc.txt >> mbinv.rateperdistance.lst 

cd $INVDIR/rateperhour
echo "#LIST" > mbinv.rateperhour.lst
ls $INVDIR/rateperhour/mbinv.HOTELLING.nc.txt >> mbinv.rateperhour.lst

cd $INVDIR/ratepervehicle
echo "#LIST" > mbinv.ratepervehicle.lst
ls $INVDIR/ratepervehicle/mbinv.VPOP.nc.txt >> mbinv.ratepervehicle.lst 

mkdir $INVDIR/rateperprofile
cd $INVDIR/rateperprofile
ln -s ../ratepervehicle/mbinv.ratepervehicle.lst mbinv.rateperprofile.lst

# Return to original directory
cd $here
rm -rf $tmpfile

echo "Installation completed successfully."
echo " "
echo "Please follow the instructions in Section 4.3 of the SMOKE User's Manual"
echo "to run the nctox test case."
echo "https://www.cmascenter.org/smoke/documentation/4.6/html/"
echo " "

exit( 0 )

