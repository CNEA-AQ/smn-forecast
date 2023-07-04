#================================================
# SMOKE
#	Build & Run.
#================================================
# Dependencias:
#	* I/O API 
#	* NetCDF (libnetcdf y libnetcdff)
#------------------------------------------------
# CAMINO OPCIONAL: descargar de la web de CMAS:
# - smoke_v48.Linux2_x86_64ifort.tar.gz # src de SMOKE
# - smoke_v48.nctox.data.tar.gz		# datos para testcase
# - smoke_install_v48.csh		# script para "instalar"
# Ejecutar:
#	csh smoke_install_v48.csh
# va a descomprmir los .tar.gz, donde vienen los exes compilados (con intel!)
#------------------------------------------------
# Build
git clone https://github.com/CEMPD/SMOKE
cd SMOKE
SMK_HOME=`pwd`
ln -s ~/m/libs/ioapi-3.2-20200828 ioapi    #link simbolco a donde instale ioapi
mkdir Linux2_x86_64gfort  #carpeta donde van a ir los *.o *.mod y exes
mkdir bin	          #carpeta donde van a ir los ejecutables
cd src
#Editar "Makeinclude":
```
BASEDIR=${SMK_HOME}/src
INCDIR=${BASEDIR}/inc
OBJDIR=${SMK_HOME}/${BIN}

IOBASE=${SMK_HOME}/ioapi
IODIR=${IOBASE}/ioapi
IOBIN=${IOBASE}/${BIN}
IOINC=${IODIR}/fixed_src

INSTDIR=${SMK_HOME}/bin
```
#>(!) No olvidar descomentar flags para GNU compilers.

#make, make install:
make FC=gfortran BIN=Linux2_x86_64gfort SMK_HOME=${SMK_HOME}
make install FC=gfortran BIN=Linux2_x86_64gfort SMK_HOME=${SMK_HOME}

#(!) Modificar ../ioapi/ioapi/fixed_src/PARMS3.EXT linea 104 (sintaxis vieja corte de linea array)
#(!) Modificar                     emqa/wrrephdr.f linea 287 (string de distinta longitud)
#(!) Si se está compilando con liberias locales, y ld no encuentra alguna librería: revisar que no quede en el Makeniclude (de smoke ó ioapi) ningun flag -Bstatic ó -Static
#(!) Tuve que agregar de flag al compilador de fortran "-fallow-argument-mismatch".

#------------------------------------------------
# Run:
#(todo en csh)

#---------------
#TEST CASE:

cd $SMK_HOME

tar -xzvf smoke_v48.nctox.data.tar.gz	#descomprimo tar con datos para testcase

setenv SMK_HOME /home/rama/CMAS/SMOKE

source $SMK_HOME/assings/ASSIGNS.nctox.cmaq.cb05_soa.us12-nc

cd $SCRIPTS/run #(change to the run scripts directory)

smk_area_nctox.csh 		#(invoke the stationary area run script)
smk_bg_nctox.csh 		#(invoke the BEIS3 biogenic run script)
smk_nonroad_nctox.csh 		#(invoke the nonroad mobile run script)
smk_point_nctox.csh 		#(invoke the point run script)
smk_rateperdistance_nctox.csh 	#(invoke the MOVES mobile sources on-roadway rate-per-distance (RPD) run script for all processes)
smk_ratepervehicle_nctox.csh 	#(invoke the MOVES mobile sources off-network rate-per-vehicle (RPV) run script)
smk_rateperprofile_nctox.csh 	#(invoke the MOVES mobile sources off-network rate-per-profile (RPP) run script)
smk_rateperhour_nctox.csh 	#(invoke the MOVES mobile sources off-network rate-per-hour (RPH) run script)
smk_mrgall_nctox.csh 		#(invoke the all-sources merge script)

#To verify that the nctox scripts have run correctly, go to the log file directory and look for errors by using:
cd $LOGS (change to the log file directory for the test case)
grep ERROR *

#If there are no errors, the next step is to run a QA script to be sure that the answers match.
cd $SCRIPTS/install 	#(change to the install directory)
check_smk_install 	#(invoke the smoke install quality assurance script)

#Los archivos de salida están en: 
ls $SMK_HOME/data/run_nctox/output/

#---------------
#Run with EDGAR data:
#To process pre-gridded data from different modeling domains, the following steps must be taken.
#1. Set the IMPORT_GRDNETCDF_YN to Y to read the native netCDF annual inventory files.
setenv IMPORT_GRDNETCDF_YN Y

#2. Grdmat is used to convert the data from the input data grid (e.g. lat-lon for EDGAR) to the output modeling grid. Set the IOAPI_GRIDNAME_1 to specify the output grid.
Grdmat

#3. SMOKE reads each grid cell in the gridded inventory file as a unique source that will have a unique source ID. The lower left corner of inventory grid will be source 1, cell (2,1) will be source 2, cell (3,1) source 3, etc

#4. For each gridded inventory file, a Geographical Code Mask Input file (GRIDMASK) is provided by UNC for EDGAR inventory files that assigns a GEOGRID_LEVEL2 value to each grid cell.

#5. The ARINV file is used to read in a list of gridded inventory files and to specify the SCC associated with each gridded inventory file

#6. Set the NETCDF_POL_UNIT to kg m-2 -s

#7. Set the NETCDF_INV_YEAR to the year of the inventory (e.g. 2010)

#8. To apply control factors to adjust emissions by GEOCODE_LEVEL and SCC use the CONTROL packet in the GCNTL file

#9. To apply vertical profiles to gridded inventory data using the SMOKE Utility Program Layalloc. The EDGAR inventory includes fixed vertical profiles by sector that define the percentages of emissions

#=================================================
#DATA DEL MANUAL:
#---------------------
# Inventarios:
# * "data types" : (polluts) ejemplo CO, NH4, Hg
# * "source cateogories": point, non-point, non-road mobile, on-road mobile, wildfire, biogenic.
# * "file formats": orl, ff10 (ascii, tipo csv) y netcdf.
#---------------------

# Tipos de inventarios:
#  + "inventarios criterio": que contienen polluts criterio de EPA: (CO,NOx,VOC,TOG)
#  + "inventarios de particulado": contienen NH3, SO2, PM,10, PM2.5

# Inventarios qeue SMOKE puede procesar:
#	* NEI (National Emission Inventory) of Hazardous Air Pollutants (HAPs)

# Inventory source categories:
# Los inventarios dividen las emisiones en categorías:
# + Non-point (Stationary Area)
# + Non-road mobile
# + On-road mobile: (MOVES)
# + Point sources
# + Wildfire
# + Biogenic: Biogenic emissions Landcover Database (BELD3)
#
#	 A su vez cada cateogria tiene:
#		+ caracteristicas
#		+ atributos de las fuentes
#		+ valores de datos
#		+ atributos de los datos


# Inventory file formats
# Formatos para cada categoría:
# + Non-point (Stationary Area):   ORL ó FF10 (Flat file 10), son ascii sep por "," ó ";" con un header.
# + Non-road mobile: 
# + On-road mobile: (MOVES)
# + Point sources
# + Wildfire
# + Biogenic: Biogenic emissions Landcover Database (BELD3)

#Importar inventarios:
#
# Smkinven: para importar data de fuentes antropogénicas
# Normbeis3: para importar data de fuentes biogenicas
#
