#================================================
# SILAM: Build & Run
#================================================
#Dependencias:
#   * Fortran Compilers (gfortran)
#   * MPI Library (OpenMPI)
#   * HDF5 Library
#   * NetCDF Library
#   * GRIB_API (usa ECCODES en su lugar)
#   * Otros: (LAPACK/BLAS, proj)
#------------------------------------------------

#(0) Instalar dependencias:
sudo apt install make python gfortran libeccodes-dev libnetcdf-dev libnetcdff-dev liblapack-dev libblas-dev libbz2-dev libproj-dev

#(1) Traer Repositorio:
git clone https://github.com/fmidev/silam-model

#------------------------------------------------
#(2) Compilar
cd silam-model/source/
make gnu
make

#Deberia crearse un binario en:
#$ ../bin/silam_v5_7pub.gnu

#Si hay problemas con alguna libreria cambiar flags en ../build/options.gnu y volver a intentar.
#En ubuntu hay un error: "grib_api.mod not found" que en el readme.md del repo explican como resolver.
#------------------------------------------------
#(3) RUN

#Archivos de control:

 + *CONTROL FILE*: sets the user-defined parameters of the run;
 	+ *source term file:* describes the emission sources: referred from the control file;
 	+ *output configuration file:* description of the output setup, referred from the control file;
 	+ *internal model setup:* referred from the control file, sets the internal model features, usually read-only or fully invisible for users;
 		+ *standard cocktails file:* defines the standard cocktails that can be used in the source description; referred from the internal setup file. Users are free to create their own cocktails, adding to the existent file;
	+ GRIB or NetCDF code table (depending of the type of files): mandatory, invisible for
users, referred from the internal setup file.

#Archivos opcionales:
Depending of the configuration of the run, there are different files that should be included in the setup configuration:
	- nuclide data file: for radioactive simulations, invisible for users, referred from the internal setup file;
	- nuclide decay data file: for radioactive simulations, invisible for users, referred from the	internal setup file;
	- land-use data: for chemical simulations of biogenic emissions, invisible for users, referred from the internal setup file;
	- optical properties: for chemical and aerosol simulation, describes the optical properties of substances, invisible for users, referred from the internal setup file;
	- chemical properties: describes the chemical properties of the species available in SILAM, invisible for users, referred from the internal setup file;









