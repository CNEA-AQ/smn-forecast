#!/bin/bash

#---------------------------------------------------------------------------------
#	Architecture file for compiling and running CHIMERE	
#	Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------
## Modules
module purge

module load intel/intel_2015.3.187
module load ompi3.1.4_intel_2015.3.187 
module load netcdf4.4_intel_2015.3.187_parallel
module load hdf5_1.8.17_intel_2015.3.187
module load miscLibs_intel_2015.3.187
module load eccodes-2.23_intel_2015.3.187

export MiDIR=/home/ramiroespada
export MiLIBS=$MiDIR/libs_intel_2015.3.187
#---------------------------------------------------------------------------------
# 	Compilers
#---------------------------------------------------------------------------------
export my_compilerF90=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/ifort #`which ifort`	  # Path to Fortran 90 compiler
export my_compilerC=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/icc     #`which icc	`		  # Path to C compiler
export my_compilerCpp=/share/apps/intel/composer_xe_2015.3.187/bin/intel64/icpc  #`which icpc	`		# Path to C++ compiler

#---------------------------------------------------------------------------------
# 	MPI - parallel execution of chimere
#---------------------------------------------------------------------------------
export  my_mpiframe=ompi  		                            	# implementaion of MPI norm [ ompi / ccrt ] TO REMOVE
export  my_mpibin=$MiLIBS/ompi/bin    			# Path to MPI binary directory
export  my_mpirun=$MiLIBS/ompi/bin/mpirun   # Path to mpirun to execute parallel job in MPI
export  my_mpif90=$MiLIBS/ompi/bin/mpifort  # Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=$MiLIBS/ompi/bin/mpicc    # Wrapper to my_compilerC to link with MPI library
export  my_mpilib=$MiLIBS/ompi/lib    			# Path to MPI libraries directory
export  my_mpiinc=$MiLIBS/ompi/include    	# Path to MPI include files directory

#---------------------------------------------------------------------------------
# 	HDF5  - parallel version	
#---------------------------------------------------------------------------------
export my_hdflib=$MiLIBS/hdf5/lib		# Path to HDF5 parallel library directory
export my_hdfinc=$MiLIBS/hdf5/include	# Path to HDF5 parallel include files directory

#---------------------------------------------------------------------------------
# 	NETCDF-C  - link with HDF5 parallel 
#---------------------------------------------------------------------------------
export my_netcdfCbin=$MiLIBS/pnetcdf/bin 		# Path to NETCDF-C (linked with HDF5 parallel) binaries directory 
export my_netcdfClib=$MiLIBS/pnetcdf/lib		# Path to NETCDF-C (linked with HDF5 parallel) library directory

#---------------------------------------------------------------------------------
# 	NETCDF-Fortran  - link with HDF5 parallel and NETCDF-C
#---------------------------------------------------------------------------------
export my_netcdfF90bin=$MiLIBS/pnetcdf/bin      # PATH to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) binaries  directory
export my_netcdfF90lib=$MiLIBS/pnetcdf/lib		# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) library  directory
export my_netcdfF90inc=$MiLIBS/pnetcdf/include	# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) include files  directory

#---------------------------------------------------------------------------------
# 	GRIB  - link with jasper 
#---------------------------------------------------------------------------------
#export my_griblib=$MiLIBS/grib_api-1.15.0-Source/build/lib     	# Path to GRIB library directory
#export my_gribinc=$MiLIBS/grib_api-1.15.0-Source/build/include 	# Path to GRIB include files directory
#eccodes
export my_griblib=$MiLIBS/eccodes-2.23-Source/build/lib64     	# Path to GRIB library directory
export my_gribinc=$MiLIBS/eccodes-2.23-Source/build/include 	# Path to GRIB include files directory

export my_jasperlib=$MiLIBS/grib2/lib 	                # Path to JASPER library directory
export my_jasperinc=$MiLIB/grib2/include                # Path to JASPER include files directory
#---------------------------------------------------------------------------------
# 	BLITZ
#---------------------------------------------------------------------------------
export my_blitzinc=$MiLIBS/blitz/build/include #/blitz		 # Path to BLITZ include files 
#---------------------------------------------------------------------------------
# 	Utilities	
#---------------------------------------------------------------------------------
export my_make=make 	                        # Path to make 
export my_awk=awk		                          # Path to awk
export my_ncdump=$MiLIBS/pnetcdf/bin/ncdump		# Path to ncdump

#---------------------------------------------------------------------------------
# 	Makefile header needed to compile CHIMERE and WRF 
#	     - with this architecture configuration - 	
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.ifort-64-ompi   		            	# Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name=configure.wrf.ifort             	# Makefile header to compile WRF in config_wrf directory
export configure_wps_file_name=configure_ifort.wps          	# Makefile header to compile WPS in config_wps directory

#---------------------------------------------------------------------------------
#	Export of Shared Library to be available at run time 	
#---------------------------------------------------------------------------------
export LD_LIBRARY_PATH=${my_hdflib}:${my_netcdfF90lib}:${my_netcdfClib}:${my_griblib}:${my_mpilib}:${my_mpilib}/openmpi:$LD_LIBRARY_PATH
