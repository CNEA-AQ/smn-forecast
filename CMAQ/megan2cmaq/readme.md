# MEGAN Emissions

> `megan2cmaq` is a fortran program that grids all the data needed to run biogenic emissions inside CMAQ model.

## Dependencies:
 +  Fortran GNU compiler
 +  NetCDF library
 +  GDAL/OGR programs (`sudo apt install gdal-bin`)

## Get MEGAN input data:

Data required:
+ Leaf Area Index / Vegetation Cover Fraction (LAIv)
+ Growth Form (fraction)			WRF classes:
   - crops	crop30s_reorder.nc		  - herbs
   - grass	gras30s_reorder.nca		
   - shurbs	shrb30s_reorder.nc		  - shurbs
   - tree	tree30s_reorder.nc		  - trees: broadleaf & needleleaf
    + tropical
    + boradleaf
    + needleleaf
+ Ecotype


+ BDSNP (*optional*): for soil NO algorithm.
   - Fertilizer
   - Land Fraction
   - Climate data	
   - Nitrogen deposition data

> \* Optional.


Output files:
	+
	+
	+


## Build
Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your `nc-config --libdir` and `nc-config --includedir`.

`> make`

If the compilation is successful, the executable `megan2cmaq.exe` should be created.

## Run

Edit the 'example.inp' that contains the following variables:

```
&control
  start_date="2019-01-01",       		!%Y-%m-%d (YYYY-MM-DD)
    end_date="2019-01-01",       		!%Y-%m-%d (YYYY-MM-DD)
  inp_data_directory = './megan_inp_data/'	!path to finn files directory
  griddesc_path="./GRIDDESC",   		!path to the GRIDDESC file
/
```

Note that the start_date, end_date, finn_data_directory, and griddesc_path variables must be adjusted to match the appropriate values for your system.

Then execute `megan2cmaq.exe`:

`> ./megan2cmaq.exe < example.inp ` 

Please feel free to contact the developer if you have any issues or suggestions.

