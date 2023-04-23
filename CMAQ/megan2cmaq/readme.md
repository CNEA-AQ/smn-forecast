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

+ CTS `MEGAN_CTS` (*Canopy Type Fractions*) file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains canopy fraction information for six canopy types in one variable, CTS, which is nondimensional and ranges from 0-100. The vegetation types are: needleleaf trees, tropical forest trees, temperate broadleaf trees, shrubs, herbaceous plants, and crops.

+ LDF `MEGAN_LDF` (*Light Dependence Fractions*) file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains nondimensional light dependence fractions for 4 of the 19 MEGAN chemical species.

+ EF `MEGAN_EFS`  (emission factors). The MEGAN_EFS file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains emission factors for the 19 MEGAN chemical species.

+ LAI `MEGAN_LAI` (Leaf Area Index). The MEGAN_LAI file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains leaf area index that is separate from LAI values used in the rest of CMAQ. By default MEGAN will use this file for LAI, but users can choose to use the LAI values that are read in from MCIP files by setting the environmental variable USE_MEGAN_LAI to N in their run script.


+ BDSNP_AFILE: arid flag. Used by: CCTM online MEGAN biogenics emissions' BDSNP soil nitrogen model option. The BDSNP_AFILE file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor for use with the BDSNP soil nitrogen option. It identifies climatically arid grid cells with 1s and 0s.
+ BDSNP_NAFILE: nonarid flag. Used by: CCTM online MEGAN biogenics emissions' BDSNP soil nitrogen model option. The BDSNP_NAFILE file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor for use with the BDSNP soil nitrogen option. It identifies climatically non-arid grid cells with 1s and 0s.
+ BDSNP_LFILE: landfile type. Used by: CCTM online MEGAN biogenics emissions' BDSNP soil nitrogen model option. The BDSNP_LFILE file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor for use with the BDSNP soil nitrogen option. It assigns each grid cell to one of 24 land types.
+ BDSNP_FFILE: fertilizer reservoir. Used by: CCTM online MEGAN biogenics emissions' BDSNP soil nitrogen model option. The BDSNP_FFILE file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor for use with the BDSNP soil nitrogen option. It contains daily fertilizer information in ng N/m2 using 366 variables.
+ BDSNP_NFILE: nitrogen deposition. Used by: CCTM online MEGAN biogenics emissions' BDSNP soil nitrogen model option. The BDSNP_NFILE file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor for use with the BDSNP soil nitrogen option. It contains monthly average nitrogen deposition values in ng/m2/s using 12 variables.


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

