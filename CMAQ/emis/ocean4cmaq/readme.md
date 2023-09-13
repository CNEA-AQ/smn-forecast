# Ocean4CMAQ

> Ocean data preprocessing for the CMAQ model.

## Dependencies
+ Fortran GNU compiler
+ NetCDF library
+ GDAL/OGR library

## Input-data

All the data required to run this pre-processor is freely available from [Ocean global data](https://drive.google.com/drive/folders/1j3efDayKluycRKPnqBvLZd33fxnbhRQq?usp=sharing).

Data required:
 + Global land shapefile (GPKG).
 + Global monthly chloropyll (CHLO) concentration.
 + Global monthly dimethil-sulfide (DMS) concentration.


## Build

None.

## Run

```
./ocean4cmaq.sh
```

## Planned future improvements

