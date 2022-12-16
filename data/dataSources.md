# Fuentes de datos

> Fuentes de datos necesarios para ejecución de modelos.

Buscar:
- buscar variables relevantes.
- encontrar fuente de dato.
- desarrollar forma de descargarlo del servidor original en tiempo real.
- procesarlo.
- mejorarlo? (mejor resol o representatividad)

---
### Meteorología:

- NCAR UCAR
	+ Global Forecast System [GFS](https://www.nco.ncep.noaa.gov/pmb/products/gfs/):
		- [1°   4/day: 00, 06, 12, 18UTC  +00, (+03, +06 precipitation fields)](https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-003-1.0-degree/analysis/)
		- [0.5° 4/day: 00, 06, 12, 18UTC  +00, (+03, +06 precipitation fields)](https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-004-0.5-degree/analysis/)
		- NCEP server forecast data for HYSPLIT:: ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/hysplit/prod/
		- NCEP NOMADS (Hyplit): https://nomads.ncep.noaa.gov/pub/data/nccf/com/hysplit/prod/hysplit.20221028/
		- NOAA Operational Model Archive and Distribution System (NOMADS): https://nomads.ncep.noaa.gov/
	+ Final Analysis (FNL) / Global Data Assimilation System (GDAS): Incorpora datos observacionales, pero está listo 1 hora más tarde que GFS.
- ECMWF:
	+ Integrated Forecast System (IFS)

- ARL: Air Resources Laboratory
	+ forecast data: ftp.arl.noaa.gov/pub/forecast/
	+ archive data: ftp.arl.noaa.gov:/pub/archives/ 

---

### Emisiones:	
> Muchos inventarios se pueden obtener en: [ECCAD]( https://eccad3.sedoo.fr/ ) (Emissions of atmospheric Compound and Compilation of Ancillary Data)


#### Antropogénicas
- EDGAR: (Emissions Database for Global Atmospheric Research)[https://edgar.jrc.ec.europa.eu/]
- EDGAR-HTAP  
- RETRO

#### Fuegos:
- FINN (Fire INventory from NCAR) Se procesan con fire_emis para wrf.
	- https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/	(near real time)
	- https://www.acom.ucar.edu/Data/fire/				(hasta 2021)

- GFED (Global Fire Emissions Database)

- GFAS (Global Fire Assimilation System)

#### Biogénicas
-  MEGAN
-  BEIS (Biogenic Emission Inventory System)

#### Cenizas
- VOLC_SO2

#### Sea-salt
  SEAC4RS  

  GOCART
  OLSON2 

---
### Global CTMS (para ICON & BCON)

- CAMS: https://ads.atmosphere.copernicus.eu/cdsapp
	+ Tiene forecast, se baja grib que se puede convertir a netcdf con `cdo -f nc copy <file>.grib file.nc`
- GEOS-5: https://www.nccs.nasa.gov/services/data-collections/coupled-products/geos5-forecast
- SILAM
- Testiar si para CMAQ pueden usarse corridas anteriores como ICON y BCON!

---
### Otras variables:

- Subgrid Orography information for Gravity Wave Drag (OROGWD)
- Variance of Subgrid Scale Orography (VAR-SSO), para drag orográfico.
- Topography (DEM)
- Solpe Index (Islope)

- Crop (growing season, harvest date, planting date, crop type:(soy,wheat,cotton,corn,crop) )
- Erodability (para dust emiss)
- ClayFraction (para dust emiss)
- Albedo
- Snow Albedo (MODIS)
- greenfrac fpar (MODIS)
- LAI
- groundwater: transmisivity e-folding depth, recharge, riverbed elevation, water table depth
- lake depth
- LANDUSE 
- SOIL TYPE (top & bottom)
- Urban fraction (urbfrac)

Variables que no vi, pero que pienso que pueden ser relevantes:
- Sea Surface Temperature?
- Soil moisture?
- Surface reflectance.
- Solar/Actinic flux
