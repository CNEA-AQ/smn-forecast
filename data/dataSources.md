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
		- NCEP NOMADS (Hyplit):                   https://nomads.ncep.noaa.gov/pub/data/nccf/com/hysplit/prod/hysplit.20221028/
		- NOAA Operational Model Archive and Distribution System (NOMADS): https://nomads.ncep.noaa.gov/
	+ Final Analysis (FNL) / Global Data Assimilation System (GDAS): Incorpora datos observacionales, pero está listo 1 hora más tarde que GFS.
- ECMWF:
	+ Integrated Forecast System (IFS)

- ARL: Air Resources Laboratory
	+ forecast data: ftp.arl.noaa.gov/pub/forecast/
	+ archive data: ftp.arl.noaa.gov:/pub/archives/ 

Los archivos están en
y su nombre es "hysplit.t{cycle}z.{name}" donde "cycle" representa la UTC inicial del forecast (e.g. 00,06,12,18), y "name":
	gfsf     :: 1-deg 3P +240h (814 Mb) global forecast at one-degree resolution at 3 hour intervals at pressure levels out to +240 hours.
	gfslrf   :: 1-deg 6P +384h (251 Mb) long-range global forecast at one-degree resolution at 6 hour intervals at pressure levels and from forecast hours +240 to +384.
	gfs0p25f :: 0.25-deg 3S +189h (2500 Mb) global forecast at quarter-degree resolution at 3 hour intervals at hybrid levels out to +189 hours. The complete forecast is split over multiple files. The forecast start time must be 	selected. There are only 21 hours per file. Files with forecast hours beyond +84 are not available on NOMADS.	

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
	- https://rda.ucar.edu/datasets/ds312.9/ (rda, es con usuario)
- GFED (Global Fire Emissions Database)

- GFAS (Global Fire Assimilation System)

#### Biogénicas
-  MEGAN: 
	- https://bai.ess.uci.edu/megan/data-and-code/megan21 
-  BEIS (Biogenic Emission Inventory System). Solo funciona para US.

#### Cenizas
- VOLC_SO2

#### Sea-salt
  SEAC4RS  

---
### Global CTMS (para ICON & BCON)

- CAMS: https://ads.atmosphere.copernicus.eu/cdsapp
	+ Tiene forecast, se baja grib que se puede convertir a netcdf con `cdo -f nc copy <file>.grib file.nc`
- GEOS-5: https://www.nccs.nasa.gov/services/data-collections/coupled-products/geos5-forecast
- SILAM:
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
