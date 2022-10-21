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
	+ Global Forecast System [GFS](https://www.nco.ncep.noaa.gov/pmb/products/gfs/)
	+ Final Analysis (FNL) / Global Data Assimilation System (GDAS): Incorpora datos observacionales, pero está listo 1 hora más tarde que GFS.
- ECMWF:
	+ Integrated Forecast System (IFS)

---
### Emisiones:

- ECCAD: (Emissions of atmospheric Compound and Compilation of Ancillary Data)[https://eccad3.sedoo.fr/]

- EDGAR: (Emissions Database for Global Atmospheric Research)[https://edgar.jrc.ec.europa.eu/]

  biogenic_emissions
  MEGAN
  EDGAR-HTAP  
  EDGAR 
  EDGARV4 
  RETRO
  fires_data
  GFEDv2-8days  
  VOLC_SO2
  SEAC4RS  
  GOCART
  OLSON2 

---
### Global CTMS (para ICON & BCON)

- CAMS
- GEOS-5
- SILAM

---
### Otras variables:

- Subgrid Orography information for Gravity Wave Drag (OROGWD)
- Variance of Subgrid Scale Orography (VAR-SSO), para drag orográfico.
- Topography(DEM)
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
