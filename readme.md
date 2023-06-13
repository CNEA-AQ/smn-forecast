# SMN AQ-forecast.

> Repositorio con scripts y documentación para la implementación de ensamble de modelos de química atmosférica para pronóstico de calidad del aire para el Servicio Meteorologico Nacional (SMN).

### Especificaciones generales:
   + Período de pronóstico: 48hs.
   + Resolución temporal (``dt``): 1 hora.
   + Resolución espacial (``dx=dy``): 20km <!--4 km-->
   + Niveles verticales (``nvert``): 45, tope a 10 hPa. (**Revisar!**)
   + Anidado: No.
   + Condiciones iniciales (``ICON``): A determinar.
   + Condiciones de borde  (``BCON``): A determinar. Silam?
   + Mecanísmo químico (gases): A determinar.
   + Mecanísmo químico (aerosoles): A determinar.

#### Dominio y grilla:
   + [namelist.wps](./smn-wrf/namelist.wps)

#### Meteorlogía, parametrizaciones y configuración:
   + [namelist.input)](./smn-wrf/namelist.input)

#### Contaminantes “target”
   - Monóxido de carbono (CO)
   - Óxidos de nitrógeno (NOx)
   - Dióxido de azufre (SO2)
   - Material Particulado < 10um (PM10)
   - Material Particulado < 2.5um (PM25)
   - Ozono (O3)
   - Compuestos Organicos Volátiles (COVs)
   - *Black-carbon* (BC)

#### Modelos a implementar:
   + [Chimere](./CHIMERE)
   + [CMAQ](./CMAQ)
   + [Silam](./SILAM)
   + ~[WRF-Chem](./WRF)~ (descartado por que no se puede correr offline)
   + ~[MUSICA](./MUSICA)~ (descartado)
   + [EMEP](./EMEP)

---

#### "TO-DO" list:

[x] Compilar modelos 
   - [x] CMAQ
   - [x] CHIMERE
   - [x] EMEP
   - [x] SILAM
   - [x] ~WRF-CHEM~  (descartado, no se puede correr desacoplado)
   - [x] ~MUSICA~    (descartado, no se puede correr desacoplado)
[ ] Probar funcionamiento de modelos
   - [x] CMAQ
   - [ ] SILAM
   - [ ] CHIMERE
   - [ ] EMEP
   - [x] ~WRF-CHEM~ (descartado, no se correr desacoplado)
   - [x] ~MUSICA~   (descartado, no se correr desacoplado)
[ ] ICON : empieza con 0 ó Global CTM (silam?), y luego usa la simulacion anterior.
[ ] BCON : empieza con 0 ó Global CTM (silam?), y luego usa la simulacion anterior.
[ ] Conseguir inventarios
   - [x] Antropogénicas : CAMS+PAPILA
   - [x] Biogénicas     : MEGAN (UCI-BAI)
   - [x] Fuegos         : FINN (NCAR)
   - [x] Polvo	        : (incluido en CTM) erodability, landuse, lai, greenfrac
   - [x] Sea-salt       : (incluido en CTM) surface temp, salinty
[ ] Mecanísmo químico
    - Fase Gaseosa
      + [ ] CBM (CB4, CB5)
      + [ ] RADM
      + [ ] RACM
      + [ ] MOZART
    - Aerosoles
      + [ ] GOCART
      + [ ] MOCART
      + [ ] MOSAIC
      + [ ] DMAT/VBS
      + [ ] LMDz/MACC
      + [ ] OPIA/SORGAM
      + [ ] Aero5 (CMAQ)
[ ] Parametrizaciónes de PBL
   - [ ] Bulk-Richardson
   - [ ] Yonsei University (YSU) 
[ ] Meteorología
   - [x]  GFS
   - [ ]  GFS-SMN
   - [ ]  ECMWF
