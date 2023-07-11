DSET ^.nc
DTYPE NETCDF
UNDEF -999998980358144.
XDEF lon      40 LINEAR  -10.0000    0.2500
YDEF lat      40 LINEAR   60.0000    0.2500
ZDEF height     9 LEVELS     12.5     50.0    125.0    275.0    575.0   1150.0   2125.0   3725.0   5725.0
VARS    9
dz=>dz    9  z  Layer thickness [m]
temp_2m=>temp_2m                                               0   t,y,x 2m temperature [K] 
U_wind_10m=>U_wind_10m                                         0   t,y,x U-component of 10m wind [m/s] 
V_wind_10m=>V_wind_10m                                         0   t,y,x V-component of 10m wind [m/s] 
prec_rate=>prec_rate                                           0   t,y,x precipitation rate [kg/m2s] 
ems_passive_gas=>ems_passive_gas                               0   t,y,x Emission intensity [massunits/s] passive_gas
cnc_passive_gas=>cnc_passive_gas                               9  t,z,y,x  Concentration in air passive_gas
dd_passive_gas=>dd_passive_gas                                 0   t,y,x Cocktail dry deposition passive_gas
wd_passive_gas=>wd_passive_gas                                 0   t,y,x Cocktail wet deposition passive_gas
ENDVARS
