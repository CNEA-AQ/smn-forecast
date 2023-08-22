#!/usr/bin/env python3

import zipfile
import csv
import datetime as dt
import sys
import numpy as np
import os.path
from dateutil.parser import parse
import pytz
import gzip
import unicodedata

from toolbox import MyTimeVars
csv.field_size_limit(sys.maxsize)

from toolbox import timeseries

year=2022
pollutants="CO PM10 PM25 O3 NOx NO NO2 SO2".split()
inp_dir="datafiles"
station_file="stations.csv"

def stations_from_stationfile(fname,sep=";"):
    with open(fname) as csvfile:
        stations = {}
        for l in csvfile:
            a=l.split(sep)
            code = a[0].strip()
            lon  = float(a[1])
            lat  = float(a[2])
            alt  = float(a[3])
            typ  = a[4].strip()
            are  = a[5].strip()
            stations[code] = timeseries.Station(code, code, lon, lat, alt, typ, are)
    return stations

stations = stations_from_stationfile(station_file,sep=";")

nst=len(stations)

codelist=sorted(stations.keys())
stlist=[]
stidx={}
for i,stcode in enumerate(codelist):
  stidx[stcode]=i
  stlist.append(stations[stcode])

#if mpirank == 0:
#    with  open("stDump.out", 'wt') as f:
#        for code in sorted(stations.keys()):
#            stations[code].toFile(f,";")

tstart=dt.datetime(year,1,1,0,0,0)
ndays=365
daylist=[tstart + dt.timedelta(days=d) for d in range(ndays)]
hrlist=[tstart + dt.timedelta(hours=h) for h in range(ndays*24)]
ntimes=len(hrlist)
timeidx={}
for i,t in enumerate(hrlist):
    timeidx[t]=i

duration=dt.timedelta(hours=1)

for ipol, pol in enumerate(pollutants):
    picklename="obs_%s.nc"%(pol,)
    
    if os.path.exists(picklename):
        print ("%s exists. Keeping it.."%(picklename,))
        continue
    valmatr = np.float32(np.nan) * np.empty((ntimes,nst),dtype=np.float32)

    polstr=pol  
    if pol=="PM25":
        polstr="PM2.5"
    elif pol=="NOX":
        polstr="NOX as NO2"

    for statid in list(stations.keys()):
        time_lut={}
        csvname="%s/%s_%d_%s.csv"%(inp_dir,statid,year,pol)
        if os.path.isfile(csvname):
            ist = stidx[statid]
            #print(istat,statid,pol)
            with open(csvname) as csvfile:
                reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
                for d in reader:
                   try: 
                       tstr=d["date"]
                       if tstr in time_lut:
                           time = time_lut[tstr]
                       else:
                           try:
                             time = parse(tstr).astimezone(pytz.utc).replace(tzinfo=None)
                           except:
                               time = None
                           time_lut[tstr] = time
                       if time == None:
                             print ("Wrong datetime row", tstr, d); continue
                       try:
                           it = timeidx[time]
                       except:
                           continue
                       try: 
                           val =float(d["conc"])
                       except:
                           val=np.float32(np.nan)
                           continue
                       if val >= 0: 
                         valmatr[it, ist] = val
                       else:
                         print ("Wrong value", d)
                   except AttributeError:
                           print (d); raise
                   except ValueError:
                           print (d); raise

    tsmatr = MyTimeVars.TsMatrix(hrlist,stlist,['val'],valmatr,["EUug/m3"])
    tsmatr.to_nc(picklename)

