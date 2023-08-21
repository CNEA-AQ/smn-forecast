#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Script for reading EEA's AQ-E-Reporting data. To use, download the zip files into a
# common directory, but don't extract. Then run the script on the directory for each species.
# 
# Timezones are handled.

# Inconsitencies in the 2013 dataset:
#
# AirQualityStationEoICode starts with STA- for some countries in the timeseries file but not in metadata
# station classification (background/traffic/etc) is missing altogether for Norway
# timezones missing some Italian stations

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
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpisize = comm.size
    mpirank = comm.Get_rank()
except:
    mpisize = 1
    mpirank = 0
    comm = None
    pass

from toolbox import timeseries

year=2022
pollutants="CO PM10 PM25 O3 NOX NO NO2 SO2".split()

def stations_from_stationfile(fname_set9):
    with open(fname_set9) as csvfile:
        bakstations = {}
        for l in csvfile:
            a=l.split()
            code = a[0]
            lon  = float(a[1])
            lat  = float(a[2])
            alt  = float(a[3])
            typ  = float(a[4])
            are  = float(a[5])
            bakstations[code] = timeseries.Station(code, code, lon, lat, alt, typ, are)
    return bakstations

stations = stations_from_stationfile("PanEuropean_metadata.csv.gz", backstations,delim='\t')

if mpirank == 0:
    with  open("stDump.out", 'wt') as f:
        for code in sorted(stations.keys()):
            stations[code].toFile(f,";")

#sys.exit()
tstart=dt.datetime(year,1,1,0,0,0)  ###End time for averaging
ndays=365
daylist=[ tstart + dt.timedelta(days=d) for d in range(ndays)]
hrlist=[ tstart + dt.timedelta(hours=h) for h in range(ndays*24)]
ntimes=len(hrlist)
timeidx={}
for i,t in enumerate(hrlist):
    timeidx[t]=i


nst=len(stations)

codelist=sorted(stations.keys())
stlist=[]
stidx={}
for i,stcode in enumerate(codelist):
  stidx[stcode]=i
  stlist.append(stations[stcode])
#stidx["SE0094A"] = stidx["SE0094R"] ##FIXME Hack


duration=dt.timedelta(hours=1)

for ipol, pol in enumerate(pollutants):
    if pol=="OX":
        continue
 #   break
    if ipol%mpisize != mpirank:
        continue
    picklename="obs_%s.nc"%(pol,)
    print("Rank %d doing %s"%(mpirank,picklename))

    if os.path.exists(picklename):
        print ("%s exists. Keeping it.."%(picklename,))
        continue

    polstr=pol  #POLUTANT STRING IN FILE
    if pol=="PM25":
        polstr="PM2.5"
    elif pol=="NOX":
        polstr="NOX as NO2"

    valmatr = np.float32(np.nan) * np.empty((ntimes,nst),dtype=np.float32)
    time_lut={}
    csvname="%d-tmp/%s_%d.gz"%(year,pol,year)
    with gzip.open(csvname, 'rt') as csvfile:
      line=0
      reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
      try:
        for d in reader:
           try: 
               line +=1 

               valid = d["Validity"]
               if not valid:
                   continue
               if  valid== '0' or   valid == 'Validity' or valid.startswith('-'):
                    continue
               if d["Verification"] == "0":
                    continue
               if d["AveragingTime"] == 'day' :
                    continue
               try:
                   code=d["AirQualityStationEoICode"]
                   if not code in stidx:
                       code = "S"+code##E164905 -> SE164905
                   ist = stidx[code] 
               except:
                   print ("Wrong station code", code, d)
                   continue
               
               tstr=d["DatetimeEnd"]
               if tstr in time_lut:
                   time = time_lut[tstr]
               else:
                   try:
                     time = parse(tstr).astimezone(pytz.utc).replace(tzinfo=None)
                   except:
                       time = None
                   time_lut[tstr] = time
               if time == None:
                     print ("Wrong DatetimeEnd", tstr, d)
                     continue
               
               tstr=d["DatetimeBegin"]
               if tstr in time_lut:
                   time_beg = time_lut[tstr]
               else:
                   try:
                     time_beg = parse(tstr).astimezone(pytz.utc).replace(tzinfo=None)
                   except:
                     time_beg = None
                     time_lut[tstr] = time_beg
               if time_beg == None:
                     print ("Wrong DatetimeBegin", tstr, d)
                     continue
               try:
                   it = timeidx[time]
               except:
                   #print "Wrong time",time
                   continue

               if d["AirPollutant"] != polstr:
                    print ("Wrong pollutant", d)
                    continue
               if time - time_beg != dt.timedelta(hours=1):
                    print ("Wrong time interval",time, time_beg, d)
                    continue

                    continue
               if d["AveragingTime"]  != "hour":
                    print ("Wrong time AveragingTime", d["AveragingTime"], d)
                    continue
               val = float(d["Concentration"])
               unit= d["UnitOfMeasurement"]

               if unit == "mg/m3":
                 val *= 1000
               elif not unit in  "µg/m3 ugNO2/m3".split():
                 if unit != "ng/m3" or pol != "NOX" or code[0:2] != "TR":  ##Turkish stations report ng/m3 for NOX
                     print ("Wrong unit", d)
                     continue
               
               if val >= 0: 
    #               print "Val ok", pol,code,time, val, it, ist 
                 valmatr[it, ist] = val
               elif val > -2:
                 valmatr[it, ist] = 0 ## pretend it is stil valid zero reading
                 ## Quiet this print ("Assumed zero value", d)
               elif polstr == 'CO' and val > -100:
                 valmatr[it, ist] = 0 ## pretend it is stil valid zero reading
               else:
                 print ("Wrong value", d)
           except AttributeError:
                   print (d)
                   raise
           except ValueError:
                   print (d)
                   raise

      except MemoryError:
          print ("Caught exception in ", csvfile)
          print ("line (approx)", line)
          print (d)
          pass
      except:
          print ("file:", csvname,  "line", line)
          raise

# #    #Bugfixes for values
      if pol=="CO": #Belgians still report milligrramms as micrograms
#         for code in "BA0032A BA0031A".split(): ## Switch units in Nov-dec
#             ist=stidx[code]
#             valmatr[valmatr[:,ist] > 50000,  ist] /= 1000
# 
         for code in "BA0037A BA0045A".split():
              #try:
                  ist=stidx[code]
                  valmatr[:,  ist] /= 1000. 
              #except:
              #    pass
# 
#         for code in "CY0002R".split(): #These guys report a lot of zeros
#              ist=stidx[code] 
#              valmatr[valmatr[:,ist] < 10,  ist] = np.nan
# 
# 
#         for code in "BA0049A BA0041A".split(): #These guys switched units
#              try:
#                  ist=stidx[code]
#                  valmatr[valmatr[:,ist] < 10,  ist] *= 1000. 
#              except:
#                  pass
# ##        
#         for code in "GB0839A GE0001A".split(): #Just nonsense
#           try:
#             ist=stidx[icode]  
#             valmatr[:,  ist] = np.nan

#          except:
#            pass
##
##
##        for code in "ES1573A MK0044A DEUB030 IT2008A ES1425A".split(): ## Few HUGE values
##             try:
##                 ist=stidx[code]
##                 iBAD = valmatr[:,  ist] > 10000.
##                 valmatr[iBAD,  ist] = np.nan
##             except:
##                 pass
##



#    if pol=="SO2":# 2017
#         idx = valmatr[:,  stidx["XK0002A"]]>4000 #9998 reading, clear outlier
#         valmatr[idx,  stidx["XK0002A"]] = np.nan
#    elif pol=="CO": #2017
#         idx = valmatr[:,  stidx["XK0002A"]]>1000 #station with wrong units
#         valmatr[:,  stidx["XK0002A"]] *= 1000
#         valmatr[idx,  stidx["XK0002A"]] = np.nan
#    el
#    if pol=="NO":
#         ist=stidx["MK0047A"] #2017 ## clear outlier
#         valmatr[:,  ist] = np.nan
#         ist=stidx["MK0048A"]  ## clear outlier
#         valmatr[valmatr[:, ist]>1000,  ist] = np.nan
#    if pol=="NO2":
#         ist=stidx["IT1486A"] # Some outliers
#         valmatr[valmatr[:,  ist]>1000,  ist] = np.nan
#    elif pol=="O3": #2018
#         for code in "BA0001A BA0029A BA0040A BA0042A BA0043A".split():
#             ist=stidx[code] #2018 or should they be divided by 1000?
#             valmatr[:,  ist] = np.nan
#    elif pol=="PM10":
#         ist=stidx["BA042A"]  # 2018 
#         valmatr[:,  ist] = np.nan
#
    # Bosnian stations report nonsense for 2018
#    for code in "BA0001G BA0029A BA0037A BA0040A BA0042A BA0043A".split():
#         ist=stidx[code] #2018 or should they be divided by 1000?
#         valmatr[:,  ist] = np.nan

    tsmatr = MyTimeVars.TsMatrix(hrlist,stlist,['val'],valmatr,["EUug/m3"])
    tsmatr.to_nc(picklename)


#
# Pickles to DB
#
#one_hour=dt.timedelta(hours=1)
#for ipol, pol in enumerate(pollutants):
#    if ipol%mpisize != mpirank:
#        continue
#    poltmp=pol
#    if pol=="OX":
#        poltmp="O3"
#    print "Getting obs_pickle_tsv_%s.gz"%(poltmp,)
#    with gzip.open("obs_pickle_tsv_%s.gz"%(poltmp,),'r') as zfile:
#           hrlist,stlist,valmatr = pickle.load(zfile)
#    if pol=="OX":
#        poltmp="NO2"
#        print "Getting obs_pickle_tsv_%s.gz"%(poltmp,)
#        with gzip.open("obs_pickle_tsv_%s.gz"%(poltmp,),'r') as zfile:
#               hrlist,stlist,valmatr1= pickle.load(zfile)
#        valmatr += valmatr1*48./46.
#
#
#    durations=[one_hour]*len(hrlist)
#    print "Dumping to DB"
#
#    for ist in range(len(stlist)):
#            series=timeseries.Timeseries(valmatr[:,ist], hrlist, durations, stations[stlist[ist]], pol, timeseries.alwaysValid)
#            db.put_series(series, db_timestamp, 'obs', 'europe') #Separate cases for pollutants
#db.release()
#



sys.exit()



class InputError(Exception):
    pass

input_encoding = 'utf-8'
micro = 'µ'.decode('utf-8').encode(input_encoding)
def get_unit_conv(unit):
    try:
        amount, volume = unit.split('/')
    except ValueError:
         
        print ('Failed to convert unit "%s"' % unit)
        raise
    if not volume == 'm3':
        raise ValueError('Cannot handle unit: %s' % unit)
    conv = {'ng':1e-12, 'ug':1e-9, 'mg':1e-6, 'g':1e-3, 'µg'.decode('utf-8'):1e-9}
    try:
        return conv[amount]
    except KeyError:
        return conv[amount.decode('utf-8')]


