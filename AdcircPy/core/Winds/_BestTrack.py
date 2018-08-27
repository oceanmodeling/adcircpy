import os
from datetime import datetime, timedelta
import numpy as np
_hurdat2 = os.path.dirname(os.path.abspath(__file__)) + '/hurdat2.txt'

def _init_hurdat2(self, hurricane_id, start_date=None, end_date=None):
  hurr = dict()
  with open(_hurdat2, 'r') as f:
    for line in f:
      line = line.split(',')
      if hurricane_id in line:
        line = [x.strip() for x in line]
        self._hurricane_id = hurricane_id
        self._name = line[1]
        self._hurricane_number = int(line[2])
        self._datetime = list()
        self._record_identifier = list()
        self._development_level = list()
        self._lat = list()
        self._lon = list()
        self._max_sustained_winds = list()
        self._min_pressure = list()
        self._wind_data = dict()
        self._wind_data['34']=dict()
        self._wind_data['50']=dict()
        self._wind_data['64']=dict()
        self._wind_data['34']['NEQ']=list()
        self._wind_data['34']['SEQ']=list()
        self._wind_data['34']['SWQ']=list()
        self._wind_data['34']['NWQ']=list()
        self._wind_data['50']['NEQ']=list()
        self._wind_data['50']['SEQ']=list()
        self._wind_data['50']['SWQ']=list()
        self._wind_data['50']['NWQ']=list()
        self._wind_data['64']['NEQ']=list()
        self._wind_data['64']['SEQ']=list()
        self._wind_data['64']['SWQ']=list()
        self._wind_data['64']['NWQ']=list()
        if start_date is None:
          _date = datetime.min
          start_date = datetime.min
        else:
          _date = start_date
        if end_date is None:
          end_date = datetime.max
        _len=0
        while _date <= end_date:
          data=f.readline().split(',')
          _len=len(data)
          if _len==4:
            break
          _new_date = datetime.strptime(data[0]+data[1],'%Y%m%d %H%M')
          if _new_date >= start_date and _new_date<=end_date:
            self._datetime.append(_new_date)
            self._record_identifier.append(data[2])
            self._development_level.append(data[3])
            self._lat.append(data[4])
            self._lon.append(data[5])
            self._max_sustained_winds.append(float(data[6]))
            self._min_pressure.append(float(data[7]))
            self._wind_data['34']['NEQ'].append(float(data[8]))
            self._wind_data['34']['SEQ'].append(float(data[9]))
            self._wind_data['34']['SWQ'].append(float(data[10]))
            self._wind_data['34']['NWQ'].append(float(data[11]))
            self._wind_data['50']['NEQ'].append(float(data[12]))
            self._wind_data['50']['SEQ'].append(float(data[13]))
            self._wind_data['50']['SWQ'].append(float(data[14]))
            self._wind_data['50']['NWQ'].append(float(data[15]))
            self._wind_data['64']['NEQ'].append(float(data[16]))
            self._wind_data['64']['SEQ'].append(float(data[17]))
            self._wind_data['64']['SWQ'].append(float(data[18]))
            self._wind_data['64']['NWQ'].append(float(data[19]))
  self._generate_best_track_data()


def _generate_best_track_data(self):
  self._best_track=list()
  def __isotach_data(isotach):
    #1 : basin
    string = "AL,"
    #2 : hurricane_no/hurricane_id
    string+= " {},".format(self._hurricane_id[2:4])
    #3 : record datetime
    string+= " {},".format(self._datetime[i].strftime('%Y%m%d%H'))
    #4 : landfall tag
    string+= "{:>3},".format(self._record_identifier[i])
    #5 : BestTrack tag
    string+= " BEST,"
    #6 : deltatime tag
    string+= "{:>4},".format(int((self._datetime[i]-_record_start_time).total_seconds()/3600))
    #7 : lat
    _lat = self._lat[i].split('.')
    _lat = _lat[0].strip()+_lat[1]
    string+= "{:>5},".format(_lat)
    #8 : lon
    _lon = self._lon[i].split('.')
    _lon = _lon[0].strip()+_lon[1]
    string+= "{:>6},".format(_lon)
    #9 : max sustained_winds 
    string+= "{:>4},".format(int(self._max_sustained_winds[i]))
    #10 : min pressure
    string+= "{:>5},".format(int(self._min_pressure[i]))
    #11 : development level
    string+= "{:>3},".format(self._development_level[i])
    #12 : isotach
    string+= "{:>4},".format(isotach)
    #13 : quadrant
    string+= "{:>4},".format('')
    #14 : radius of isotach NEQ
    string+= "{:>5},".format(int(self._wind_data[isotach]['NEQ'][i]))
    #15 : radius of isotach SEQ
    string+= "{:>5},".format(int(self._wind_data[isotach]['SEQ'][i]))
    #16 : radius of isotach SWQ
    string+= "{:>5},".format(int(self._wind_data[isotach]['SWQ'][i]))
    #17 : radius of isotach NWQ
    string+= "{:>5},".format(int(self._wind_data[isotach]['NWQ'][i]))
    #18 : background pressure
    string+= "{:>5},".format('1013')
    #19 : radius of last closed isobar RRP, unsued
    string+= "{:>5},".format('')
    #20 : Rmax, unused
    string+= "{:>5},".format('')
    #21 : gusts, unused
    string+= "{:>4},".format('')
    #22 : eye diamer, unused
    string+= "{:>4},".format('')
    #23 : subregion, unsused
    string+= "{:>4},".format('')
    #24 : maxseas, unsused
    string+= "{:>4},".format('')
    #25 : initials, unused
    string+= "{:>4},".format('')
    #26 : direction, **aswip**
    string+= "{:>3},".format('')
    #27 : speed, **aswip**
    string+= "{:>4},".format('')
    #28 : stormname
    string+= "{:^12},".format(self._name)
    #29 : record number, **aswip**
    string+= "{:>4},".format('')
    #30 : number of isotachs reported
    string+= "{:>5},".format(_cnt)
    #31 : use NEQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #32 : use SEQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #33 : use SWQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #34 : use NWQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    self._best_track.append(string)
  
  for i in range(len(self._datetime)):
    for j, isotach in  enumerate(['34', '50', '64']):
      _isotach=list()
      for _isot in list(self._wind_data.keys()):
        for quadrant in list(self._wind_data[_isot].keys()):
          _isotach.append(self._wind_data[_isot][quadrant][i])
      _34_iso = _isotach[0:4]
      _50_iso = _isotach[4:8]
      _64_iso = _isotach[8:12]
      _cnt=0
      if any(_34_iso) > 0:
        _cnt+=1
      if any(_50_iso) > 0:
        _cnt+=1
      if any(_64_iso) > 0:
        _cnt+=1
      if _cnt==0:
        continue
      if any(_34_iso)>0 and isotach=='34':
        try: _record_start_time
        except: _record_start_time = self._datetime[i]
        __isotach_data(isotach)
      if any(_50_iso)>0 and isotach=='50':
        try: _record_start_time
        except: _record_start_time = self._datetime[i]
        __isotach_data(isotach)
      if any(_64_iso)>0 and isotach=='64':
        try: _record_start_time
        except: _record_start_time = self._datetime[i]
        __isotach_data(isotach)

def _dump(self, path):
  with open(path,'w') as f:
    for line in self._best_track:
      f.write(line)
      if line != self._best_track[-1]:
        f.write('\n')

def _printf(self):
  for line in self._best_track:
    print(line)
