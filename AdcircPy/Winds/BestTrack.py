import os
from datetime import datetime
from collections import defaultdict

class BestTrack(object):
  _hurdat2_url="https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2017-050118.txt"
  def __init__(self, hurdat2_id, start_date=None, end_date=None):
    self.start_date = start_date
    self.end_date   = end_date
    self.hurdat2_id = hurdat2_id
    self._init_vars()
    self._init_hurdat2()
    self._cleanup_record()
    self._init_dates()
    self._init_datetime_indexes()
    self._generate_best_track()

  def dump(self, path):
    with open(path,'w') as f:
      for line in self.best_track:
        f.write(line)
        if line != self.best_track[-1]:
          f.write('\n')

  def printf(self):
    for line in self.best_track:
      print(line)

  def _init_vars(self):
    self.datetime = list()
    self.record_identifier = list()
    self.development_level = list()
    self.lat = list()
    self.lon = list()
    self.max_sustained_winds = list()
    self.min_pressure = list()
    self.wind_data = dict()
    self.wind_data['34']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}
    self.wind_data['50']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}
    self.wind_data['64']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}
  
  def __process_line(self):
    self.datetime.append(datetime.strptime(self.data[0]+self.data[1],'%Y%m%d %H%M'))
    self.record_identifier.append(self.data[2])
    self.development_level.append(self.data[3])
    self.lat.append(self.data[4])
    self.lon.append(self.data[5])
    self.max_sustained_winds.append(float(self.data[6]))
    self.min_pressure.append(float(self.data[7]))
    self.wind_data['34']['NEQ'].append(float(self.data[8]))
    self.wind_data['34']['SEQ'].append(float(self.data[9]))
    self.wind_data['34']['SWQ'].append(float(self.data[10]))
    self.wind_data['34']['NWQ'].append(float(self.data[11]))
    self.wind_data['50']['NEQ'].append(float(self.data[12]))
    self.wind_data['50']['SEQ'].append(float(self.data[13]))
    self.wind_data['50']['SWQ'].append(float(self.data[14]))
    self.wind_data['50']['NWQ'].append(float(self.data[15]))
    self.wind_data['64']['NEQ'].append(float(self.data[16]))
    self.wind_data['64']['SEQ'].append(float(self.data[17]))
    self.wind_data['64']['SWQ'].append(float(self.data[18]))
    self.wind_data['64']['NWQ'].append(float(self.data[19]))

  def _init_hurdat2(self):
    hurr = dict()
    with open(os.path.dirname(os.path.abspath(__file__))+'/hurdat2.txt', 'r') as self.f:
      for line in self.f:
        line = line.split(',')
        if self.hurdat2_id in line:
          line = [x.strip() for x in line]
          self.event_name = line[1]
          self.event_number = int(line[2])
          while True:
            self.data=self.f.readline().split(',')
            if len(self.data)==4:
              break
            self.__process_line()

  def _cleanup_record(self):
    """
    ADCIRC requires 6-hourly forcing only, but HURDAT2 contains some other
    values in between. This function removes those values so that the delta
    time between record entries remains 6 hourly.
    """
    idxs_to_delete = list()
    for i, dates in enumerate(self.datetime):
      if self.datetime[i].hour not in [0, 6, 8, 12] or self.datetime[i].minute!=0:
        idxs_to_delete.append(i)
    for index in sorted(idxs_to_delete, reverse=True):
      del self.datetime[index]
      del self.record_identifier[index]
      del self.development_level[index]
      del self.lat[index]
      del self.lon[index]
      del self.max_sustained_winds[index]
      del self.min_pressure[index]
      for isotach in self.wind_data.keys():
        for quadrant in self.wind_data[isotach]:
          del self.wind_data[isotach][quadrant][index]
  
  def _init_dates(self):
    if self.start_date is not None or self.end_date is not None:
      if isinstance(self.start_date, datetime)==False:
        raise Exception("start_date must be a datetime.datetime instance.")
      if isinstance(self.end_date, datetime)==False:
        raise Exception("end_date must be a datetime.datetime instance.")
      if self.end_date>=self.start_date:
        raise Exception("start_date must be previous to end_date.")
      if self.datetime[0]>self.end_date or self.end_date>self.datetime[-1]:
        raise Exception("end_date provided is out of range with the record")
      if self.datetime[0]>self.start_date or self.start_date>self.datetime[-1]:
        raise Exception("start_date provided is out of range with the record")
    else:
      self.start_date = self.datetime[0]
      self.end_date = self.datetime[-1]

  def _init_datetime_indexes(self):
    _diff = [(x - self.start_date).total_seconds() for x in self.datetime]
    self.start_index = _diff.index(min(_diff))
    _diff = [(self.end_date - x).total_seconds() for x in self.datetime]
    self.end_index = _diff.index(min(_diff))


  def _generate_best_track(self):
    self.best_track=list()
    for self.i in range(self.start_index, self.end_index+1):
      for self._isotach in  self.wind_data.keys():
        _isotach=list()
        for _isot in list(self.wind_data.keys()):
          for quadrant in list(self.wind_data[_isot].keys()):
            _isotach.append(self.wind_data[_isot][quadrant][self.i])
        _34_iso = _isotach[0:4]
        _50_iso = _isotach[4:8]
        _64_iso = _isotach[8:12]
        self._cnt=0
        if any(_34_iso) > 0:
          self._cnt+=1
        if any(_50_iso) > 0:
          self._cnt+=1
        if any(_64_iso) > 0:
          self._cnt+=1
        if self._cnt==0:
          continue
        if any(_34_iso)>0 and self._isotach=='34':
          self.__append_isotach_data()
        if any(_50_iso)>0 and self._isotach=='50':
          self.__append_isotach_data()
        if any(_64_iso)>0 and self._isotach=='64':
          self.__append_isotach_data()

  def __append_isotach_data(self):
    #1 : basin
    string = "AL,"
    #2 : hurricane_no/hurricane_id
    string+= " {},".format(self.hurdat2_id[2:4])
    #3 : record datetime
    string+= " {},".format(self.datetime[self.i].strftime('%Y%m%d%H'))
    #4 : landfall tag
    string+= "{:>3},".format(self.record_identifier[self.i])
    #5 : BestTrack tag
    string+= " BEST,"
    #6 : deltatime tag
    string+= "{:>4},".format(int((self.datetime[self.i]-self.datetime[0]).total_seconds()/3600))
    #7 : lat
    _lat = self.lat[self.i].split('.')
    _lat = _lat[0].strip()+_lat[1]
    string+= "{:>5},".format(_lat)
    #8 : lon
    _lon = self.lon[self.i].split('.')
    _lon = _lon[0].strip()+_lon[1]
    string+= "{:>6},".format(_lon)
    #9 : max sustained_winds 
    string+= "{:>4},".format(int(self.max_sustained_winds[self.i]))
    #10 : min pressure
    string+= "{:>5},".format(int(self.min_pressure[self.i]))
    #11 : development level
    string+= "{:>3},".format(self.development_level[self.i])
    #12 : isotach
    string+= "{:>4},".format(self._isotach)
    #13 : quadrant
    string+= "{:>4},".format('')
    #14 : radius of isotach NEQ
    string+= "{:>5},".format(int(self.wind_data[self._isotach]['NEQ'][self.i]))
    #15 : radius of isotach SEQ
    string+= "{:>5},".format(int(self.wind_data[self._isotach]['SEQ'][self.i]))
    #16 : radius of isotach SWQ
    string+= "{:>5},".format(int(self.wind_data[self._isotach]['SWQ'][self.i]))
    #17 : radius of isotach NWQ
    string+= "{:>5},".format(int(self.wind_data[self._isotach]['NWQ'][self.i]))
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
    string+= "{:^12},".format(self.event_name)
    #29 : record number, **aswip**
    string+= "{:>4},".format('')
    #30 : number of isotachs reported
    string+= "{:>5},".format(self._cnt)
    #31 : use NEQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #32 : use SEQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #33 : use SWQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    #34 : use NWQ isotach flag, **aswip**
    string+= "{:>2},".format('')
    self.best_track.append(string)
