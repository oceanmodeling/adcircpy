import os
from datetime import datetime
import numpy as np

class HURDAT2(object):
  """
   hurdat2_url="https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2017-050118.txt"
  """
  def __init__(self, hurdat2_id, start_date=None, end_date=None):
    self._start_date = start_date
    self._end_date   = end_date
    self._hurdat2_id = hurdat2_id
    self._init_vars()
    self._parse_hurdat2()
    self._init_dates()

  @property
  def start_date(self):
    return self._datetime[0]

  @property
  def end_date(self):
    return self._datetime[-1]

  @property
  def datetime(self):
    return self._datetime

  @property
  def record_identifier(self):
    return self._record_identifier

  @property
  def development_level(self):
    return self._development_level

  @property
  def latitude(self):
    return self._latitude
  @property
  def longitude(self):
    return self._longitude

  @property
  def max_sustained_wind(self):
    return self._max_sustained_wind
  
  @property
  def min_pressure(self):
    return self._min_pressure

  @property
  def wind_data(self):
    return self._wind_data

  @property
  def hurdat2_id(self):
    return self._hurdat2_id
  
  def get_ATCF(self):
    pass

  def _init_vars(self):
    self._datetime=list()
    self._record_identifier=list()
    self._development_level=list()
    self._latitude=list()
    self._longitude=list()
    self._max_sustained_winds=list()
    self._min_pressure=list()
    self._wind_data=dict()
    self._wind_data['34']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}
    self._wind_data['50']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}
    self._wind_data['64']={'NEQ':list(), 'SEQ':list(), 'SWQ':list(), 'NWQ':list()}

  def __remove_invalid_isotach_data(self):
    indexes=list()
    for i,_ in enumerate(self.datetime):
      isotach_data=list()
      for isotach in self.wind_data.keys():
        for quadrant in self.wind_data[isotach].keys():
          isotach_data.append(self.wind_data[isotach][quadrant][i])
      if all(isotach == 0 for isotach in isotach_data):
        indexes.append(i)
      elif all(isotach == -999 for isotach in isotach_data):
        indexes.append(i)
    for index in sorted(indexes, reverse=True):
      del self._datetime[index]
      del self._record_identifier[index]
      del self._development_level[index]
      del self._latitude[index]
      del self._longitude[index]
      del self._max_sustained_winds[index]
      del self._min_pressure[index]
      for isotach in self.wind_data.keys():
        for quadrant in self.wind_data[isotach].keys():
          del self.wind_data[isotach][quadrant][index]

  def _parse_hurdat2(self):
    found=False
    with open(os.path.dirname(os.path.abspath(__file__))+'/hurdat2.txt', 'r') as f:
      for line in f:
        line = line.split(',')
        if self._hurdat2_id in line:
          line = [x.strip() for x in line]
          self._event_name = line[1]
          self._event_number = int(line[2])
          while True:
            line=f.readline().split(',')
            if len(line)==4:
              found=True
              break
            self._datetime.append(datetime.strptime(line[0]+line[1],'%Y%m%d %H%M'))
            self._record_identifier.append(line[2])
            self._development_level.append(line[3])
            self._latitude.append(line[4])
            self._longitude.append(line[5])
            self._max_sustained_winds.append(float(line[6]))
            self._min_pressure.append(float(line[7]))
            self._wind_data['34']['NEQ'].append(float(line[8]))
            self._wind_data['34']['SEQ'].append(float(line[9]))
            self._wind_data['34']['SWQ'].append(float(line[10]))
            self._wind_data['34']['NWQ'].append(float(line[11]))
            self._wind_data['50']['NEQ'].append(float(line[12]))
            self._wind_data['50']['SEQ'].append(float(line[13]))
            self._wind_data['50']['SWQ'].append(float(line[14]))
            self._wind_data['50']['NWQ'].append(float(line[15]))
            self._wind_data['64']['NEQ'].append(float(line[16]))
            self._wind_data['64']['SEQ'].append(float(line[17]))
            self._wind_data['64']['SWQ'].append(float(line[18]))
            self._wind_data['64']['NWQ'].append(float(line[19]))
    if found==False:
      raise Exception('The storm id provided ({}) was not found in HURDAT2 database!'.format(self._hurdat2_id))
    self.__remove_invalid_isotach_data()
  
  def __init_start_date(self):
    if self._start_date is not None:
      if isinstance(self._start_date, datetime)==False:
        raise Exception("start_date must be a datetime.datetime instance.")
  
  def __init_end_date(self):
    if self._end_date is not None:
      if isinstance(self._end_date, datetime)==False:
        raise Exception("end_date must be a datetime.datetime instance.")

  def __init_datetime_indexes(self):
    if self._start_date is not None:
      _diff = [np.abs((x - self._start_date).total_seconds()) for x in self._datetime]
      start_index = _diff.index(min(_diff))
    else:
      start_index=0
    if self._end_date is not None:
      _diff = [np.abs((self._end_date - x).total_seconds()) for x in self._datetime]
      end_index = _diff.index(min(_diff))
    else:
      end_index=-1
    self._datetime = self._datetime[start_index:end_index]
    self._record_identifier = self._record_identifier[start_index:end_index]
    self._development_level = self._development_level[start_index:end_index]
    self._latitude = self._latitude[start_index:end_index]
    self._longitude = self._longitude[start_index:end_index]
    self._max_sustained_winds = self._max_sustained_winds[start_index:end_index]
    self._min_pressure = self._min_pressure[start_index:end_index]
    for isotach in self._wind_data.keys():
      for quadrant in self._wind_data[isotach].keys():
        self._wind_data[isotach][quadrant] = self._wind_data[isotach][quadrant][start_index:end_index]

  def _init_dates(self):
    self.__init_start_date()
    self.__init_end_date()
    self.__init_datetime_indexes()
    if self._start_date is None:
      self._start_date=self._datetime[0]
    if self._end_date is None:
      self._end_date=self._datetime[-1]
    if self._end_date<=self._start_date:
      raise Exception("start_date must be previous to end_date.")
    if self._start_date<self._datetime[0] or self._start_date>self._datetime[-1]:
      raise Exception("start_date provided is out of range with the record.\n"+
                      "Valid date ranges are from {} to {}".format(self._datetime[0], self._datetime[-1]))
    if self._end_date<self._datetime[0] or self._end_date>self._datetime[-1]:
      raise Exception("end_date provided is out of range with the record.\n"+
                      "Valid date ranges are from {} to {}".format(self._datetime[0], self._datetime[-1]))


