from AdcircPy.COOPS import _COOPS
from AdcircPy.COOPS import _HarmonicConstituents

class COOPS(dict, _COOPS._REST):
  def __init__(self, **kwargs):
    _COOPS._REST.__init__(self)
    dict.__init__(self, **kwargs)

  @staticmethod
  def from_station_list(stations, start_date, end_date):
    return _COOPS._from_station_list(stations, start_date, end_date)

class HarmonicConstituents(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)
  
  @staticmethod
  def from_station_list(stations):
    return _HarmonicConstituents._from_station_list(stations)

