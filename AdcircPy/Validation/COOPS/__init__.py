from AdcircPy.Validation.COOPS.TidalStations import TidalStations
from AdcircPy.Validation.COOPS import _HarmonicConstituents
    

class HarmonicConstituents(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)
  
  @staticmethod
  def from_station_list(stations):
    return _HarmonicConstituents._from_station_list(stations)

