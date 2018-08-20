from AdcircPy.core.Validation import COOPS
from AdcircPy.core.Validation.COOPS import _HarmonicConstituents

class HarmonicConstituents(COOPS):
  def __init__(self, **kwargs):
    COOPS.__init__(self, **kwargs)
  
  @staticmethod
  def from_stations_list(stations):
    return _HarmonicConstituents._from_stations_list(stations)
