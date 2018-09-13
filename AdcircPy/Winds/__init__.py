from collections import defaultdict
from AdcircPy.Winds import _BestTrack

class BestTrack(object):
  _hurdat2_url="https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2017-050118.txt"
  def __init__(self, hurricane_id, start_date=None, end_date=None):
    self._init_hurdat2(hurricane_id, start_date, end_date)
    
  @property
  def start_time(self):
    return self._datetime[0]

  @property
  def end_time(self):
    return self._datetime[-1]
 
  @property
  def name(self):
    return self._name

  def dump(self, path):
    _BestTrack._dump(self, path)
    
  def printf(self):
    _BestTrack._printf(self)

  def _init_hurdat2(self, hurricane_id, start_date=None, end_date=None):
    _BestTrack._init_hurdat2(self, hurricane_id, start_date, end_date)

  def _generate_best_track_data(self):
    _BestTrack._generate_best_track_data(self)
