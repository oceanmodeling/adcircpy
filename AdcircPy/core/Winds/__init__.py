from collections import defaultdict

from AdcircPy.core.Winds import _BestTrack

class BestTrack(object):
    _hurdat2_url="https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2017-050118.txt"
    
    def __init__(self, hurricane_id, start_date=None, end_date=None):
        self._init_hurdat2(hurricane_id, start_date, end_date)
    
    def _init_hurdat2(self, hurricane_id, start_date=None, end_date=None):
        _BestTrack._init_hurdat2(self, hurricane_id, start_date, end_date)

    def dump(self, path):
        _BestTrack._dump(self, path)
        

