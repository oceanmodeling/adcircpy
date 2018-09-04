from collections import OrderedDict
from AdcircPy.Tides import _TidalDB

class TidalDB(OrderedDict):
  def __init__(self):
    self._init_constituents()
    self._units='rad/sec'


  def _init_constituents(self):
    _TidalDB._init_constituents(self)