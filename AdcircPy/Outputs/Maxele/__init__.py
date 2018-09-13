from AdcircPy.Outputs import _OutputSurface
from AdcircPy.Outputs.Maxele import _Maxele

class Maxele(_OutputSurface):
  def __init__(self, **kwargs):
    _OutputSurface.__init__(self, **kwargs)

  @staticmethod
  def from_ascii(path, fort14):
    _Maxele._from_ascii(path, fort14)
    
  @staticmethod
  def from_netcdf(path, fort14=None):
    return _Maxele._from_netcdf(path, fort14)