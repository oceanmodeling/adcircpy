from AdcircPy.Mesh import AdcircMesh
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid
from AdcircPy.Outputs import _Outputs
# from AdcircPy.Outputs import __Stations
from AdcircPy.Outputs import _ElevationStations
# from AdcircPy.Outputs import _VelocityStations
from AdcircPy.Outputs import _HarmonicConstituentsStations
from AdcircPy.Outputs import __OutputSurface
from AdcircPy.Outputs import _Maxele
from AdcircPy.Outputs import _HarmonicConstituentsSurface

class Outputs(object):
  def __init__(self, path, fort14=None, fort15=None, fort13=None, datum=None, epsg=None):
    self._path   = path
    self._fort14 = fort14
    self._fort15 = fort15
    self._fort13 = fort13
    self.datum   = datum
    self.epsg    = epsg

  @staticmethod
  def read_outputs(path, **kwargs):
    """ """
    return _Outputs.read_outputs(path, **kwargs)

  def _get_netcdf_output(self):
    """ """
    return _Outputs._get_output(self)

  def _open_file(self):
    """ """
    return _Outputs._open_file(self)

  def _read_netcdf(self):
    """ """
    return _Outputs._read_netcdf(self)

  def _read_ascii_type(self):
    """ """
    return _Outptuts._read_ascii_type(self)
    
class _netCDF4(object):
  def __init__(self, Dataset, var, **kwargs):
    self.Dataset = Dataset
    self.__slice = 0
    self.__var   = var

  @property
  def _slice(self):
    return self.__slice

  @_slice.setter
  def _slice(self, __slice):
    self.__slice = __slice
    self.values = self.Dataset[self._var][self._slice,:,:]

  @property
  def _var(self):
    return self.__var

  @_var.setter
  def _var(self, name):
    if name not in self.Dataset.keys():
      raise Exception('Invalid netCDF variable. Options are {}'.format(self.Dataset.keys()))
    self.__var = name
    self.values = self.Dataset[self._var][self._slice,:,:]

class _OutputSurface(UnstructuredGrid, _netCDF4):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, **kwargs)
    _netCDF4.__init__(self, kwargs.pop('Dataset', None), **kwargs)

  def make_animation(self, **kwargs):
    return __OutputSurface.make_animation(self, **kwargs)

class _Stations(dict, _netCDF4):
  def __init__(self, **kwargs):
    _netCDF4.__init__(self, Dataset = kwargs.pop('Dataset', None),
                            var = kwargs.pop('var', None), **kwargs)
    dict.__init__(self, **kwargs)

class ElevationStations(_Stations):
  def __init__(self, **kwargs):
    _Stations.__init__(self, **kwargs)

  def make_plot(self, station, **kwargs):
    return _ElevationStations._make_plot(self, station, **kwargs)

  @staticmethod
  def from_netcdf(path):
    return _ElevationStations._from_netcdf(path)

class HarmonicConstituentsStations(_Stations):
  def __init__(self, **kwargs):
    _Stations.__init__(self, **kwargs)

  @staticmethod
  def from_fort51(path, fort14, fort15):
    return _HarmonicConstituentsStations._from_fort51(path, fort14, fort15)

class Maxele(_OutputSurface):
  def __init__(self, **kwargs):
    _OutputSurface.__init__(self, **kwargs)

  @staticmethod
  def from_ascii(path, fort14):
    _Maxele._from_ascii(path, fort14)
    
  @staticmethod
  def from_netcdf(path, fort14=None):
    return _Maxele._from_netcdf(path, fort14)

class VelocityStations(_Stations):
  def __init__(self, **kwargs):
    _Stations.__init__(self, **kwargs)

  def make_plot(self, station, **kwargs):
    return _VelocityStations._make_plot(self, station, **kwargs)

  @staticmethod
  def from_netcdf(path):
    return _VelocityStations._from_netcdf(path)



