from AdcircPy.Mesh import AdcircMesh
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid
from AdcircPy.Outputs import _Outputs
from AdcircPy.Outputs import _OutputStations
from AdcircPy.Outputs import _ElevationStations
# from AdcircPy.Outputs import _VelocityStations
from AdcircPy.Outputs import _HarmonicConstituentsStations
from AdcircPy.Outputs import _OutputSurface
from AdcircPy.Outputs import _Maxele
from AdcircPy.Outputs import _HarmonicConstituentsSurface

class Outputs(object):
  def __init__(self, path, fort14=None, fort15=None, datum='MSL', epsg=4326):
    self._path   = path
    self._fort14 = fort14
    self._fort15 = fort15
    self._datum  = datum
    self._epsg    = epsg

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

  def _read_ascii(self, **kwargs):
    """ """
    return _Outputs._read_ascii(self, **kwargs)

  def _parse_harmonic_constituent_output(self, f):
    """  """
    return _Outputs._parse_harmonic_constituent_output(self, f)

  def _parse_ascii_output(self, f):
    """  """
    return _Outputs._parse_ascii_output(self, f)
    
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

class _ASCII(object):
  def __init__(self, f, **kwargs):
    self._f = f
    self.__slice = 0.
  
  @property
  def _slice(self):
    return self.__slice

  @_slice.setter
  def _slice(self, __slice):
    self.__slice = __slice
    values = list()
    for i in range(self.x.size):
      values.append(float(f.readline().split()[1]))
    self.values = np.ma.masked_equal(values, -99999.)

class _OutputSurface(AdcircMesh, _netCDF4, _ASCII):
  def __init__(self, x, y, values, elements, **kwargs):
    AdcircMesh.__init__(self, x, y, values, elements, **kwargs)
    if 'Dataset' in kwargs.keys():
      _netCDF4.__init__(self, kwargs.pop('Dataset'), **kwargs)
    elif 'f' in kwargs.keys():
      _ASCII.__init__(self, kwargs.pop('f'))

  @staticmethod
  def from_ascii(path, fort14, **kwargs):
    return _OutputSurface._from_ascii(path, fort14, **kwargs)
    
  def make_animation(self, **kwargs):
    return _OutputSurface.make_animation(self, **kwargs)

class _OutputStations(dict):
  def __init__(self, time, **station_data):
    self.time = time
    dict.__init__(self, **station_data)

class ElevationStations(_OutputStations):
  def __init__(self, time, **station_data):
    _OutputStations.__init__(self, time, **station_data)
  
  @staticmethod
  def from_netcdf(path):
    return _ElevationStations._from_netcdf(path)

class HarmonicConstituentsStations(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)

  @staticmethod
  def from_fort51(path, fort14, fort15):
    return _HarmonicConstituentsStations._from_fort51(path, fort14, fort15)

class Maxele(_OutputSurface):
  def __init__(self, **kwargs):
    _OutputSurface.__init__(self, **kwargs)

  # @staticmethod
  # def from_ascii(path, fort14):
  #   _Maxele._from_ascii(path, fort14)
    
  @staticmethod
  def from_netcdf(path, fort14=None):
    return _Maxele._from_netcdf(path, fort14)

class VelocityStations(_OutputStations):
  def __init__(self, **kwargs):
    _OutputStations.__init__(self, **kwargs)

  def make_plot(self, station, **kwargs):
    return _VelocityStations._make_plot(self, station, **kwargs)

  @staticmethod
  def from_netcdf(path):
    return _VelocityStations._from_netcdf(path)



