from AdcircPy.Mesh import UnstructuredGrid
from AdcircPy.Outputs import _OutputSurfaceExtrema as _OutputSurfaceExtrema_
from AdcircPy.Outputs import _OutputFactory as _OutputFactory_

class _OutputFactory(object):
  def __init__(self, path, fort14=None, fort15=None, datum='MSL', epsg=None):
    self._path=path
    self._fort14=fort14
    self._fort15=fort15
    self._datum=datum
    self._epsg=epsg

  def __new__(cls, path, fort14=None, fort15=None, datum='MSL', epsg=None):
    cls.__init__(cls, path, fort14, fort15, datum, epsg)
    return cls._read_file(cls)
  
  def _read_file(self):
    return _OutputFactory_._read_file(self)
  
  def _read_netcdf(self):
    """ """
    return _OutputFactory_._read_netcdf(self)

  def _read_ascii(self, **kwargs):
    """ """
    return _OutputFactory_._read_ascii(self, **kwargs)

class _OutputSurfaceExtrema(UnstructuredGrid):
  def __init__(self, x, y, elements, values, times, **kwargs):
    self._times = times
    UnstructuredGrid.__init__(self, x, y, elements, values, **kwargs)
  
  @staticmethod
  def from_ascii(path, fort14, **kwargs):
    return _OutputSurfaceExtrema_._from_ascii(path, fort14, **kwargs)

  def make_plot(self, **kwargs):
    _OutputSurfaceExtrema_._make_plot(self, **kwargs)

class _AsciiOutputSurfaceTimeseries(UnstructuredGrid):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, x, y, values, elements, **kwargs)
    self._f = f
    self._slice = 0.

  # def __new__(cls):
  #   return cls.__init__(cls, cls.x, cls.y, cls.values, cls.elements, cls.f)
  
  @staticmethod
  def from_ascii(path, fort14, **kwargs):
    return _OutputSurface_._from_ascii(path, fort14, **kwargs)

  @property
  def slice(self):
    return self._slice

  @slice.setter
  def slice(self, _slice):
    self._slice = _slice
    values = list()
    for i in range(self.x.size):
      values.append(float(f.readline().split()[1]))
    self.values = np.ma.masked_equal(values, -99999.)
    
class _NetCDFOutputSurface(UnstructuredGrid):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, x, y, values, elements, **kwargs)
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

  @staticmethod
  def from_ascii(path, fort14, **kwargs):
    return _OutputSurface_._from_ascii(path, fort14, **kwargs)
    
  def make_animation(self, **kwargs):
    return _OutputSurface_._make_animation(self, **kwargs)

class _NetCDFOutputSurface(UnstructuredGrid):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, x, y, values, elements, **kwargs)
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

  @staticmethod
  def from_ascii(path, fort14, **kwargs):
    return _OutputSurface_._from_ascii(path, fort14, **kwargs)
    
  def make_animation(self, **kwargs):
    return _OutputSurface_._make_animation(self, **kwargs)

