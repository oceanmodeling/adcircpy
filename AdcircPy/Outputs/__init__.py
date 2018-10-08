from AdcircPy.Model import UnstructuredMesh
from AdcircPy.Outputs import _OutputFactory as _OutputFactory_
from AdcircPy.Outputs import _ScalarSurfaceExtrema as _ScalarSurfaceExtrema_
from AdcircPy.Outputs import _ScalarSurfaceTimeseries as _ScalarSurfaceTimeseries_
from AdcircPy.Outputs import Maxele as _Maxele

class _OutputFactory(object):
  """
  Private class called by AdcircPy.read_output() that returns
  the appropriate subclass belonging to the given output file.
  Supports ASCII and NetCDF outputs.
  Fortran binary outputs are not supported.
  """

  def __new__(cls, path, fort14=None, fort15=None, datum='LMSL', epsg=None, datum_grid=None):
    return _OutputFactory_.__new__(cls, path, fort14, fort15, datum, epsg, datum_grid)

  @staticmethod
  def is_ncfile(path):
    """ Tests if input file is NetCDF format. """
    return _OutputFactory_.is_ncfile(path)

  @staticmethod
  def netcdf_factory(path, fort14=None, fort15=None, datum='LMSL', epsg=None, datum_grid=None):
    """ Will determine which type of netcdf output and returns an instance of the appropritate class. """
    return _OutputFactory_.netcdf_factory(path, fort14, fort15, datum, epsg, datum_grid)

  @staticmethod
  def ascii_factory(path, fort14=None, fort15=None, datum='LMSL', epsg=None, datum_grid=None):
    """ Will determine which type of ascii output and returns an instance of the appropritate class. """
    return _OutputFactory_.ascii_factory(path, fort14, fort15, datum, epsg, datum_grid)

class _ScalarSurfaceExtrema(UnstructuredMesh):
  """
  Private subclass that represents a general Output Surface Extrema, which is instantiated
  by the _OutputFactory class.
  """
  def __init__(self, x, y, elements, values, times, **kwargs):
    self._times = times
    UnstructuredMesh.__init__(self, x, y, elements, values, **kwargs)
      
  @classmethod
  def from_file(cls, path, fort14, **kwargs):
    return _ScalarSurfaceExtrema_.from_file(cls, path, fort14, **kwargs)

  def make_plot(self, **kwargs):
    _ScalarSurfaceExtrema_.make_plot(self, **kwargs)

class _ScalarSurfaceTimeseries(UnstructuredMesh):
  """
  Private subclass that 
  """
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
    self._f = f
    self._slice = 0.

  @staticmethod
  def from_file(path, fort14, **kwargs):
    return _ScalarSurfaceTimeseries_._from_ascii(path, fort14, **kwargs)

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

class _VectorSurfaceExtrema(UnstructuredMesh):
  """
  Private subclass that 
  """
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
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
    
class _VectorSurfaceTimeseries(UnstructuredMesh):
  pass

class _NetCDFScalarSurfaceTimeseries(_ScalarSurfaceTimeseries):
  """ Contains reusable methods for the NetCDF scalar surface timeseries outputs """
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
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

class _NetCDFVectorSurfaceTimeseries(_VectorSurfaceTimeseries):
  """ Contains reusable methods for the NetCDF vector surface timeseries outputs """
  pass

class Maxele(UnstructuredMesh):
  def __init__(self, x, y, elements, values, times, **kwargs):
    self._times = times
    UnstructuredMesh.__init__(self, x, y, elements, values, **kwargs)
    
  @classmethod
  def from_ascii(cls, path, fort14):
    _Maxele.from_ascii(cls, path, fort14)
    
  @classmethod
  def from_netcdf(cls, path, fort14=None, fort15=None, datum='MSL', epsg=None, datum_grid=None):
    return _Maxele.from_netcdf(cls, path, fort14, fort15, datum, epsg, datum_grid)
