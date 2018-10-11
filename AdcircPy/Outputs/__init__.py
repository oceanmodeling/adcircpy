from AdcircPy.Outputs._OutputFactory import _OutputFactory
from AdcircPy.Outputs.Maxele import Maxele


# class _VectorSurfaceExtrema(UnstructuredMesh):
#   """
#   Private subclass that 
#   """
#   def __init__(self, x, y, values, elements, **kwargs):
#     UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
#     self._f = f
#     self._slice = 0.

#   # def __new__(cls):
#   #   return cls.__init__(cls, cls.x, cls.y, cls.values, cls.elements, cls.f)
  
#   @staticmethod
#   def from_ascii(path, fort14, **kwargs):
#     return _OutputSurface_._from_ascii(path, fort14, **kwargs)

#   @property
#   def slice(self):
#     return self._slice

#   @slice.setter
#   def slice(self, _slice):
#     self._slice = _slice
#     values = list()
#     for i in range(self.x.size):
#       values.append(float(f.readline().split()[1]))
#     self.values = np.ma.masked_equal(values, -99999.)
    
# class _VectorSurfaceTimeseries(UnstructuredMesh):
#   pass

# class _NetCDFScalarSurfaceTimeseries(_ScalarSurfaceTimeseries):
#   """ Contains reusable methods for the NetCDF scalar surface timeseries outputs """
#   def __init__(self, x, y, values, elements, **kwargs):
#     UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
#     self.Dataset = Dataset
#     self.__slice = 0
#     self.__var   = var

#   @property
#   def _slice(self):
#     return self.__slice

#   @_slice.setter
#   def _slice(self, __slice):
#     self.__slice = __slice
#     self.values = self.Dataset[self._var][self._slice,:,:]

#   @property
#   def _var(self):
#     return self.__var

#   @_var.setter
#   def _var(self, name):
#     if name not in self.Dataset.keys():
#       raise Exception('Invalid netCDF variable. Options are {}'.format(self.Dataset.keys()))
#     self.__var = name
#     self.values = self.Dataset[self._var][self._slice,:,:]

#   @staticmethod
#   def from_ascii(path, fort14, **kwargs):
#     return _OutputSurface_._from_ascii(path, fort14, **kwargs)
    
#   def make_animation(self, **kwargs):
#     return _OutputSurface_._make_animation(self, **kwargs)

# class _NetCDFVectorSurfaceTimeseries(_VectorSurfaceTimeseries):
#   """ Contains reusable methods for the NetCDF vector surface timeseries outputs """
#   pass

