class ScalarSurfaceTimeseries(UnstructuredMesh):
  """
  Private subclass that 
  """
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredMesh.__init__(self, x, y, values, elements, **kwargs)
    # self._f = f
    # self._slice = 0.

  # @classmethod
  # def from_ascii(cls, path, fort14, **kwargs):
  #   return _ScalarSurfaceTimeseries_._from_ascii(path, fort14, **kwargs)

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
