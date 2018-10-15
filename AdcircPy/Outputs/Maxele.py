from netCDF4 import Dataset
from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs.ScalarSurfaceExtrema import ScalarSurfaceExtrema

class Maxele(ScalarSurfaceExtrema):
  def __init__(self, x, y, elements, values, times, **kwargs):
    super(Maxele, self).__init__(x, y, elements, values, times, **kwargs)

  @classmethod
  def from_netcdf(cls, path, fort14=None, datum='MSL', epsg=None, datum_grid=None):
    nc = Dataset(path)
    if 'zeta_max' not in nc.variables.keys():
      raise Exception('Not a maxele file!')
    if fort14 is not None:
      if isinstance(fort14, AdcircMesh)==False:
        fort14 = AdcircMesh.from_fort14(fort14)
    return cls(nc['x'][:],
               nc['y'][:],
               nc['element'][:]-1,
               nc['zeta_max'][:],
               nc['time_of_zeta_max'],
               datum=datum,
               # **fort14.get_boundary_dictionary() # Need to add the additional infor provided by fort.14
               datum_grid=datum_grid)

  @classmethod
  def from_ascii(cls):
    raise NotImplementedError

