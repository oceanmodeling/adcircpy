from netCDF4 import Dataset
from AdcircPy import Outputs
from AdcircPy.Model import AdcircMesh

def from_netcdf(cls, path, fort14=None, fort15=None, datum='MSL', epsg=None, datum_grid=None):
  nc = Dataset(path)
  if 'zeta_max' not in nc.variables.keys():
    raise Exception('Not a maxele file!')
  if fort14 is not None:
    if isinstance(fort14, AdcircMesh)==False:
      fort14 = AdcircMesh.from_fort14(fort14, fort15)
  return cls(nc['x'][:],
             nc['y'][:],
             nc['element'][:]-1,
             nc['zeta_max'][:],
             nc['time_of_zeta_max'],
             datum=datum,
             # **fort14.get_boundary_dictionary() # Need to add the additional infor provided by fort.14
             datum_grid=datum_grid)

def from_ascii(cls):
  raise NotImplementedError

