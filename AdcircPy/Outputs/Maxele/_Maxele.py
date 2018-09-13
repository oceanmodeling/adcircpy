from netCDF4 import Dataset
from AdcircPy import Outputs
from AdcircPy.Mesh import AdcircMesh

def _from_netcdf(path, fort14=None):
  _Dataset = Dataset(path)
  if 'zeta_max' not in _Dataset.variables.keys():
    raise Exception('Not a maxele file!')
  # if fort14 is not None and isinstance(fort14, AdcircMesh)==False:
    # fort14 = AdcircMesh.from_fort14(fort14)
  kwargs = dict()
  kwargs['Dataset'] = _Dataset
  return Outputs.Maxele(_Dataset['x'][:],
                        _Dataset['y'][:],
                        _Dataset['zeta_max'][:],
                        _Dataset['element'][:]-1,
                        **kwargs)




