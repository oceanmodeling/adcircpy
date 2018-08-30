# from collections import OrderedDict
from AdcircPy import Mesh
from AdcircPy import Outputs
from AdcircPy.core.Tides import TidalDB as _TidalDB

def _from_fort51(path, fort14, fort15, **kwargs):
  TidalDB = _TidalDB()
  if fort14 is not None:
    if isinstance(fort14, Mesh.AdcircMesh):
      if isinstance(fort14.fort15, Mesh.fort15):
        pass
      else:
        raise Exception('If passing an AdcircMesh instance, fort15 must be initialized!')
    else:
      if fort15 is None:
        raise Exception('fort.51 files lack metadata and therefore both fort14 and fort15 are required for analysis.')
      fort14 = Mesh.AdcircMesh.from_fort14(fort14, fort15=fort15, **kwargs)
  _f = open(path)
  def __line():
    return _f.readline()
  _num_of_comp = int(__line())
  _components = list()
  for i in range(_num_of_comp):
    _components.append(__line().split()[-1].strip())
  _num_of_stations = int(__line())
  stations = dict()
  for _station in range(_num_of_stations):
    _id = int(__line()) - 1
    if fort15 is not None:
      _id = fort14.fort15['StationOutputs']['elevation']['id'][_id]
    else:
      _id = str(_id)
    stations[_id] = dict()
    for _component in _components:
      line = __line().split()
      amplitude = float(line[0])
      phase = float(line[1])
      stations[_id][_component] = {
              "orbital_frequency" : TidalDB[_component]["orbital_frequency"],
              "amplitude"         : amplitude,
              "phase"             : phase }
  _f.close()
  return Outputs.HarmonicConstituentsStations(**stations)
    




  
