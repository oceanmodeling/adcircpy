import numpy as np
from AdcircPy.Mesh import AdcircMesh
from AdcircPy import Outputs

def _from_ascii(path, fort14, fort15=None, datum='MSL', epsg=4326):
  f = open(path)
  line = f.readline().strip()
  line = f.readline().split()
  number_of_datasets = int(line[0])
  number_of_points = int(line[1])
  line = f.readline()
  values = list()
  for i in range(number_of_points):
    line = f.readline().split()
    if len(line)==2:
      values.append(float(line[1]))
    elif len(line)==3:
      raise NotImplementedError('Is this a velocity output? Great moment to start coding!')
      # values.append((float(line[1]), float(line[2])))
  values = np.ma.masked_equal(values, -99999.)
  if isinstance(fort14, AdcircMesh):
    pass
  else:
    fort14 = AdcircMesh.from_fort14(fort14, fort15=fort15, epsg=epsg, datum=datum)
  if number_of_datasets==1:
    f.close()
    f=None
  return Outputs.__OutputSurface(fort14.x,
                                 fort14.y,
                                 values,
                                 fort14.elements,
                                 epsg=fort14.epsg,
                                 datum=fort14.datum,
                                 f=f,
                                 ocean_boundaries = fort14.ocean_boundaries,
                                 land_boundaries =  fort14.land_boundaries,
                                 inner_boundaries = fort14.inner_boundaries,
                                 weir_boundaries =  fort14.weir_boundaries,
                                 inflow_boundaries = fort14.inflow_boundaries,
                                 outflow_boundaries = fort14.outflow_boundaries,
                                 culvert_boundaries = fort14.culvert_boundaries)

def make_animation(self, **kwargs):
  axes, idx = _fig.init_fig(self, axes, extent, title, epsg)
  start_slice = kwargs.pop("start_slice", 0)
  stop_slice  = kwargs.pop("stop_slice", len(self.timestep))
  vals = list()
  for values in self.values[start_slice:stop_slice]:
      vals.append(values[idx])
  vals = np.asarray(vals)
  vals = np.ma.masked_where(vals==-99999.0, vals)
  vmin = kwargs.pop("vmin", np.min(vals))
  vmax = kwargs.pop("vmax", np.max(vals))
  return axes

