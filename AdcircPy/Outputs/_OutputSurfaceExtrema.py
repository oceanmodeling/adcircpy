import numpy as np
import matplotlib.pyplot as plt
from AdcircPy.Mesh import AdcircMesh, _UnstructuredGrid
from AdcircPy import Outputs


def _make_plot(self, **kwargs):
  _UnstructuredGrid._make_plot(self, show=True)
  plt.tripcolor(self._Tri, self._times)
  plt.show()



# def _from_ascii(path, fort14, fort15=None, datum='MSL', epsg=4326):
#   f = open(path)
#   line = f.readline().strip()
#   line = f.readline().split()
#   number_of_datasets = int(line[0])
#   number_of_points = int(line[1])
#   line = f.readline()
#   values = list()
#   for i in range(number_of_points):
#     line = f.readline().split()
#     if len(line)==2:
#       values.append(float(line[1]))
#     elif len(line)==3:
#       raise NotImplementedError('Is this a velocity output? Great moment to start coding!')
#       # values.append((float(line[1]), float(line[2])))
#   values = np.ma.masked_equal(values, -99999.)
#   if isinstance(fort14, AdcircMesh):
#     pass
#   else:
#     fort14 = AdcircMesh.from_fort14(fort14, fort15=fort15, epsg=epsg, datum=datum)
#   if number_of_datasets==1:
#     f.close()
#     f=None
#   return Outputs._OutputSurface(fort14.x,
#                                  fort14.y,
#                                  values,
#                                  fort14.elements,
#                                  epsg=fort14.epsg,
#                                  datum=fort14.datum,
#                                  f=f,
#                                  ocean_boundaries = fort14.ocean_boundaries,
#                                  land_boundaries =  fort14.land_boundaries,
#                                  inner_boundaries = fort14.inner_boundaries,
#                                  weir_boundaries =  fort14.weir_boundaries,
#                                  inflow_boundaries = fort14.inflow_boundaries,
#                                  outflow_boundaries = fort14.outflow_boundaries,
#                                  culvert_boundaries = fort14.culvert_boundaries)

