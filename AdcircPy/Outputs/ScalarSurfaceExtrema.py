import numpy as np
import matplotlib.pyplot as plt
import abc
from netCDF4 import Dataset
from AdcircPy.Mesh import TrimeshSurface
from AdcircPy import Outputs

class ScalarSurfaceExtrema(TrimeshSurface):
  """
  Subclass representing a general Output Surface Extrema, which is instantiated by the _OutputFactory class.
  """
  def __init__(self, x, y, elements, values, times, epsg=4326, vertical_datum='LMSL', nodeID=None, elementID=None, **Boundaries):
    self._times = times
    TrimeshSurface.__init__(self, x, y, elements, values, epsg, nodeID=nodeID, elementID=elementID, vertical_datum=vertical_datum, **Boundaries)

  def make_plot(self, title='Surface Extrema', axes=None, plot_times=False, **kwargs):
    if axes is None:
      fig1 = plt.figure()
      ax1 = fig1.add_subplot(111)
    else:
      ax1 = axes
    TrimeshSurface.make_plot(self, axes=ax1, title=title)
    if len(self._times)>0 and plot_times==True:
      fig2 = plt.figure()
      ax2 = fig2.add_subplot(111)
      _ax = ax2.tripcolor(self._Tri, self._times)
      ax2.set_title('Time in seconds after coldstart at which extrema happened.')
      ax2.axis('scaled')
      plt.colorbar(_ax)

  @staticmethod
  def is_netcdf(path):
    try:
      Dataset(path)
      return True
    except:
      return False
