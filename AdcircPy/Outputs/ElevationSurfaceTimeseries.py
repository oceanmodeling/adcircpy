from netCDF4 import Dataset
import numpy as np
from matplotlib.animation import FuncAnimation
from AdcircPy.Model import UnstructuredMesh

class ElevationSurfaceTimeseries(UnstructuredMesh):
  def __init__(self, Dataset, **kwargs):
    super(ElevationSurfaceTimeseries, self).__init__(Dataset['x'][:].flatten(), Dataset['y'][:].flatten(), Dataset['element'][:]-1, Dataset['zeta'][0,:].flatten(), **kwargs)
    self._Dataset = Dataset
    self._index = 0

  @property
  def values(self):
    return self._values

  def __next_value(self):
    if self._index < self._Dataset['zeta'].shape[0]:
      self._index += 1
    else:
      self._index = 0
    values = self._Dataset['zeta'][self._index,:].flatten()
    self._values = np.ma.masked_where(values==-99999.0, values)
    self._update_Tri()


  @classmethod
  def from_netcdf(cls, path, fort14=None, datum=None, epsg=None, datum_grid=None):
    if fort14 is not None:
      raise NotImplementedError('when adding a fort.14 more coding is needed.')
    return cls(Dataset(path))

  def plot_surface_animation(self, axes=None, extent=None, title='Surface level',
                            vmin=None, vmax=None, start_slice=0, end_slice=None,
                            speedup=25.0, show=False):
    self._init_fig(axes, extent, title, self.epsg)
    if vmin is None:
      values = self._Dataset['zeta'][-1,:]
      vmin = np.min(np.ma.masked_where(values==-99999.0, values))
    if vmax is None:
      values = self._Dataset['zeta'][-1,:]
      vmax = np.max(np.ma.masked_where(values==-99999.0, values))
    if end_slice == None:
      end_slice   = self._Dataset['zeta'].shape[0]
    # Need to figure out time interval between slices.
    time_interval = None
    print(vmin, vmax)
    self._values = self._Dataset['zeta'][-1,:]
    self._values = np.ma.masked_where(self._values==-99999.0, self._values)
    self._update_Tri()
    _ax = self._axes.tripcolor(self._Tri, self._values, vmin=vmin, vmax=vmax)
    # self._animation = FuncAnimation(self.plt.gcf(), self.__update_animation(_ax),
    #                                 frames = np.arange(start_slice+1, end_slice),
    #                                 # interval = time_interval
    #                                 )
    # mappable = ScalarMappable(cmap=self.wavecmap)
    # divider = make_axes_locatable(self._axes)
    # cax = divider.append_axes("right", size="2%", pad=0.5)
    # mappable.set_array([])
    # mappable.set_clim(_vmin, _vmax)
    # plt.colorbar(mappable, cax=cax, extend='both', orientation='vertical')
    
    if extent is None:
      extent = self.get_extent()
    self._axes.axis('scaled')
    self._axes.axis(extent)
    self.plt.colorbar(_ax)
    if show ==True:
      self.plt.show()

  
  def __update_animation(self, _ax):
    self.__next_value()
    return _ax.set_array(self._values)