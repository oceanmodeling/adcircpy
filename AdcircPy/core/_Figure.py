import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from AdcircPy.core._FixPointNormalize import _FixPointNormalize

class _Figure(object):
  def __init__(self):
    super(_Figure, self).__init__()

  def _init_fig(self, axes=None, figsize=None):
    self._plt = plt
    if axes is None:
      self._axes = self._plt.figure(figsize=figsize).add_subplot(111)
    else:
      self._axes=axes

  def _init_vmin(self, vmin):
    if vmin is None:
      vmin = np.ceil(np.min(self.values))
    self._vmin = vmin

  def _init_vmax(self, vmax):
    if vmax is None:
      vmax = np.ceil(np.max(self.values))
    self._vmax = vmax

  def _init_cmap_levels(self, cmap, colors):
    if cmap is None:
      self._cmap = self._plt.cm.get_cmap('viridis')
      self._levels=np.linspace(self._vmin, self._vmax, colors)
      self._col_val=0.

    elif cmap == 'topobathy':
      if self._vmax < 0.:
        self._cmap = self._plt.cm.seismic
        self._col_val = 0.
      else:
        wet_count = int(np.floor(256.*(float((self.values < 0.).sum())/float(self.values.size))))
        self._col_val = float(wet_count)/256.
        dry_count = 256. - wet_count
        colors_undersea = self._plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = self._plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        self._cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(self._vmin, 0.0, wet_count, endpoint=False)
        dlevels = np.linspace(0.0, self._vmax, dry_count)
        self._levels = np.hstack((wlevels, dlevels))
    else:
      self._cmap = self._plt.cm.get_cmap(cmap)
      self._col_val = 0.
      self._levels=np.linspace(self._vmin, self._vmax, colors)

  def _init_norm(self):
    self._norm = _FixPointNormalize(sealevel=0.0, vmax=self._vmax, vmin=self._vmin, col_val=self._col_val)

  def _init_title(self, title):
    if title is not None:
      self._axes.set_title(title)

  def _auto_show(self, show):
    if show==True:
      self._plt.show()

  def _auto_scaling(self, extent=None):
    self._axes.axis('scaled')
    if extent is not None:
      self._axes.axis(extent)

  def _init_cbar(self, cmap_extend):
    mappable = ScalarMappable(cmap=self._cmap)
    mappable.set_array([])
    mappable.set_clim(self._vmin, self._vmax)
    divider = make_axes_locatable(self._axes)
    cax = divider.append_axes("bottom", size="2%", pad=0.5)
    self._cbar = self._plt.colorbar(mappable, cax=cax, extend=cmap_extend, orientation='horizontal')
    self._cbar.set_ticks([self._vmin, self._vmin + self._col_val *(self._vmax-self._vmin), self._vmax])
    self._cbar.set_ticklabels([np.around(self._vmin, 2), 0.0, np.around(self._vmax, 2)])

  def _init_cbar_label(self, cbar_label=None):
    if cbar_label is not None:
      self._cbar.set_label(cbar_label)

##########

# def _auto_color(self):

#   elif plot_type == 'diff':
#     # mlevel = np.mean(self.values)
#     self._cmap = self._plt.cm.seismic
#     self._col_val = 0.5
#     self._levels=np.linspace(self._vmin, self._vmax, 256)
  
#   self._norm = _FixPointNormalize(sealevel=0.0, vmax=self._vmax, vmin=self._vmin, col_val=self._col_val)

# def _init_cmap_extend(self, vmin, vmax):
#   if vmin == 'auto' and vmax=='auto':
#     self._cmap_extend = 'neither'
#   if vmin > self._vmin:
#     pass

