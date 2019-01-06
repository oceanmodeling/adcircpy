from AdcircPy.DEM import _BaseDEM

class _DEMDiff(_BaseDEM):

  def __init__(self, x, y, values, geoTransform, epsg):
    super(_DEMDiff, self).__init__(x, y, values, geoTransform, epsg)

  def make_plot(self, axes=None, vmin=None, vmax=None, title=None):
    self._init_make_plot('diff', axes, vmin, vmax, figsize)
    self._ax = self._axes.pcolormesh(self.x, self.y, np.flipud(self.values), cmap=self._cmap, norm=self._norm)
    self._finalize_make_plot(title, show, r'elevation [$\Delta$m]', 'neither')