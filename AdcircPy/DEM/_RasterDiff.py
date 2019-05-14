# global imports
import numpy as np

# local imports
from AdcircPy.DEM._BaseRaster import _BaseRaster


class _RasterDiff(_BaseRaster):

    def __init__(self, x, y, values, GeoTransform, SpatialReference):
        super(_RasterDiff, self).__init__(x, y, values, GeoTransform,
                                          SpatialReference)

    def make_plot(self, axes=None, vmin=None, vmax=None, title=None,
                  figsize=(10, 10), show=True):
        self._init_make_plot('diff', axes, vmin, vmax, figsize)
        self._ax = self._axes.pcolormesh(self.x, self.y,
                                         np.flipud(self.values),
                                         cmap=self._cmap, norm=self._norm)
        self._finalize_make_plot(title, show, r'elevation [$\Delta$m]',
                                 'neither')
