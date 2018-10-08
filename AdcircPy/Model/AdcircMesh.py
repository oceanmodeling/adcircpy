import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from haversine import haversine
from AdcircPy.core._FixPointNormalize import FixPointNormalize
from AdcircPy.Model import _fort14
from AdcircPy import Model


def from_fort14(cls, fort14, datum='MSL', epsg=4326, fort13=None, fort15=None, datum_grid=None):
  kwargs = _fort14.parse_fort14(fort14)
  kwargs['datum'] = datum
  kwargs['epsg']  = epsg
  kwargs['fort13'] = fort13
  kwargs['fort15'] = fort15
  kwargs['datum_grid'] = datum_grid
  return cls(**kwargs)

def make_plot(self, extent=None, epsg=None, axes=None, title=None, total_colors=256, cbar_label='elevation [m]', **kwargs):
  self._init_fig(axes, extent, title, epsg)
  vmin = kwargs.pop("vmin", np.min(self.values[self._idx]))
  vmax = kwargs.pop("vmax", np.max(self.values[self._idx]))
  show = kwargs.pop("show", False)
  if vmax < 0.:
    levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
    mlevel = np.mean(self.values[self._idx])
    cmap   = kwargs.pop("cmap", plt.cm.seismic)
    col_val = 0.5      
  else:
    mlevel = 0.
    wet_count = int(np.floor(total_colors * \
            (float((self.values[self._idx] < 0.).sum())/self.values[self._idx].size)))
    col_val = float(wet_count)/float(total_colors)
    dry_count = total_colors - wet_count
    colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
    colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
    colors = np.vstack((colors_undersea, colors_land))
    cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
    wlevels = np.linspace(vmin, mlevel, wet_count, endpoint=False)
    dlevels = np.linspace(mlevel, vmax, dry_count)
    levels = np.hstack((wlevels, dlevels))
  norm = FixPointNormalize(sealevel=mlevel, vmax=vmax, vmin=vmin, col_val=col_val)
  self._axes.tricontourf(self.x, self.y, self.elements, self.values,
                         levels=levels, cmap=cmap, norm=norm, extend='both', **kwargs)
  self._init_cbar(cmap, vmin, vmax)
  self._cbar.set_label(cbar_label)
  self._cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
  self._cbar.set_ticklabels([np.around(vmin, 2), mlevel, np.around(vmax, 2)])
  if show == True:
    plt.show()


def interpolate_DEM(self, tile, **kwargs):
    method = kwargs.pop("method", "FVM")
    channel_polygons = kwargs.pop("channel_polygons", None)
    bar_polygons = kwargs.pop("bar_polygons", None)
    tile_extent = tile.get_extent()
    idxs_of_points_in_tile = self.get_extent_idx(tile_extent, tile.epsg)
    xyz = tile.get_xyz(epsg=self.epsg)
    for idx in idxs_of_points_in_tile:
        path = self.get_finite_volume_Path(idx)
        _idx, = np.where(np.logical_and(
                    np.logical_and(xyz[:,0]>=np.min(path.vertices[:,0]), xyz[:,0]<=np.max(path.vertices[:,0])),
                    np.logical_and(xyz[:,1]>=np.min(path.vertices[:,1]), xyz[:,1]<=np.max(path.vertices[:,1]))))
        values = xyz[_idx,:]
        if method=="FVM":
            _values = values[np.where(path.contains_points(values[:,0:2])),2]
            self.values[idx] = np.mean(_values)
        else:
            self.values[idx] = griddata((xyz[_idx,0], xyz[_idx,1]), xyz[_idx,2], (self.x[idx], self.y[idx]), method=method)
        # path.contains_points() takes time, so we try to recycle the _values variable with try/except
        if channel_polygons is not None:
            for channel_polygon in channel_polygons:
                if channel_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(channel_polygon):
                    try: 
                        _values
                    except:
                        _values = values[np.where(path.contains_points(values[:,0:2])),2]
                    self.values[idx] = np.min(_values)
    
        if bar_polygons is not None:
            for bar_polygon in bar_polygons:
                if bar_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(bar):
                    try:
                        _values
                    except:
                        _values = values[np.where(path.contains_points(values[:,0:2])),2]
                    self.values[idx] = np.max(_values)
