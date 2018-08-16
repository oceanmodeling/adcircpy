import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from haversine import haversine
# from AdcircPy.core import Figures
from AdcircPy import Mesh
from AdcircPy.Mesh import _fort14
from AdcircPy.Mesh import _fort13

def _from_fort14(fort14, datum='MSL', epsg=4326, fort13=None, fort15=None):
    kwargs = _fort14.parse_fort14(fort14)
    kwargs['datum'] = datum
    kwargs['epsg']  = epsg
    kwargs['fort13'] = fort13
    kwargs['fort15'] = fort15
    return Mesh.AdcircMesh(**kwargs)

def make_plot(self, surface='bathy', fort13=None, **kwargs):
    if fort13 is not None:
        surface_types = [ x for x in self.fort13.keys() ]
        surface_types.append('bathy')
    else:
        surface_types = ['bathy']
    
    if surface not in surface_types:
        raise ValueError("Specified surface not found! Available surfaces are {}".format(surface_types))

    if surface=='bathy':
        return Mesh._Mesh.plot_bathy(self, **kwargs)
    else:
        return #call the general surface module.



def plot_bathy(self, extent=None, epsg=None, axes=None, title=None, total_colors=256, cbar_label=r'elevation [m]', **kwargs):
    axes, idx = Figures._init_fig(self, axes, extent, title, epsg)
    vmin = kwargs.pop("vmin", np.min(self.values[idx]))
    vmax = kwargs.pop("vmax", np.max(self.values[idx]))
    show = kwargs.pop("show", False)
    if vmax < 0.:

        levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
        mlevel = np.mean(self.values[idx])
        cmap   = kwargs.pop("cmap", plt.cm.seismic)
        col_val = 0.5      
    else:
        mlevel = 0.
        wet_count = int(np.floor(total_colors * \
                (float((self.values[idx] < 0.).sum())/self.values[idx].size)))
        col_val = float(wet_count)/float(total_colors)
        dry_count = total_colors - wet_count
        colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(vmin, mlevel, wet_count, endpoint=False)
        dlevels = np.linspace(mlevel, vmax, dry_count)
        levels = np.hstack((wlevels, dlevels))
    norm = _fig.FixPointNormalize(sealevel=mlevel, vmax=vmax, vmin=vmin, col_val=col_val)
    axes.tricontourf(self.x, self.y, self.elements, self.values,
                    levels=levels, cmap=cmap, norm=norm, extend='both', **kwargs)
    cbar = _fig.init_colorbar(axes, cmap, vmin, vmax)
    cbar.set_label(cbar_label)
    cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), mlevel, np.around(vmax, 2)])
    if show == True:
        plt.show()
    return axes

def plot_shoreline(self, extent=None, epsg=None, axes=None, title=None, color='black', linewidth=0.5):
    axes, idx = Figures._init_fig(self, axes, extent, title, epsg)
    if type(self.values) is list:
        if self.bathymetry is None:
            raise TypeError("fort.14 required to plot shoreline")
        values = self.bathymetry
    else:
        values = self.values

    axes.tricontour(self.x, self.y, self.elements, values, levels=[0.], colors=color, linewidths=linewidth)
    return axes

def get_mean_value(self, extent=None, step=0):
    if extent is None:
        extent = self.get_extent()
    idx = self.get_extent_idx(extent)
    if type(self.values) == type(list()):
        values = self.values[step][idx]
    else:
        values = self.values[idx]
    return np.mean(values)

def get_values_under_Path(self, path, **kwargs):
    xin = self.x
    yin = self.y
    if len(self.values.shape) == 3:
        zin = self.values[:,0,timestep]
    else:
        zin = self.values
    idx, = np.where(np.logical_and(
        np.logical_and(xin>=np.min(path.vertices[:,0]), xin<=np.max(path.vertices[:,0])),
        np.logical_and(yin>=np.min(path.vertices[:,1]), yin<=np.max(path.vertices[:,1]))))
    return griddata((xin[idx], yin[idx]), zin[idx], (path.vertices[:,0],path.vertices[:,1]))


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