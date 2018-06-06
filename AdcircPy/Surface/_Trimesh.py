import numpy as np
from matplotlib.path import Path
from scipy.interpolate import griddata
import pyproj

def get_extent(self, **kwargs):
    epsg = kwargs.pop("epsg", self.epsg)

    if epsg != self.epsg:
        self_proj = pyproj.Proj(init="epsg:{}".format(self.epsg))
        target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
        x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    else:
        x = self.x
        y = self.y

    return [np.min(x), np.max(x), np.min(y), np.max(y)]

def get_values_at_lonlat(self, lon, lat, list_idx=0, method='linear'):
    if isinstance(self.values, list):
        values = self.values[list_idx]
    else:
        values = self.values
    if np.ma.is_masked(values):
        values = np.ma.filled(values, 0.0)
    elif np.isin(values, -99999.0).any():
        idx = np.where(np.isin(values, -99999.0))
        values[idx] = 0.0    
    if method != 'force':
        return griddata((self.x,self.y),values,(lon,lat), method=method)
    else:
        idx = np.where(~np.isin(values, 0.0))
        return griddata((self.x[idx], self.y[idx]), values[idx], (lon, lat), method='nearest')


def get_extent_idx(self, extent, epsg, **kwargs):
    # epsg   = kwargs.pop("epsg", self.epsg)
    # extent = kwargs.pop("extent", self.get_extent(epsg=epsg))
    if epsg != self.epsg:
        self_proj = pyproj.Proj(init="epsg:{}".format(self.epsg))
        target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
        x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    else:
        x = self.x
        y = self.y
    if isinstance(extent, list) or isinstance(extent, tuple):
        bound_box = np.logical_and(
                        np.logical_and(x>=extent[0], x<=extent[1]),
                        np.logical_and(y>=extent[2], y<=extent[3]))
        idx, = np.where(bound_box)
    elif isinstance(extent, Path):
        idx, = np.where(extent.contains_points(np.vstack((x,y)).T))
    return idx

def plot_trimesh(self, extent=None, axes=None, title=None, color='black', linewidth=0.5, alpha=0.4):
    axes, idx = fig._init_fig(self, axes, extent, title)
    axes.triplot(self.x, self.y, self.elements, color=color, linewidth=linewidth, alpha=alpha)
    return axes

def get_xyz(self, **kwargs):
    idx = self.get_extent_idx(self.get_extent(), self.epsg)
    return np.vstack((self.x[idx], self.y[idx], self.values[idx])).T


def get_xy(self, **kwargs):
    idx = self.get_extent_idx(self.get_extent(), self.epsg)
    return np.vstack((self.x[idx], self.y[idx])).T


def get_elements_surrounding_node(self, node_index):
    adjacent_element_idxs, = np.where(np.any(np.isin(self.elements, node_index), axis=1))
    return self.elements[adjacent_element_idxs]

def get_element_paths(self, extent=None):
    if extent is None: extent = self.get_extent()
    mask_x  = np.logical_and(self.x > extent[0], self.x < extent[1])
    mask_y  = np.logical_and(self.y > extent[2], self.y < extent[3])
    masked  = np.ma.masked_where(np.logical_and(mask_x, mask_y), self.values)
    trimask = np.all(masked.mask[self.elements], axis=1)
    idx,  = np.where(trimask) 
    paths = list()
    for i in idx:
        vertex = list()
        codes = [Path.MOVETO]
        for j in [0,1,2]:
            vertex.append((self.x[self.elements[i][j]], self.y[self.elements[i][j]]))
            codes.append(Path.LINETO)
        vertex.append(vertex[0])
        codes[-1] = Path.CLOSEPOLY
        paths.append(Path(vertex, codes))
    return paths