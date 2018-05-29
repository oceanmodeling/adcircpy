import numpy as np
from matplotlib.path import Path
from scipy.interpolate import griddata

def get_extent(self):
    return [np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]

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


def get_extent_idx(self, extent=None):
    if extent is None:
        extent = self.get_extent()
    bound_box = np.logical_and(
                    np.logical_and(self.x>=extent[0], self.x<=extent[1]),
                    np.logical_and(self.y>=extent[2], self.y<=extent[3]))
    # if remove_values is not None:
    #     if np.ma.is_masked(self.values):
    #         _values = self.values.data
    #     else:
    #         _values = self.values
    #     idx, = np.where(np.logical_and(bound_box, _values!=remove_values))
    # else:
    idx, = np.where(bound_box)
    return idx



def get_xyz(self, extent=None, radius=None):

    if isinstance(extent, Path):
        idx, = np.where(extent.contains_points(self.get_xy(), radius=radius))
    else:
        idx = self.get_extent_idx(extent)
    return np.vstack((self.x[idx], self.y[idx], self.values[idx])).T


def get_xy(self, min_x=None, min_y=None, max_x=None, max_y=None, epsg=None):
    return np.vstack((self.x[idx], self.y[idx], self.values[idx])).T

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