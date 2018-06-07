import numpy as np
from matplotlib.path import Path
from scipy.interpolate import griddata
import pyproj
from haversine import haversine

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

def get_values_at_lonlat(self, lon, lat, step=0, method='linear'):
    if isinstance(self.values, list):
        values = self.values[step]
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

def get_finite_volume_Path(self, idx, radius=None):
    midpoints = list()
    centroids = list()
    adjacent_elements = self.get_elements_surrounding_node(idx)
    # Iterate over surrounding elements
    for j, element in enumerate(adjacent_elements):
        _element = list(element)
        _element.append(_element[0])
        _midpoints = list()
        # Calculate element midpoints and centroid
        for i, point_idx in enumerate(element):         
            dx =  self.x[_element[i+1]] - self.x[_element[i]]
            dy =  self.y[_element[i+1]] - self.y[_element[i]]
            if np.isin(idx, [_element[i], _element[i+1]]):
                midpoints.append((self.x[_element[i]] + 0.5*dx, self.y[_element[i]] + 0.5*dy))
            _midpoints.append((self.x[_element[i]] + 0.5*dx, self.y[_element[i]] + 0.5*dy))
        _midpoints = np.array(_midpoints)
        centroids.append((np.mean(_midpoints[:,0]),np.mean(_midpoints[:,1])))
        _element = list(element)
        while _element[0] != idx:
            _element = list(np.roll(_element,1))
        adjacent_elements[j] = _element
    adjacent_elements = np.array(adjacent_elements)
    lateral_indexes, counts = np.unique(adjacent_elements, return_counts=True)
    vertices = list()
    for centroid in centroids:
        vertices.append(centroid)
    for midpoint in midpoints:
        vertices.append(midpoint)    
    ordered_vertices = list()
    if 1 in counts: # means that this node is at a boundary and need to be included in the polygon.
        ordered_vertices.append((self.x[idx], self.y[idx]))
        index = list(np.delete(lateral_indexes, np.where(counts>1))).pop()
        x = self.x[idx]
        y = self.y[idx]
        dx = self.x[index] - self.x[idx]
        dy = self.y[index] - self.y[idx]
        ordered_vertices.append(vertices.pop(vertices.index((x+0.5*dx,y+0.5*dy))))
    else: # means that this is node is fully surrounded by elements so the obs point is not included in polygon.
        ordered_vertices.append(vertices.pop())
    lon = ordered_vertices[-1][0]
    lat = ordered_vertices[-1][1]
    while vertices:
        diff = list()
        for vertex in vertices:
            _diff = haversine((vertex[1], vertex[0]), (lat, lon))
            diff.append(_diff)
        ordered_vertices.append(vertices.pop(diff.index(min(diff))))
        lon = ordered_vertices[-1][0]
        lat = ordered_vertices[-1][1]
    ordered_vertices.append(ordered_vertices[0])
    return Path(ordered_vertices, closed=True)