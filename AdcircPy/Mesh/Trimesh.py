from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.cm import ScalarMappable
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import pyproj
from haversine import haversine
from AdcircPy.core import _Figure

class Trimesh(_Figure):   
  def __init__(self, x, y, elements, epsg, nodeID=None, elementID=None,
              ocean_boundaries=None, land_boundaries=None, inner_boundaries=None,
              inflow_boundaries=None, outflow_boundaries=None,
              weir_boundaries=None, culvert_boundaries=None):
    super(Trimesh, self).__init__()
    self._init_inputs(x, y, elements)
    self._epsg           = epsg
    self._init_nodeID(nodeID)
    self._init_elementID(elementID)
    self._init_KDTree()
    self.ocean_boundaries    = ocean_boundaries
    self.land_boundaries     = land_boundaries
    self.inner_boundaries    = inner_boundaries
    self.inflow_boundaries   = inflow_boundaries
    self.outflow_boundaries  = outflow_boundaries
    self.weir_boundaries     = weir_boundaries
    self.culvert_boundaries = culvert_boundaries

  @property
  def x(self):
    return self._x
  
  @property
  def y(self):
    return self._y

  @property
  def elements(self):
    return self._elements
  
  @property
  def epsg(self):
    return self._epsg
  
  def get_xy(self, extent=None, epsg=None):
    if extent is None:
      extent = self.get_extent()
    if epsg is None:
      epsg = self.epsg
    idx = self.__get_extent_idx(extent, epsg)
    return np.vstack((self.x[idx], self.y[idx])).T

  def get_finite_volume_element_indexes(self, index):
    """ """
    return self.elements[np.where(np.any(np.isin(self.elements, index), axis=1))[0]]

  def get_finite_volume_elements_as_Path_list(self, index):
    """ """
    return [self.get_Path_from_element(element) for element in self.get_finite_volume_element_list(index)]

  def get_element_containing_coord(self, coord):
    x = coord[0]; y = coord[1]
    distance, node_idx = self.KDTree.query([x,y])
    elements = self.get_finite_volume_Path_list(node_idx)
    for i, element in enumerate(elements):
      if element.contains_point((x,y)):
        return self.get_finite_volume_element_list(node_idx)[i]

  def get_extent(self, epsg=None, **kwargs):
    if epsg is None:
      epsg = self.epsg
    if epsg != self.epsg:
      self_proj = pyproj.Proj(init="epsg:{}".format(self.epsg))
      target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
      x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    else:
      x = self.x
      y = self.y
    return [np.min(x), np.max(x), np.min(y), np.max(y)]
  
  def get_Path_from_element(self, element):
      return Path([[self.x[element[0]], self.y[element[0]]],
                   [self.x[element[1]], self.y[element[1]]],
                   [self.x[element[2]], self.y[element[2]]],
                   [self.x[element[0]], self.y[element[0]]]], closed=True)

  def get_elements_in_extent(self, extent):
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

  def transform_to_epsg(self, epsg):
    self_proj = pyproj.Proj(init="epsg:{}".format(self.epsg))
    target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
    x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    self._x = np.asarray(x).flatten()
    self._y = np.asarray(y).flatten()
    self.epsg = epsg

  def make_plot(self, extent=None, epsg=None, axes=None, title=None, show=False, color='black', linewidth=0.5, alpha=0.4):
    self._init_fig(axes, figsize)
    self._axes.triplot(self.x, self.y, self.elements, color=color, linewidth=linewidth, alpha=alpha)
    self._auto_scaling()
    self._init_title(title)
    self._auto_show(show)

  def _init_inputs(self, x, y, elements):
    self._x              = np.asarray(x)
    self._y              = np.asarray(y)
    self._elements       = np.asarray(elements)
    if self._x.size != self._y.size != self._elements.shape[0]:
      raise Exception('Inconsistent sizes for x, y and element arrays.')

  def _init_nodeID(self, nodeID):
    if nodeID is not None:
      nodeID = np.asarray(nodeID)
      if self._x.size != self._y.size != nodeID.size:
        raise Exception('Inconsistent nodeID array.')
    else:
      nodeID = np.arange(self._x.size)
    self._nodeID = nodeID

  def _init_elementID(self, elementID):
    if elementID is not None:
      elementID = np.asarray(elementID)
      if self._elements.shape[0] != elementID.size:
        raise Exception('Inconsistent elementID array.')
    else:
      elementID = np.arange(self._elements.shape[0])
    self._elementID = elementID

  def _init_KDTree(self):
    """  """
    self.KDTree = cKDTree(self.get_xy())

  def __get_extent_idx(self, extent=None, epsg=None):
    # epsg   = kwargs.pop("epsg", self.epsg)
    # extent = kwargs.pop("extent", self.get_extent(epsg=epsg))
    if epsg is None:
      epsg=self.epsg
    if extent is None:
      extent = self.get_extent(epsg=epsg)
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

  def build_outer_polygon(self):
    boundary_list = list()
    if self.ocean_boundaries is not None:
      for i in range(len(self.ocean_boundaries)):
        vertices = list()
        codes = [Path.MOVETO]
        for j in range(self.ocean_boundaries[i].shape[0]):
          idx = self.ocean_boundaries[i][j]
          vertices.append((self.x[idx], self.y[idx]))
          codes.append(Path.LINETO)
        path = Path(vertices, codes[:-1])
        boundary_list.append(path)

    if self.land_boundaries is not None:
      for i in range(len(self.land_boundaries)):
        vertices = list()
        codes = [Path.MOVETO]
        for j in range(self.land_boundaries[i][0].shape[0]):
          idx = self.land_boundaries[i][0][j]
          vertices.append((self.x[idx], self.y[idx]))
          codes.append(Path.LINETO)
        path = Path(vertices, codes[:-1])
        boundary_list.append(path)

    if self.inflow_boundaries is not None:
      for i in range(len(self.inflow_boundaries)):
        vertices = list()
        codes = [Path.MOVETO]
        for j in range(self.inflow_boundaries[i][0].shape[0]):
          idx = self.inflow_boundaries[i][0][j]
          vertices.append((self.x[idx], self.y[idx]))
          codes.append(Path.LINETO)
        path = Path(vertices, codes[:-1])
        boundary_list.append(path)
        
    if self.outflow_boundaries is not None:
      for i in range(len(self.outflow_boundaries)):
        vertices = list()
        codes = [Path.MOVETO]
        for j in range(self.outflow_boundaries[i][0].shape[0]):
          idx = self.outflow_boundaries[i][0][j]
          vertices.append((self.x[idx], self.y[idx]))
          codes.append(Path.LINETO)
        path = Path(vertices, codes[:-1])
        boundary_list.append(path)
    # if self.weirBoundaries is not None:
        # vertices=list()
        # for boundary in self.weirBoundaries:
            # vertices.append(np.vstack(((self.x[boundary['front_face']], self.y[boundary['front_face']]))).T)
            # vertices.append(np.vstack(((self.x[boundary['back_face']], self.y[boundary['back_face']]))).T)
        # print(vertices)
            
            # vertices=[self.]
            # for face in boundary['front_face']
            # boundary_list.append(path)
    try:
      ordered_list=[boundary_list.pop()]
    except:
      return     
    b1 = ordered_list[-1]
    b1_bottom_lon =  b1.vertices[-1,0]
    b1_bottom_lat =  b1.vertices[-1,1]
    while boundary_list:
      diff = list()    
      for boundary in boundary_list:
        eudiff = np.sqrt((b1_bottom_lon-boundary.vertices[0,0])**2. + (b1_bottom_lat-boundary.vertices[0,1])**2.)
        diff.append(eudiff)
      # append the boundary with the minimum euclidean difference, which we assume to be next in the sequence
      ordered_list.append(boundary_list.pop(diff.index(min(diff))))
      # get the new boundary limits to find the next shape in the sequence
      b1 = ordered_list[-1]
      b1_bottom_lon = b1.vertices[-1,0]
      b1_bottom_lat = b1.vertices[-1,1]
    final_list = list()
    for path in ordered_list:
      for i in range(len(path.vertices)):
        final_list.append((path.vertices[i,0],path.vertices[i,1]))
    codes = len(final_list) * [Path.LINETO]
    codes[0] = Path.MOVETO
    final_list.append(final_list[0])
    codes.append(Path.CLOSEPOLY)
    return Path(final_list, codes)

  def build_inner_polygons(self, epsg=None):
    if self.inner_boundaries is None:
      return
    if epsg is None:
      epsg = self.epsg

    if epsg != self.epsg:
      self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
      target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
    innerBoundaries = list()
    for i in range(len(self.inner_boundaries)):
      idxs =self.inner_boundaries[i][0]
      vertices = list()
      for idx in idxs:
        vertices.append((self.x[idx], self.y[idx]))
      if epsg != self.epsg:
        x = [x for x,y in vertices]
        y = [y for x,y in vertices]
        x, y = pyproj.transform(self_proj, target_proj, x, y)
        vertices = list(zip(x,y))
      innerBoundaries.append(Path(vertices, closed=True))
      # Uncomment to see a plot of the geometries as they are generated.
      # import matplotlib.pyplot as plt
      # import matplotlib.patches as patches
      # fig = plt.figure()
      # ax = fig.add_subplot(111)
      # patch = patches.PathPatch(path, facecolor='orange', lw=2)
      # ax.add_patch(patch)
      # ax.axis('scaled')
      # plt.show()
    return innerBoundaries

  def get_land_boundaries(self, epsg=None):
    if self.land_boundaries is not None:
      if epsg is None:
          epsg = self.epsg
      if epsg != self.epsg:
        self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
        target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
      boundary_list=list()
      for i in range(len(self.land_boundaries)):
        vertices = list()
        for j in range(self.land_boundaries[i][0].shape[0]):
          idx = self.land_boundaries[i][0][j]
          vertices.append((self.x[idx], self.y[idx]))
        if epsg != self.epsg:
          x = [x for x,y in vertices]
          y = [y for x,y in vertices]
          x, y = pyproj.transform(self_proj, target_proj, x, y)
          vertices = list(zip(x,y))
        boundary_list.append(Path(vertices, closed=True))
      return boundary_list


  def plot_outerBoundary(self, extent=None, axes=None, **kwargs):
    axes, idx = self._init_fig(self, axes, extent, title)
    patch = patches.PathPatch(self.outerBoundary, **kwargs)
    axes.add_patch(patch)

  def get_mean_value(self, extent=None, step=0):
    if extent is None:
      extent = self.get_extent()
    idx = self.get_extent_idx(extent)
    if type(self.values) == type(list()):
      values = self.values[step][idx]
    else:
      values = self.values[idx]
    return np.mean(values)

