from collections import defaultdict
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from scipy.interpolate import griddata
from scipy.spatial import KDTree
import pyproj
from haversine import haversine


class UnstructuredMesh(object):   
  def __init__(self, x, y, elements, values, nodeID=None,
                                     elementID=None,
                                     datum=None,
                                     epsg=None,
                                     datum_grid=None,
                                     ocean_boundaries=None,
                                     land_boundaries=None,
                                     inner_boundaries=None,
                                     weir_boundaries=None,
                                     inflow_boundaries=None,
                                     outflow_boundaries=None,
                                     culvert_boundaries=None):
    self._x           = x
    self._y           = y
    self._values      = values
    self.elements     = elements
    self.nodeID       = nodeID
    self.elementID    = elementID
    self.datum        = datum
    self.epsg         = epsg
    self._datum_grid  = datum_grid
    self.ocean_boundaries    = ocean_boundaries
    self.land_boundaries     = land_boundaries
    self.inner_boundaries    = inner_boundaries
    self.weir_boundaries     = weir_boundaries
    self.inflow_boundaries   = inflow_boundaries
    self.outflow_boundaries  = outflow_boundaries
    self.culvert_boundaries  = culvert_boundaries
    self._init_datum_grid()
    self._init_Tri()
    self._init_KDTree()

  @property
  def values(self):
    return self._values

  @property
  def x(self):
    return self._x
  
  @property
  def y(self):
    return self._y
    
  def get_mean_value(self, **kwargs):
    """   """
    epsg   = kwargs.pop("epsg", self.epsg)
    extent = kwargs.pop("extent", self.get_extent(epsg=epsg))
    step   = kwargs.pop("step", 0)
    idx = self.get_extent_idx(extent, epsg)
    if isinstance(self.values, list):
        values = self.values[step]
    else:
        values = self.values
    return np.nanmean(values[idx])


  def get_xyz(self, extent=None, epsg=None):
    if extent is None:
      extent = self.get_extent()
    if epsg is None:
      epsg = self.epsg
    idx = self._get_extent_idx(extent, epsg)
    return np.vstack((self.x[idx], self.y[idx], self.values[idx])).T

  def get_xy(self, extent=None, epsg=None):
    if extent is None:
      extent = self.get_extent()
    if epsg is None:
      epsg = self.epsg
    idx = self._get_extent_idx(extent, epsg)
    return np.vstack((self.x[idx], self.y[idx])).T

  def transform_to_epsg(self, epsg):
    self_proj = pyproj.Proj(init="epsg:{}".format(self.epsg))
    target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
    x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    self.x = np.asarray(x).flatten()
    self.y = np.asarray(y).flatten()
    self.epsg = epsg

  def _init_Tri(self):
    if np.ma.is_masked(self._values):
      trimask = np.any(self._values.mask[self.elements], axis=1)
      self._Tri = Triangulation(self.x, self.y, self.elements, trimask)
    else:
      self._Tri = Triangulation(self.x, self.y, self.elements)

  def _update_Tri(self):
    self._init_Tri()

  def _init_KDTree(self):
    self.KDTree = KDTree(self.get_xy())

  def make_plot(self, extent=None, epsg=None, axes=None, title=None, show=False, **kwargs):
    # print(self.va)
    total_colors=256
    self._init_fig(axes, extent, title, epsg)
    self._vmin = kwargs.pop("vmin", np.min(self._values))
    self._vmax = kwargs.pop("vmax", np.max(self._values))
    self._cmap = kwargs.pop("cmap", "jet")
    self._levels = kwargs.pop("levels", np.linspace(self._vmin, self._vmax, total_colors))
    self._axes.tricontourf(self._Tri, self._values, levels=self._levels, cmap=self._cmap, extend='both')
    self._init_cbar(self._cmap, self._vmin, self._vmax)
    cbar_label=None
    if cbar_label is not None:
      self._cbar.set_label(cbar_label)
    self._cbar.set_ticks([self._vmin,
                    self._vmin+(1./4.)*(self._vmax-self._vmin),
                    self._vmin+(1./2.)*(self._vmax-self._vmin),
                    self._vmin+(3./4.)*(self._vmax-self._vmin),
                    self._vmax])
    self._cbar.set_ticklabels([np.around(self._vmin, 2), 
                            np.around(self._vmin+(1./4.)*(self._vmax-self._vmin), 2),
                            np.around(self._vmin+(1./2.)*(self._vmax-self._vmin), 2),
                            np.around(self._vmin+(3./4.)*(self._vmax-self._vmin), 2),
                            np.around(self._vmax, 2)])
    
    if show == True:
      plt.show()

    return self._axes

  def plot_velocity(self, **kwargs):
      raise NotImplementedError("Coming soon!")
      start_timestep = kwargs.pop('start_timestep', 0)
      stop_timestep  = kwargs.pop('stop_timestep', len(self.time))
      if axes is None:                
          fig = plt.figure()
          axes  = fig.add_subplot(111)
      ax = axes.quiver(self.x, self.y, self.u[0,:], self.v[:,0], vmin=vmin, vmax=vmax, **kwargs)
      axes.axis('scaled')
      def update(i):
         ax.set_array(self.Dataset['zs'][i,:-1,:-1].ravel())
         return ax
      
      anim = FuncAnimation(fig, update, frames=np.arange(start_timestep+1, stop_timestep), interval=interval)
      if colorbar==True:
          plt.colorbar(ax)
      if show is True:
          plt.show()
      return anim

  def __sub__(self, other, step=0):
    if isinstance(self.values, list):
      self_values = self.values[step]
    else:
      self_values = self.values

    if isinstance(other.values, list):
      other_values = other.values[step]
    else:
      other_values = other.values

    if np.ma.is_masked(self_values):
      self_values = np.ma.filled(self_values, 0.)

    if np.ma.is_masked(other.values):
      other_values = np.ma.filled(other_values, 0.)
    values = self_values - other_values
    # values = np.ma.masked_equal(values, 0.)
    kwargs = self.get_dict()
    kwargs['values'] = values
    return Surface.SurfaceDifference(**kwargs)


  def plot_diff(self, extent=None, epsg=None, axes=None, vmin=None, vmax=None, title=None, **kwargs):
    axes, idx = _fig.init_fig(self, axes, extent, title, epsg)
    if vmin is None:
      vmin = np.min(self.values[idx])
    if vmax is None:
      np.max(self.values[idx])

    cmap = plt.get_cmap(kwargs.pop("cmap", "seismic"))
    levels = kwargs.pop("levels", np.linspace(np.min(self.values[idx]), np.max(self.values[idx]), 256))
    norm = _fig.FixPointNormalize(sealevel=0, vmax=vmax, vmin=vmin, col_val=0.5)
    if np.ma.is_masked(self.values):
      trimask = np.any(self.values.mask[self.elements], axis=1)
      Tri = Triangulation(self.x, self.y, self.elements, trimask)
      axes.tricontourf(Tri, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    else:
      axes.tricontourf(self.x, self.y, self.elements, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    cbar = self._init_cbar(axes, cmap, vmin, vmax)
    cbar.set_ticks([vmin, vmin + 0.5*(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    cbar.set_label(r'elevation [$\Delta$ m]')
    return axes

  def get_contours(self, levels, **kwargs):
      epsg = kwargs.pop("epsg", self.epsg)

      if epsg != self.epsg:
          self_proj   = pyproj.Proj(init='epsg:{}'.format(self.epsg))
          target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
          x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)

      else:
          x = self.x
          y = self.y

      if np.ma.is_masked(self.values):
          mask = np.where(self.values.mask)
          trimask = np.any(np.in1d(self.elements, mask).reshape(-1, 3), axis=1)
          Tri = Triangulation(x, y, self.elements, trimask)
      else:
          Tri = Triangulation(x, y, self.elements)
      fig = plt.figure()
      axes = fig.add_subplot(111)
      ax = axes.tricontour(Tri, self.values, levels=levels)
      plt.close(fig)

      contours = defaultdict(list)
      for i, LineCollection in enumerate(ax.collections):
          for Path in LineCollection.get_paths():
              contours[ax.levels[i]].append(Path)
      return dict(contours)

  def get_values_at_xy(self, x, y, step=0, method='linear'):
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
          return griddata((self.x,self.y),values,(x,y), method=method)
      else:
          idx = np.where(~np.isin(values, 0.0))
          return griddata((self.x[idx], self.y[idx]), values[idx], (x, y), method='nearest')

  def get_finite_volume_element_list(self, index):
      return self.elements[np.where(np.any(np.isin(self.elements, index), axis=1))[0]]

  def get_finite_volume_Path_list(self, index):
      return [self.get_Path_from_element(element) for element in self.get_finite_volume_element_list(index)]

  def get_element_containing_coord(self, coord):
      x = coord[0]; y = coord[1]
      distance, node_idx = self.KDTree.query([x,y])
      elements = self.get_finite_volume_Path_list(node_idx)
      for i, element in enumerate(elements):
          if element.contains_point((x,y)):
              return self.get_finite_volume_element_list(node_idx)[i]

  def _init_fig(self, axes=None, extent=None, title=None, epsg=None):
    self.plt = plt
    if axes is None:                
      fig = self.plt.figure()
      axes  = fig.add_subplot(111)
    if title is not None:
      axes.set_title(title)
    if extent is None:
      extent = self.get_extent()
    if epsg is None:
      epsg = self.epsg
      idx  = self._get_extent_idx(extent, epsg)
      axes.axis('scaled')
      axes.axis(extent) 
    self._axes=axes
    self._idx=idx

  def _init_cbar(self, cmap, vmin, vmax):
    divider = make_axes_locatable(self._axes)
    cax     = divider.append_axes("bottom", size="2%", pad=0.5)
    mappable = ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(vmin, vmax)
    self._cbar = plt.colorbar(mappable, cax=cax, extend='both', orientation='horizontal')

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

  def _get_extent_idx(self, extent, epsg, **kwargs):
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
      self._init_fig(axes, extent, title)
      self._axes.triplot(self.x, self.y, self.elements, color=color, linewidth=linewidth, alpha=alpha)


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

  def _get_finite_volume_interp(self, idx, radius=None):
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

  def _init_datum_grid(self):
    """
    The way this is implemented is not self consistent.
    Vertical datum grids have to be handled with care until
    some standards can be set. Meshes should actually be referenced
    to an equipotential surface, and not to space variable datums.
    This is because at time t0 in the model run, the water surface 
    elevation must sit on an equipotential surface. Datums derived from
    tidal measurements are not generally equipotentials.
    """
    if self._datum_grid is not None:
      # if self.datum is not None:
      with open(self._datum_grid, 'r') as f: 
        f.readline().rstrip()
        line = f.readline().strip('\n').split(': ')
        line = line.pop().split('to')
        original_vdatum = line[0].split(':')[1].strip(' ')
        target_vdatum = line[1].split(':')[1].strip(' ')
        if original_vdatum == self.datum:
          NP = int(f.readline().split()[1])
          values = list()
          for k in range(NP):
            values.append(float(f.readline().split()[3]))
          self.original_mesh_values = self.values
          self.original_mesh_datum = self.datum
          self._values += np.asarray(values)
          self.datum = target_vdatum

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
