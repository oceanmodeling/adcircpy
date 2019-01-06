import numpy as np
from matplotlib.tri import Triangulation
from AdcircPy.Mesh import Trimesh

class TrimeshSurface(Trimesh):
  def __init__(self, x, y, elements, values, epsg, nodeID=None, elementID=None, vertical_datum=None, datum_grid=None, **Boundaries):
    super(TrimeshSurface, self).__init__(x, y, elements, epsg, nodeID, elementID, **Boundaries)
    self._init_values(values)
    self._init_Tri()
    self._vertical_datum = vertical_datum
    self._init_datum_grid(datum_grid)

  @property
  def values(self):
    return self._values

  @property
  def vertical_datum(self):
    return self._vertical_datum

  @property
  def Tri(self):
    return self._Tri
  

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

  def get_values_under_Path(self, Path, **kwargs):
    xin = self.x
    yin = self.y
    if len(self.values.shape) == 3:
      zin = self.values[:,0,timestep]
    else:
      zin = self.values
    idx, = np.where(np.logical_and(
        np.logical_and(xin>=np.min(Path.vertices[:,0]), xin<=np.max(Path.vertices[:,0])),
        np.logical_and(yin>=np.min(Path.vertices[:,1]), yin<=np.max(Path.vertices[:,1]))))
    return griddata((xin[idx], yin[idx]), zin[idx], (Path.vertices[:,0],Path.vertices[:,1]))

  def make_plot(self, extent=None, vmin=None, vmax=None, colors=256, cmap='jet', axes=None, figsize=None, title=None, show=False, Trimesh=False, color='black', linewidth=1.0, alpha=0.5):
    self._init_fig(axes, figsize)
    self._init_vmin(vmin)
    self._init_vmax(vmax)
    self._init_cmap_levels(cmap, colors)
    self._axes.tricontourf(self.Tri, self.values, levels=self._levels)
    if Trimesh==True:
      self._axes.triplot(self.x, self.y, self.elements, color=color, linewidth=linewidth, alpha=alpha)
    if extent is None:
      extent=self.get_extent()
    self._auto_scaling(extent)
    self._init_title(title)
    self._init_cbar('neither')
    self._auto_show(show)

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

  def _init_values(self, values):
    self._values = np.asarray(values)
    if self.values.size != self.x.size != self.y.size:
      raise Exception('Provided values array does not match ')

  def _init_Tri(self):
    if np.ma.is_masked(self._values):
      trimask = np.any(self._values.mask[self.elements], axis=1)
      self._Tri = Triangulation(self.x, self.y, self.elements, trimask)
    else:
      self._Tri = Triangulation(self.x, self.y, self.elements)

  def _fv_interp_node(self, idx, radius=None):
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

  # def _init_cmap(self, cmap):
  #   if cmap is None:
  #     self._cmap = self._plt.cm.get_cmap('viridis')
  #   else:
  #     self._cmap = self._plt.cm.get_cmap(cmap)
  
  # def _init_levels(self, colors):
  #   self._levels=np.linspace(self._vmin, self._vmax, colors)

  def _init_datum_grid(self, datum_grid):
    """
    The way this is implemented is not self consistent.
    Vertical datum grids have to be handled with care until
    some standards can be set. Meshes should actually be referenced
    to an equipotential surface, and not to space variable datums.
    This is because at time t0 in the model run, the water surface 
    elevation must sit on an equipotential surface. Datums derived from
    tidal measurements are not generally equipotentials.
    """
    if datum_grid is not None:
      # if self.datum is not None:
      with open(datum_grid, 'r') as f: 
        f.readline().rstrip()
        line = f.readline().strip('\n').split(': ')
        line = line.pop().split('to')
        original_vdatum = line[0].split(':')[1].strip(' ')
        target_vdatum = line[1].split(':')[1].strip(' ')
        if original_vdatum == self._vertical_datum:
          NP = int(f.readline().split()[1])
          values = list()
          for k in range(NP):
            values.append(float(f.readline().split()[3]))
          self.original_mesh_values = self.values
          self.original_mesh_datum = self._vertical_datum
          self._values += np.asarray(values)
          self._vertical_datum = target_vdatum

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
