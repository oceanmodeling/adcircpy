from matplotlib.colors import Normalize
from AdcircPy.core.UnstructuredGrid import _UnstructuredGrid, _Boundaries

class Boundaries(object):
  def __init__(self, **kwargs):
    self.ocean_boundaries    = kwargs.pop("ocean_boundaries", None)
    self.land_boundaries     = kwargs.pop("land_boundaries", None)
    self.inner_boundaries    = kwargs.pop("inner_boundaries", None)
    self.weir_boundaries     = kwargs.pop("weir_boundaries", None)
    self.inflow_boundaries   = kwargs.pop("inflow_boundaries", None)
    self.outflow_boundaries  = kwargs.pop("outflow_boundaries", None)
    self.culvert_boundaries  = kwargs.pop("culvert_boundaries", None)

  def get_land_boundaries(self, **kwargs):
    return _Boundaries.get_land_boundaries(self, **kwargs)

  def build_outer_polygon(self, **kwargs):
    return _Boundaries.build_outer_polygon(self, **kwargs)

  def build_inner_polygons(self, **kwargs):
    return _Boundaries.build_inner_polygons(self, **kwargs)

  def plot_outerBoundary(self):
    return _Boundaries.plot_outerBoundary(self)

class UnstructuredGrid(object):   
  def __init__(self, x, y, values, elements, **kwargs):
    self.x          = x
    self.y          = y
    self.values = values
    self.elements   = elements
    self.nodeID     = kwargs.pop("nodeID", None)
    self.elementID  = kwargs.pop("elementID", None)
    self.datum      = kwargs.pop("datum", None)
    self.epsg       = kwargs.pop("epsg", None)
    self.Boundaries = Boundaries(**kwargs)

  def __sub__(self, other):
    """
    Used in the case where two different grids are subtracted directly.
    Grids must have the same dimensions.
    Example:
        grid1 = adcpy.read_grid("fort.14.1")
        grid2 = adcpy.read_grid("fort.14.2")
        diff = grid2 - grid1
        diff.make_plot(show=True)
    """
    return _UnstructuredGrid.get_difference(self, other)

  def get_mean_value(self, **kwargs):
    """   """
    return _UnstructuredGrid.get_mean_value(self, **kwargs)

  def get_contours(self, levels, **kwargs):
    """   """
    return _UnstructuredGrid.get_contours(self, levels, **kwargs)

  def get_values_under_Path(self, Path, **kwargs):
    """
    Input is a matplotlib.path.Path instance.
    Returns array of values under specified path.
    """
    return _UnstructuredGrid.get_values_under_Path(self, Path, **kwargs)

  def make_plot(self, **kwargs):
    """  """
    return _UnstructuredGrid.make_plot(self, **kwargs)

  def get_dict(self):
    """  """
    return _UnstructuredGrid.get_dict(self)
  
  def get_values_at_xy(self, x, y, **kwargs):
    """
    Returns numpy array of values at coordinates or list of coordinates give by lon lat
    """
    return _UnstructuredGrid.get_values_at_xy(self, x, y, **kwargs)
      
  def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
    """
    """
    return _UnstructuredGrid.rasterize_to_geoTransform(self, geoTransform, shape, **kwargs)
  
  def _init_fig(self, axes=None, extent=None, title=None, epsg=None):
    """    """
    return _UnstructuredGrid._init_fig(self, axes, extent, title, epsg)

  def _init_cbar(self, axes, cmap, vmin, vmax):
    """ """
    return _UnstructuredGrid._init_cbar(self, axes, cmap, vmin, vmax)

  def _get_finite_volume_interp(self, point_idx):
    """
    Fetches the Path for the DEM interpolator.
    """
    return _UnstructuredGrid._get_finite_volume_interp(self, point_idx)

  def plot_trimesh(self, **kwargs):
    """ """
    return _UnstructuredGrid.plot_trimesh(self, **kwargs)    

  def get_xyz(self, extent=None, **kwargs):
    """
    Returns numpy array of dimension [D,3] with columns <x, y, values>
    """
    return _UnstructuredGrid.get_xyz(self, **kwargs)
  
  def get_xy(self, extent=None, epsg=None):
    """
    Returns numpy array of dimension [D,2] with columns <x, y>
    """
    return _UnstructuredGrid.get_xy(self, extent, epsg)

  def get_elements_from_exent(self, extent):
    """ """
    return _UnstructuredGrid.get_elements_from_extent(self, extent)

  def get_elements_surrounding_index(self, index):
    """ """
    return _UnstructuredGrid.get_elements_surrounding_index(self, index)
    
  def get_Path_from_element_indexes(self, elements):
    """ """
    return _UnstructuredGrid._get_Path_from_element_indexes(self, elements)


  def _get_extent_idx(self, extent, epsg, **kwargs):
    """
    Finds the indices of the mesh nodes inside a bounding box.
    kwargs:
        extent : list of the form [min_x, max_x, min_y, max_y]
    """
    return _UnstructuredGrid._get_extent_idx(self, extent, epsg, **kwargs)

  def get_extent(self, **kwargs):
    """ """
    return _UnstructuredGrid.get_extent(self, **kwargs)

  def transform_to_epsg(self, epsg):
    """ """
    _UnstructuredGrid.transform_to_epsg(self, epsg)


class _SurfaceDifference(UnstructuredGrid):
  def make_plot(self, **kwargs):
    return _UnstructuredGrid.plot_diff(self, **kwargs)

class _FixPointNormalize(Normalize):
  """ 
  Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
  Subclassing Normalize to obtain a colormap with a fixpoint 
  somewhere in the middle of the colormap.
  This may be useful for a `terrain` map, to set the "sea level" 
  to a color in the blue/turquise range.
  """
  def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.5, clip=False):
    # sealevel is the fix point of the colormap (in data units)
    self.sealevel = sealevel
    # col_val is the color value in the range [0,1] that should represent the sealevel.
    self.col_val = col_val
    Normalize.__init__(self, vmin, vmax, clip)

  def __call__(self, value, clip=None):
    x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
    return np.ma.masked_where(value.mask, np.interp(value, x, y))