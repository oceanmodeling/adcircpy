from collections import OrderedDict
from scipy.spatial import KDTree
from AdcircPy.Model import UnstructuredMesh as _UnstructuredMesh
from AdcircPy.Model import AdcircMesh as _AdcircMesh
from AdcircPy.Model import _AdcircRun as _AdcircRun__AdcircRun
from AdcircPy.Model import _fort14
from AdcircPy.Model import _fort13


class UnstructuredMesh(object):   
  def __init__(self, x, y, elements, values,
                                     nodeID=None,
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
    return _UnstructuredMesh.get_mean_value(self, **kwargs)

  def get_contours(self, levels, **kwargs):
    """   """
    return _UnstructuredMesh.get_contours(self, levels, **kwargs)

  def get_values_under_Path(self, Path, **kwargs):
    """
    Input is a matplotlib.path.Path instance.
    Returns array of values under specified path.
    """
    return _UnstructuredMesh.get_values_under_Path(self, Path, **kwargs)

  def make_plot(self, **kwargs):
    """  """
    return _UnstructuredMesh.make_plot(self, **kwargs)
  
  def get_values_at_xy(self, x, y, **kwargs):
    """
    Returns numpy array of values at coordinates or list of coordinates give by lon lat
    """
    return _UnstructuredMesh.get_values_at_xy(self, x, y, **kwargs)

  def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
    """  """
    return _UnstructuredMesh.rasterize_to_geoTransform(self, geoTransform, shape, **kwargs)

  def get_xyz(self, extent=None, **kwargs):
    """
    Returns numpy array of dimension [D,3] with columns <x, y, values>
    """
    return _UnstructuredMesh.get_xyz(self, **kwargs)
  
  def get_xy(self, extent=None, epsg=None):
    """
    Returns numpy array of dimension [D,2] with columns <x, y>
    """
    return _UnstructuredMesh.get_xy(self, extent, epsg)

  def get_elements_from_exent(self, extent):
    """ """
    return _UnstructuredMesh.get_elements_from_extent(self, extent)

  def get_finite_volume(self, index):
    """ """
    return _UnstructuredMesh.get_finite_volume(self, index)
    
  def get_Path_from_element(self, element):
    """ """
    return _UnstructuredMesh.get_Path_from_element(self, element)
  
  def get_extent(self, **kwargs):
    """ """
    return _UnstructuredMesh.get_extent(self, **kwargs)

  def get_element_containing_coord(self, coord):
    return _UnstructuredMesh.get_element_containing_coord(self, coord)

  def get_finite_volume_Path_list(self, index):
    return _UnstructuredMesh.get_finite_volume_Path_list(self, index)

  def get_finite_volume_element_list(self, index):
    return _UnstructuredMesh.get_finite_volume_element_list(self, index)
  
  def transform_to_epsg(self, epsg):
    """ """
    _UnstructuredMesh.transform_to_epsg(self, epsg)

  def plot_trimesh(self, **kwargs):
    """ """
    return _UnstructuredMesh.plot_trimesh(self, **kwargs) 
  
  def plot_outerBoundary(self):
    return _UnstructuredMesh.plot_outerBoundary(self)
  
  def _init_fig(self, axes=None, extent=None, title=None, epsg=None):
    """    """
    _UnstructuredMesh._init_fig(self, axes, extent, title, epsg)

  def _init_cbar(self, cmap, vmin, vmax):
    """ """
    _UnstructuredMesh._init_cbar(self, cmap, vmin, vmax)

  def _init_Tri(self):
    _UnstructuredMesh._init_Tri(self)

  def _init_KDTree(self):
    self.KDTree = KDTree(self.get_xy())

  def _init_datum_grid(self):
    _UnstructuredMesh._init_datum_grid(self)

  def _get_finite_volume_interp(self, point_idx):
    """
    Fetches the Path for the DEM interpolator.
    """
    return _UnstructuredMesh._get_finite_volume_interp(self, point_idx)

  def _get_extent_idx(self, extent, epsg, **kwargs):
    """
    Finds the indices of the mesh nodes inside a bounding box.
    kwargs:
        extent : list of the form [min_x, max_x, min_y, max_y]
    """
    return _UnstructuredMesh._get_extent_idx(self, extent, epsg, **kwargs)

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
    return _UnstructuredMesh.get_difference(self, other)

class AdcircMesh(UnstructuredMesh):
  def __init__(self, x, y, elements, values,
                                     fort13=None,
                                     fort15=None,
                                     description=None,
                                     nodeID=None,
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
    UnstructuredMesh.__init__(self, x, y, elements, values, nodeID, elementID,
                              datum, epsg, datum_grid, ocean_boundaries,
                              land_boundaries, inner_boundaries, weir_boundaries,
                              inflow_boundaries, outflow_boundaries, culvert_boundaries)
    self._init_fort13(fort13)

  @classmethod
  def from_fort14(cls, fort14, datum='MSL', epsg=4326, fort13=None, fort15=None, datum_grid=None):
    return _AdcircMesh.from_fort14(cls, fort14, datum, epsg, fort13, fort15, datum_grid)

  def generate_tidal_run(self, start_time, end_time, constituents=None, **kwargs):
    return TidalRun(self, start_time, end_time, constituents, **kwargs)

  def make_plot(self, **kwargs):
    return _AdcircMesh.make_plot(self, **kwargs)

  def interpolate_DEM(self, DEM, **kwargs):
    _AdcircMesh.interpolate_DEM(self, DEM, **kwargs)
  
  def write_fort14(self, path):
    _fort14.write_fort14(self, path)

  def _init_fort13(self, path):
    _fort13._init_fort13(self, path)


class _AdcircRun(object):

  def __init__(self, start_time, end_time, constituents, spinup_date=None, **kwargs):
    self.init_TidalForcing(start_time, end_time, constituents, spinup_date)
    self.init_TPXO()
    self.init_fort15(**kwargs)

  def init_TidalForcing(self, start_time, end_time, constituents, spinup_date=None):
    _AdcircRun__AdcircRun.init_TidalForcing(self, start_time, end_time, constituents, spinup_date)

  def init_TPXO(self):
    _AdcircRun__AdcircRun.init_TPXO(self)

  def init_fort15(self):
    _AdcircRun__AdcircRun.init_fort15(self)

class TidalRun(_AdcircRun):
  def __init__(self, AdcircMesh, start_time, end_time, constituents, **kwargs):
    self.AdcircMesh = AdcircMesh
    super(TidalRun, self).__init__(start_time, end_time, constituents, **kwargs)
    
  @classmethod
  def from_fort14(self, fort14_path, start_time, end_time):
    pass

  def dump(self, directory):
    _TidalRun.dump(self, directory)


class StormRun(object):
  def __init__(self, AdcircMesh, **kwargs):
    self.AdcircMesh = AdcircMesh



