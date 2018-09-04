from collections import OrderedDict
from AdcircPy.Mesh import _Mesh
from AdcircPy.Mesh import _fort13
from AdcircPy.Mesh import _fort14
from AdcircPy.Mesh import _fort15
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid

class UnstructuredGrid(object):   
  def __init__(self, x, y, values, elements, **kwargs):
    self.x           = x
    self.y           = y
    self._values     = values
    self.elements    = elements
    self.nodeID      = kwargs.pop("nodeID", None)
    self.elementID   = kwargs.pop("elementID", None)
    self.datum       = kwargs.pop("datum", None)
    self.epsg        = kwargs.pop("epsg", None)
    self.ocean_boundaries    = kwargs.pop("ocean_boundaries", None)
    self.land_boundaries     = kwargs.pop("land_boundaries", None)
    self.inner_boundaries    = kwargs.pop("inner_boundaries", None)
    self.weir_boundaries     = kwargs.pop("weir_boundaries", None)
    self.inflow_boundaries   = kwargs.pop("inflow_boundaries", None)
    self.outflow_boundaries  = kwargs.pop("outflow_boundaries", None)
    self.culvert_boundaries  = kwargs.pop("culvert_boundaries", None)

  @property
  def values(self):
    return self._values

  @values.setter
  def values(self, values):
    self._values = values

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

  def get_land_boundaries(self, **kwargs):
    return _Boundaries.get_land_boundaries(self, **kwargs)

  def build_outer_polygon(self, **kwargs):
    return _Boundaries.build_outer_polygon(self, **kwargs)

  def build_inner_polygons(self, **kwargs):
    return _Boundaries.build_inner_polygons(self, **kwargs)

  def plot_outerBoundary(self):
    return _Boundaries.plot_outerBoundary(self)

class AdcircMesh(UnstructuredGrid):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, x, y, values, elements,
                              epsg=kwargs.pop('epsg', 4326), **kwargs)
    self.description = kwargs.pop('description', None)
    self._boundary_forcing=OrderedDict()
    self._init_TPXO()
    self._init_fort15(kwargs.pop('fort15', None))
    self._init_fort13(kwargs.pop('fort13', None))

  @staticmethod
  def from_fort14(fort14, datum='MSL', epsg=4326, fort13=None, fort15=None):
    return _Mesh._from_fort14(fort14, datum, epsg, fort13, fort15)

  @property
  def boundary_forcing(self):
    return self._boundary_forcing

  @boundary_forcing.setter 
  def boundary_forcing(self, constituent_list):
    _Mesh._set_boundary_forcing(self, constituent_list) 

  def make_plot(self, surface='bathy', **kwargs):
    return _Mesh.make_plot(self, surface, **kwargs)
    # return _Mesh.plot_bathy(self, surface, **kwargs)

  def interpolate_DEM(self, DEM, **kwargs):
    _Mesh.interpolate_DEM(self, DEM, **kwargs)
  
  def write_fort14(self, path):
    _fort14.write_fort14(self, path)

  def _init_fort15(self, fort15):
    _fort15._init_fort15(self, fort15)

  def _init_fort13(self, fort13):
    _fort13._init_fort13(self, fort13)

  def _init_TPXO(self):
    _Mesh._init_TPXO(self)




class fort13(dict):
  def __init__(self, **kwargs):
    # self.primitive_weighting_in_continuity_equation = kwargs.pop("primitive_weighting_in_continuity_equation", None)
    # self.surface_submergence_state = kwargs.pop("surface_submergence_state", None)
    # self.quadratic_friction_coefficient_at_sea_floor = kwargs.pop("quadratic_friction_coefficient_at_sea_floor", None) 
    # self.surface_directional_effective_roughness_length = kwargs.pop("surface_directional_effective_roughness_length", None)
    # self.surface_canopy_coefficient = kwargs.pop("surface_canopy_coefficient", None)
    # self.bridge_pilings_friction_paramenters = kwargs.pop("bridge_pilings_friction_paramenters", None)
    # self.mannings_n_at_sea_floor = kwargs.pop("mannings_n_at_sea_floor", None)
    # self.chezy_friction_coefficient_at_sea_floor = kwargs.pop("chezy_friction_coefficient_at_sea_floor", None)
    # self.sea_surface_height_above_geoid = kwargs.pop("sea_surface_height_above_geoid", None)
    # self.bottom_roughness_length = kwargs.pop("bottom_roughness_length", None)
    # self.wave_refraction_in_swan = kwargs.pop("wave_refraction_in_swan", None)
    # self.average_horizontal_eddy_viscosity_in_sea_water_wrt_depth = kwargs.pop("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth", None)
    # self.elemental_slope_limiter = kwargs.pop("elemental_slope_limiter", None)
    # self.advection_state = kwargs.pop("advection_state", None)
    # self.initial_river_elevation = kwargs.pop("initial_river_elevation", None)
    dict.__init__(self, **kwargs)

class fort15(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)

  def generate_equilibrium_arguments(self, tpxo_path, start_date, end_date):
    return _fort15._generate_equilibrium_arguments(self, start_date, end_date)
