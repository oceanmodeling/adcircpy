from collections import OrderedDict
from AdcircPy.Mesh import _Mesh
from AdcircPy.Mesh import _fort13
from AdcircPy.Mesh import _fort14
from AdcircPy.Mesh import _fort15
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid

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
