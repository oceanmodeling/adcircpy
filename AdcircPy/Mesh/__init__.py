from AdcircPy.Mesh import _Mesh
from AdcircPy.Mesh import _fort13
from AdcircPy.Mesh import _fort14
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid

class AdcircMesh(UnstructuredGrid):
  def __init__(self, **kwargs):
    UnstructuredGrid.__init__(self, **kwargs)
    self.datum       = kwargs.pop('datum', None)
    self.description = kwargs.pop('description', None)
    self.fort13      = kwargs.pop('fort13', None)
    if ~isinstance(self.fort13, fort13) and self.fort13 is not None:
      self.fort13    = fort13.from_fort13(self.fort13)
    self.fort15      = kwargs.pop('fort15', None)
    if ~isinstance(self.fort15, fort15) and self.fort15 is not None:
      self.fort15    = fort15.from_fort15(self.fort13)

  @staticmethod
  def from_fort14(fort14, datum='MSL', epsg=4326, fort13=None):
    return _Mesh._from_fort14(fort14, datum=datum, epsg=epsg, fort13=fort13)

  def make_plot(self, surface='bathy', **kwargs):
    return _Mesh.make_plot(self, surface, **kwargs)
    # return _Mesh.plot_bathy(self, surface, **kwargs)

  def interpolate_DEM(self, DEM, **kwargs):
    _Mesh.interpolate_DEM(self, DEM, **kwargs)
  
  def write_fort14(self, path):
    _fort14.write_fort14(self, path)


class fort13(object):
  def __init__(self, **kwargs):
    self.primitive_weighting_in_continuity_equation = kwargs.pop("primitive_weighting_in_continuity_equation", None)
    self.surface_submergence_state = kwargs.pop("surface_submergence_state", None)
    self.quadratic_friction_coefficient_at_sea_floor = kwargs.pop("quadratic_friction_coefficient_at_sea_floor", None) 
    self.surface_directional_effective_roughness_length = kwargs.pop("surface_directional_effective_roughness_length", None)
    self.surface_canopy_coefficient = kwargs.pop("surface_canopy_coefficient", None)
    self.bridge_pilings_friction_paramenters = kwargs.pop("bridge_pilings_friction_paramenters", None)
    self.mannings_n_at_sea_floor = kwargs.pop("mannings_n_at_sea_floor", None)
    self.chezy_friction_coefficient_at_sea_floor = kwargs.pop("chezy_friction_coefficient_at_sea_floor", None)
    self.sea_surface_height_above_geoid = kwargs.pop("sea_surface_height_above_geoid", None)
    self.bottom_roughness_length = kwargs.pop("bottom_roughness_length", None)
    self.wave_refraction_in_swan = kwargs.pop("wave_refraction_in_swan", None)
    self.average_horizontal_eddy_viscosity_in_sea_water_wrt_depth = kwargs.pop("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth", None)
    self.elemental_slope_limiter = kwargs.pop("elemental_slope_limiter", None)
    self.advection_state = kwargs.pop("advection_state", None)
    self.initial_river_elevation = kwargs.pop("initial_river_elevation", None)

  @staticmethod
  def from_fort13(path):
    return _fort13._from_fort13(path)

class fort15(object):
  def __init__(self, **kwargs):
    self.RUNDES           = kwargs.pop("RUNDES", None)
    self.RUNID            = kwargs.pop("RUNID", None)
    self.NFOVER           = kwargs.pop("NFOVER", 0)
    self.NABOUT           = kwargs.pop("NABOUT", 1)
    self.NSCREEN          = kwargs.pop("NSCREEN", 100)
    self.IHOT             = kwargs.pop("IHOT", 0)
    self.ICS              = kwargs.pop("ICS", 0)
    self.IHOT             = kwargs.pop("IHOT", 0)
    self.ICS              = kwargs.pop("IHOT", 0)
    self.IM               = kwargs.pop("IHOT", 0)
    # self.IDEN 
    # self.NOLIBF
    # self.NOLIFA
    # self.NOLICA
    # self.NOLICAT
    self.reference_date   = kwargs.pop("reference_date", None)
    self.reference_time   = kwargs.pop("reference_time", None)
    self.RNDAY            = kwargs.pop("RNDAY", None)
    self.DRAMP            = kwargs.pop("DRAMP", None)
    self.nodal_attributes = kwargs.pop("nodal_attributes", None)
    self.constituent_list = kwargs.pop("constituent_list", None)
    self.boundary_TPXO    = kwargs.pop("boundary_TPXO", None)
    self._mode            = kwargs.pop("mode", None)
    self.eta_stations     = kwargs.pop("eta_stations", None)
    self.met_stations     = kwargs.pop("met_stations", None)
    self.vel_stations     = kwargs.pop("vel_stations", None)
    self.eta_start        = kwargs.pop("eta_start", None)
    self.eta_stop         = kwargs.pop("eta_stop", None)
    self.eta_steps        = kwargs.pop("eta_steps", None)
    self.vel_start        = kwargs.pop("vel_start", None)
    self.vel_stop         = kwargs.pop("vel_stop", None)
    self.vel_steps        = kwargs.pop("vel_steps", None)
    self.hotstart_steps   = kwargs.pop("hotstart_steps", None)
    self.output_filepath  = kwargs.pop("output_filepath", None)

  @staticmethod
  def from_fort15(path):
    return _fort15._from_fort15(path)



  # def generate_forcing_from_TPXO(self, Mesh):
  #     return _Fort15.generate_forcing_from_TPXO(self, Mesh)
