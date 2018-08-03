from AdcircPy.Mesh import _Mesh
from AdcircPy.Mesh import _fort13
from AdcircPy.Mesh import _fort14
from AdcircPy.core import Surface


class Mesh(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)
        self.datum       = kwargs.pop("datum", "MSL")
        self.description = kwargs.pop("description", None)
        self.fort13      = kwargs.pop("fort13", None)
        # self.NodalAttributes = NodalAttributes(**kwargs)

    @staticmethod
    def init_from_files(fort14, datum='MSL', epsg=4326, fort13=None):
        return _Mesh.init_from_files(fort14, datum, epsg, fort13)

    def make_plot(self, surface='bathy', **kwargs):
        return _Mesh.make_plot(self, surface, **kwargs)
        # return _Mesh.plot_bathy(self, surface, **kwargs)

    def interpolate_DEM(self, DEM, **kwargs):
        _Mesh.interpolate_DEM(self, DEM, **kwargs)
    
    def write_fort14(self, path):
        _fort14.write_fort14(self, path)



# class NodalAttributes(object):
#     def __init__(self, **kwargs):
#         self.primitive_weighting_in_continuity_equation = kwargs.pop("primitive_weighting_in_continuity_equation", None)
#         self.surface_submergence_state = kwargs.pop("surface_submergence_state", None)
#         self.quadratic_friction_coefficient_at_sea_floor = kwargs.pop("quadratic_friction_coefficient_at_sea_floor", None) 
#         self.surface_directional_effective_roughness_length = kwargs.pop("surface_directional_effective_roughness_length", None)
#         self.surface_canopy_coefficient = kwargs.pop("surface_canopy_coefficient", None)
#         self.bridge_pilings_friction_paramenters = kwargs.pop("bridge_pilings_friction_paramenters", None)
#         self.mannings_n_at_sea_floor = kwargs.pop("mannings_n_at_sea_floor", None)
#         self.chezy_friction_coefficient_at_sea_floor = kwargs.pop("chezy_friction_coefficient_at_sea_floor", None)
#         self.sea_surface_height_above_geoid = kwargs.pop("sea_surface_height_above_geoid", None)
#         self.bottom_roughness_length = kwargs.pop("bottom_roughness_length", None)
#         self.wave_refraction_in_swan = kwargs.pop("wave_refraction_in_swan", None)
#         self.average_horizontal_eddy_viscosity_in_sea_water_wrt_depth = kwargs.pop("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth", None)
#         self.elemental_slope_limiter = kwargs.pop("elemental_slope_limiter", None)
#         self.advection_state = kwargs.pop("advection_state", None)
#         self.initial_river_elevation = kwargs.pop("initial_river_elevation", None)