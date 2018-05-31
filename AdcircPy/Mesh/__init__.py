from AdcircPy.Mesh import _Trimesh
from AdcircPy.Mesh import _Boundaries
from AdcircPy.Mesh import _SurfaceDifference
from AdcircPy.Mesh import _Mesh
from AdcircPy.Mesh import _fort14
from AdcircPy.Mesh import _raster
from AdcircPy.Datum import Datum
    
class Trimesh(object):
    def __init__(self, **kwargs):
        self.x         = kwargs.pop("x", None)
        self.y         = kwargs.pop("y", None)
        self.nodeID    = kwargs.pop("nodeID", None)
        self.elements  = kwargs.pop("elements",  None)
        self.elementID = kwargs.pop("elementID", None)

    def plot_trimesh(self, **kwargs):
        return _Trimesh.plot_trimesh(self, **kwargs)    

    def get_xyz(self, extent=None, **kwargs):
        """
        Returns numpy array of dimension [D,3] with columns <x, y, values>
        """
        return _Trimesh.get_xyz(self, **kwargs)
    
    def get_xy(self, extent=None, **kwargs):
        """
        Returns numpy array of dimension [D,2] with columns <x, y>
        """
        return _Trimesh.get_xy(self, extent, **kwargs)

    def get_elements_as_Path_list(self, extent=None):
        return _Trimesh.get_element_as_Path_list(self, extent)

    def get_elements_surrounding_node(self, node_index):
        return _Trimesh.get_elements_surrounding_node(self, node_index)
      
    def get_extent_idx(self, extent):
        """
        Finds the indices of the mesh nodes inside a bounding box.
        kwargs:
            extent : list of the form [min_x, max_x, min_y, max_y]
        """
        return _Trimesh.get_extent_idx(self, extent)

    def get_extent(self):
        return _Trimesh.get_extent(self)


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

class Surface(Trimesh, Boundaries):




    def get_values_under_Path(self, Path, **kwargs):
        """
        Input is a matplotlib.path.Path instance.
        Returns array of values under specified path.
        """
        return _Trimesh.get_values_under_Path(self, Path, **kwargs)

class NodalAttributes(Surface):
    def __init__(self, **kwargs):
        pass

class Mesh(Surface, Datum):
    def __init__(self, **kwargs):
        Trimesh.__init__(self, **kwargs)
        Boundaries.__init__(self, **kwargs)
        self.values = kwargs.pop("values", None)
        self.description = kwargs.pop("description", None)

    @staticmethod
    def init_from_fort14(fort14, datum='MSL', epsg=4326):
        return _Mesh.init_from_fort14(fort14, datum, epsg)
      
    def reproject(self, epsg):
        return _Mesh.reproject(self, epsg)

    def make_plot(self, **kwargs):
        return _Mesh.plot_bathy(self, **kwargs)

    def interpolate_DEM(self, DEM, **kwargs):
        _Mesh.interpolate_DEM(self, DEM, **kwargs)
    
    def get_values_at_lonlat(self, lon, lat, **kwargs):
        """
        Returns numpy array of values at coordinates or list of coordinates give by lon lat
        """
        return _Mesh.get_values_at_lonlat(self, lon, lat, **kwargs)
        
    def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
        return _raster.rasterize_to_geoTransform(self, geoTransform, shape, **kwargs)
    
    def write_fort14(self, path):
        _fort14.write_fort14(self, path)
    
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
        return _SurfaceDifference.get_difference(self, other)
        
class SurfaceDifference(ScalarSurface):
    def make_plot(self, **kwargs):
        return _SurfaceDifference.plot_diff(self, **kwargs)

