from AdcircPy.Mesh import _Mesh
from AdcircPy.Surface import Surface
    
class Mesh(Surface):
    def __init__(self, **kwargs):
        Trimesh.__init__(self, **kwargs)
        Boundaries.__init__(self, **kwargs)
        self.epsg        = kwargs.pop("epsg", None)
        self.datum       = kwargs.pop("datum", "MSL")
        self.description = kwargs.pop("description", None)

    @staticmethod
    def init_from_fort14(fort14, datum='MSL', epsg=4326):
        return _Mesh.init_from_fort14(fort14, datum, epsg)

    def make_plot(self, **kwargs):
        return _Mesh.plot_bathy(self, **kwargs)

    def interpolate_DEM(self, DEM, **kwargs):
        _Mesh.interpolate_DEM(self, DEM, **kwargs)
    
    def write_fort14(self, path):
        _Mesh.write_fort14(self, path)

class NodalAttributes(Surface):
    def __init__(self, **kwargs):
        pass
   


