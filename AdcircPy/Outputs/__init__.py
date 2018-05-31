from AdcircPy.Mesh    import Trimesh, Boundaries, ScalarSurface
from AdcircPy.Datum   import Datum
from AdcircPy.Outputs import _Outputs

class Outputs(object):
    @staticmethod
    def read_outputs(path, **kwargs):
        return _Outputs.read_outputs(path, **kwargs)

        
class ScalarOutput(ScalarSurface):
    def __init__(self, **kwargs):      
        Trimesh.__init__(self, **kwargs)
        Boundaries.__init__(self, **kwargs)
        Datum.__init__(self, **kwargs)
        self.values = kwargs.pop("values", None)

class VectorOutput(object):
    pass
        
        
class Maxele(ScalarOutput):
    
    def __init__(self, **kwargs):      
        Trimesh.__init__(self, **kwargs)
        Boundaries.__init__(self, **kwargs)
        Datum.__init__(self, **kwargs)
        self.values = kwargs.pop("values", None)

    def make_plot(self, **kwargs):
        return Surface.make_plot(self, **kwargs)

class Velocity(Trimesh, Boundaries):
    pass

class TidalAmplitudes(Trimesh, Boundaries, Datum):
    """docstring for ClassName"""
    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg
        
class Elevation(Trimesh, Boundaries, Datum):
    """docstring for ClassName"""
    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg

