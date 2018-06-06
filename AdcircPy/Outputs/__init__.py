from AdcircPy.Surface import Surface
from AdcircPy.Outputs import _Outputs

# Factory class
class Outputs(object):
    @staticmethod
    def read_outputs(path, **kwargs):
        return _Outputs.read_outputs(path, **kwargs)

# Generic outputs classes
class SurfaceTimeseries(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)

    def make_animation(self, **kwargs):
        return _SurfaceTimeseries.make_animation(self, **kwargs)

class StationTimeseries(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)

class SurfaceExtrema(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)

# Specific output classes derived from generic output classes
class Maxele(SurfaceExtrema):
    def __init__(self, **kwargs):
        SurfaceExtrema.__init__(self, **kwargs)
        self.datum = kwargs.pop("datum", "MSL")