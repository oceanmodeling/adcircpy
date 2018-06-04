from AdcircPy.Surface import Surface
from AdcircPy.Datum import _DatumGrid


class DatumGrid(object):

    def __init__(self, **kwargs):
        self.x = kwargs.pop("x", None)
        self.y = kwargs.pop("y", None)
        self.values = kwargs.pop("values", None)
        self._type = kwargs.pop("_type", None)

    @staticmethod
    def msl_to_navd88(Datum_grid):
        return _DatumGrid.msl_to_navd88(Datum_grid)

    def convert(self, Mesh, **kwargs):
        return _DatumGrid.convert(self, Mesh, **kwargs)


class VDatum(object):
    """
    This class calls the VDatum Java executables over a subprocess pipe.
    You must specify the directory your copy of VDatum using the argument vdatumdir.
    VDatum requires to have Java 8 installed. Other version of Java will not work.
    """
    @staticmethod
    def convert(Mesh, **kwargs):
        return _vdatum.convert(Mesh, **kwargs)
