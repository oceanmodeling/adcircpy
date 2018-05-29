from AdcircPy.Datum import _Datum


class Datum(object):
    def __init__(self, **kwargs):
        self.epsg                = kwargs.pop("epsg", 4326)
        self.datum               = kwargs.pop("datum", "MSL")

    def msl_to_navd88(self, datum_grid):
    	_Datum.msl_to_navd88(self, datum_grid)

    def navd88_to_msl(self):
    	_Datum.navd88_to_msl(self)