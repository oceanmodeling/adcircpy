from AdcircPy.Surface import Surface
from AdcircPy.Outputs import _Outputs
from AdcircPy.Outputs import _ElevationStations
from AdcircPy.Outputs import _HarmonicConstituents

# Factory class
class Outputs(object):
    def __init__(self, path, fort14=None, datum=None, epsg=None):
        self._path   = path
        self._fort14 = fort14
        self.datum  = datum
        self.epsg   = epsg
   
    @staticmethod
    def read_outputs(path, **kwargs):
        return _Outputs.read_outputs(path, **kwargs)

    def _get_netcdf_output(self):
        return _netcdf._get_output(self)

    def _open_file(self):
        return _Outputs._open_file(self)

    def _read_netcdf(self):
        return _netcdf._read_netcdf(self)

    def _read_ascii(self):
        return _ascii._read_ascii(self)

    def _check_netcdf(self):
        _netcdf._check_netcdf(self)

# Generic outputs classes
class StationTimeseries(object):
    def __init__(self, x, y, zeta, time, station_name, Dataset=None):
        self.x            = x
        self.y            = y
        self.time         = time
        self.station_name = station_name
        self.Dataset      = Dataset

class ElevationStations(StationTimeseries):
    def __init__(self, x, y, zeta, time, station_name, Dataset=None):
        StationTimeseries.__init__(self, x, y, time, station_name, Dataset)
        self.zeta         = zeta

    @staticmethod
    def from_netcdf(path):
        return _ElevationStations._from_netcdf(path)

    def make_plot(self, station, **kwargs):
        return _ElevationStations._make_plot(self, station, **kwargs)

class HarmonicConstituents(StationTimeseries):
    def __init__(self, x, y, zeta, time, station_name, Dataset=None):
        StationTimeseries.__init__(self, x, y, time, station_name, Dataset)
        self.amplitude = amplitude


    @staticmethod
    def from_fort51(path):
        return _HarmonicConstituents._from_fort51(path)


    


class SurfaceTimeseries(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)

    def make_animation(self, **kwargs):
        return _SurfaceTimeseries.make_animation(self, **kwargs)

class SurfaceExtrema(Surface):
    def __init__(self, **kwargs):
        Surface.__init__(self, **kwargs)

# Specific output classes derived from generic output classes
class Maxele(SurfaceExtrema):
    def __init__(self, **kwargs):
        SurfaceExtrema.__init__(self, **kwargs)
        self.datum = kwargs.pop("datum", "MSL")

class _ncOutput(Outputs):
    def __init__(self, path, fort14=None):
        Outputs.__init__(self, path, fort14)
        self.Dataset = Dataset(self.path)
        self._get__type()

    def _get__type(self):
        for key in self.Dataset.variables.keys():
            print(key)
