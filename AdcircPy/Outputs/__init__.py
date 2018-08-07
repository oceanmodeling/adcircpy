from AdcircPy.core.UnstructuredGrid import UnstructuredGrid
from AdcircPy.Outputs import _Outputs
from AdcircPy.Outputs import _ElevationStations
from AdcircPy.Outputs import _HarmonicConstituents

# Factory class
class Outputs(object):
    def __init__(self, path, fort14=None, fort15=None, fort13=None, datum=None, epsg=None):
        self._path   = path
        self._fort14 = fort14
        self._fort15 = fort15
        self._fort13 = fort13
        self.datum   = datum
        self.epsg    = epsg
 
   
    @staticmethod
    def read_outputs(path, **kwargs):
        return _Outputs.read_outputs(path, **kwargs)

    def _get_netcdf_output(self):
        return _netcdf._get_output(self)

    def _open_file(self):
        return _Outputs._open_file(self)

    def _read_netcdf(self):
        return _netcdf._read_netcdf(self)

    def _read_ascii_type(self):
        return _ascii._read_ascii_type(self)

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

class HarmonicConstituents(object):
    def __init__(self, data_dict):
        for constituent in data_dict.keys():
            setattr(self, constituent, data_dict[constituent])

    @staticmethod
    def from_ascii(path, fort14=None, fort15=None, datum=None, epsg=None):
        return _HarmonicConstituents._from_ascii(path, fort14=fort14, fort15=fort15, datum=datum, epsg=epsg)


class HarmonicConstituentTimeseries(HarmonicConstituents):
    def __init__(self, harmonic_constituents):
        HarmonicConstituents.__init__(self, )

    @staticmethod
    def from_ascii(path, fort15=None, station_data_dict=None):
        return _HarmonicConstituentTimeseries._from_ascii(path, fort15=fort15, datum=datum)

class HarmonicConstituentSurface(HarmonicConstituents):
    def __init__(self, harmonic_constituents):
        HarmonicConstituents.__init__(self, )

    @staticmethod
    def from_ascii(path, fort14=None, datum=None, epsg=None):
        return _HarmonicConstituents._from_ascii(path, fort14=fort14, datum=datum, epsg=epsg)

class SurfaceTimeseries(UnstructuredGrid):
    def __init__(self, **kwargs):
        UnstructuredGrid.__init__(self, **kwargs)

    def make_animation(self, **kwargs):
        return _SurfaceTimeseries.make_animation(self, **kwargs)

class SurfaceExtrema(UnstructuredGrid):
    def __init__(self, **kwargs):
        UnstructuredGrid.__init__(self, **kwargs)

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
