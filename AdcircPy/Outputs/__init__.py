from AdcircPy.Mesh import AdcircMesh
from AdcircPy.core.UnstructuredGrid import UnstructuredGrid
from AdcircPy.Outputs import _Outputs
from AdcircPy.Outputs import _Maxele
from AdcircPy.Outputs import _ElevationStations
from AdcircPy.Outputs import _HarmonicConstituentsSurfaceTimeseries

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

class StationTimeseries(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)
    self.Dataset = kwargs.pop('Dataset', None)

class ElevationStationTimeseries(StationTimeseries):
  def __init__(self, x, y, zeta, time, station_name, Dataset=None):
    StationTimeseries.__init__(self, x, y, time, station_name, Dataset)
    self.zeta         = zeta

  @staticmethod
  def from_netcdf(path):
    return _ElevationStations._from_netcdf(path)

  def make_plot(self, station, **kwargs):
    return _ElevationStations._make_plot(self, station, **kwargs)

class HarmonicConstituentsStationTimeseries(StationTimeseries):
    def __init__(self, data_dict):
        for constituent in data_dict.keys():
            setattr(self, constituent, data_dict[constituent])

    @staticmethod
    def from_ascii(path, fort14=None, fort15=None, datum=None, epsg=None):
        return _HarmonicConstituents._from_ascii(path, fort14=fort14, fort15=fort15, datum=datum, epsg=epsg)

class ScalarOutputSurface(UnstructuredGrid):
  def __init__(self, x, y, values, elements, **kwargs):
    UnstructuredGrid.__init__(self, x, y, values, elements, **kwargs)
    self.Dataset = kwargs.pop('Dataset', None)

class Maxele(ScalarOutputSurface):
  def __init__(self, x, y, values, elements, **kwargs):
    ScalarOutputSurface.__init__(self, x, y, values, elements, **kwargs)

  @staticmethod
  def from_ascii(path, fort14):
    _Maxele._from_ascii(path, fort14)
    
  @staticmethod
  def from_netcdf(path, fort14=None):
    return _Maxele._from_netcdf(path, fort14)
    
class OutputSurfaceTimeseries(ScalarOutputSurface):
  def __init__(self, **kwargs):
    ScalarOutputSurface.__init__(self, **kwargs)
    self.values = kwargs.pop('values')

  def make_animation(self, **kwargs):
    return _SurfaceTimeseries.make_animation(self, **kwargs)

class HarmonicConstituentsSurfaceTimeseries(OutputSurfaceTimeseries):
  def __init__(self, **kwargs):
    OutputSurfaceTimeseries.__init__(self, **kwargs)

  @staticmethod
  def from_ascii(path, fort14=None, datum=None, epsg=None):
    return _HarmonicConstituentsSurfaceTimeseries._from_ascii(path, fort14, datum, epsg)
