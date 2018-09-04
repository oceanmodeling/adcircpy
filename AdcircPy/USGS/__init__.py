from AdcircPy.USGS import _HighWaterMarks

class USGS(dict):
  def __init__(self, **kwargs):
    self.epsg   = kwargs.pop("epsg", None)
    dict.__init__(self, **kwargs)

class HighWaterMarks(USGS):
  def __init__(self, **kwargs):
    USGS.__init__(self, **kwargs)
    
  @staticmethod
  def from_event_name(eventName):
    return _HighWaterMarks.from_event_name(eventName)

  @staticmethod
  def from_csv(path):
    return _HighWaterMarks.from_csv(path)

  @staticmethod
  def get_event_list():
    return _HighWaterMarks.get_event_list()


  def filter(self, **kwargs):
    return _HighWaterMarks.filter(self, **kwargs)

  def clip_from_shapefile(self, path, **kwargs):
    return _HighWaterMarks.clip_from_shapefile(self, path, **kwargs)

  def get_xy(self):
    return _HighWaterMarks.get_xy(self)

  def get_values(self):
    return _HighWaterMarks.get_values(self)
    
  def get_xyz(self):
    return _HighWaterMarks.get_xyz(self)
    
  def get_environments(self):
    return _HighWaterMarks.get_environments(self)
        
  def get_counties(self):
    return _HighWaterMarks.get_counties(self)

  def get_from_extent(self, extent, epsg):
    return _HighWaterMarks.get_from_extent(self, extent, epsg)
        
  def export_shapefile(self, path):
    _HighWaterMarks.export_shapefile(self, path)

  def export_to_PostGIS(self, dbname, **kwargs):
    _HighWaterMarks.export_to_PostGIS(self, dbname, **kwargs)

  def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None,**kwargs):
    return _HighWaterMarks.make_plot(self, axes, vmin, vmax, extent, epsg, **kwargs)

  def _parse_args(self):
    _HighWaterMarks._parse_args(self)

  @staticmethod
  def _get_event_id(eventName):
    return _HighWaterMarks._get_event_id(eventName)
