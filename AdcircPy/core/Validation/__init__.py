from AdcircPy.core.Validation import _USGSHighWaterMarks

class USGSHighWaterMarks(dict):
  def __init__(self, **kwargs):
    self.epsg   = kwargs.pop("epsg", None)
    dict.__init__(self, **kwargs)
    
  @staticmethod
  def from_event_name(eventName):
    return _USGSHighWaterMarks.from_event_name(eventName)

  @staticmethod
  def from_csv(path):
    return _USGSHighWaterMarks.from_csv(path)

  @staticmethod
  def get_event_list():
    return _USGSHighWaterMarks.get_event_list()

  @staticmethod
  def _get_event_id(eventName):
    return _USGSHighWaterMarks._get_event_id(eventName)

  def filter(self, excellent=False, good=False, fair=False, poor=False, riverine=False, non_still_water=False, return_count=False, copy=True):
    return _USGSHighWaterMarks.filter(self, excellent, good, fair, poor, riverine, non_still_water, return_count, copy)

  def clip_from_shapefile(self, path, **kwargs):
    return _USGSHighWaterMarks.clip_from_shapefile(self, path, **kwargs)

  def get_xy(self):
    return _USGSHighWaterMarks.get_xy(self)

  def get_values(self):
    return _USGSHighWaterMarks.get_values(self)
    
  def get_xyz(self):
    return _USGSHighWaterMarks.get_xyz(self)
    
  def get_environments(self):
    return _USGSHighWaterMarks.get_environments(self)
        
  def get_counties(self):
    return _USGSHighWaterMarks.get_counties(self)

  def get_from_extent(self, extent, epsg):
    return _USGSHighWaterMarks.get_from_extent(self, extent, epsg)
        
  def export_shapefile(self, path):
    _USGSHighWaterMarks.export_shapefile(self, path)

  def export_to_PostGIS(self, dbname, **kwargs):
    _USGSHighWaterMarks.export_to_PostGIS(self, dbname, **kwargs)

  def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None,**kwargs):
    return _USGSHighWaterMarks.make_plot(self, axes, vmin, vmax, extent, epsg, **kwargs)

  def _parse_args(self):
    _USGSHighWaterMarks._parse_args(self)
